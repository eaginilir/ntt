#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// CUDA错误检查宏
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " << cudaGetErrorString(error) << std::endl; \
            exit(1); \
        } \
    } while(0)

void fRead(int *a, int *b, int *n, int *p, int input_id){
    std::string strin = "D:\\study\\concurrent processing\\NTT\\ntt\\nttdata\\" + std::to_string(input_id) + ".in";
    std::ifstream fin;
    fin.open(strin, std::ios::in);
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++){
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++){   
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id){
    std::string strout = "D:\\study\\concurrent processing\\NTT\\ntt\\nttdata\\" + std::to_string(input_id) + ".out";
    std::ifstream fin;
    fin.open(strout, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++){
        int x;
        fin>>x;
        if(x != ab[i]){
            std::cout<<"wrong"<<std::endl;
            return;
        }
    }
    std::cout<<"correct"<<std::endl;
    return;
}

void fWrite(int *ab, int n, int input_id){
    std::string strout = "D:\\study\\concurrent processing\\NTT\\ntt\\files\\" + std::to_string(input_id) + ".out";
    std::ofstream fout;
    fout.open(strout, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}

// GPU快速幂函数
__device__ uint64_t power_gpu(int base, int exp, int mod) {
    uint64_t res = 1;
    while (exp > 0) {
        if (exp & 1) {
            res = (1LL * res * base) % mod;
        }
        base = (1LL * base * base) % mod;
        exp >>= 1;
    }
    return res;
}

// CPU快速幂函数
uint64_t power(int base, int exp, int mod) {
    uint64_t res = 1;
    while (exp > 0) {
        if (exp & 1) {
            res = (1LL * res * base) % mod;
        }
        base = (1LL * base * base) % mod;
        exp >>= 1;
    }
    return res;
}

// 预处理旋转因子的核函数
__global__ void precompute_twiddle_factors(int *twiddle_factors, int n, int p, int inv) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int len = 2;
    int offset = 0;
    
    // 计算当前线程对应的len和索引
    int temp_idx = idx;
    while (len <= n && temp_idx >= 0) {
        int group_size = n / len;
        if (temp_idx < group_size) {
            int wn = power_gpu(3, (p - 1) / len, p);
            if (inv == -1) {
                wn = power_gpu(wn, p - 2, p);
            }
            twiddle_factors[offset + temp_idx] = wn;
            return;
        }
        temp_idx -= group_size;
        offset += group_size;
        len <<= 1;
    }
}

// GPU位反转核函数
__global__ void bit_reverse_gpu(int *a, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    int j = 0;
    int temp_idx = idx;
    for (int bit = n >> 1; bit > 0; bit >>= 1) {
        j = (j << 1) | (temp_idx & 1);
        temp_idx >>= 1;
    }
    
    if (idx < j) {
        int temp = a[idx];
        a[idx] = a[j];
        a[j] = temp;
    }
}

// GPU NTT核函数
__global__ void ntt_kernel(int *a, int *twiddle_factors, int n, int len, int p, int inv) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int group_id = idx / (len / 2);
    int pos_in_group = idx % (len / 2);
    int group_start = group_id * len;
    
    if (group_start >= n) return;
    
    // 计算旋转因子索引
    int twiddle_offset = 0;
    int temp_len = 2;
    while (temp_len < len) {
        twiddle_offset += n / temp_len;
        temp_len <<= 1;
    }
    
    int wn = twiddle_factors[twiddle_offset + (n / len) - 1];
    int w = power_gpu(wn, pos_in_group, p);
    
    int u = a[group_start + pos_in_group];
    int v = (1LL * w * a[group_start + pos_in_group + len/2]) % p;
    
    a[group_start + pos_in_group] = (u + v) % p;
    a[group_start + pos_in_group + len/2] = (u - v + p) % p;
}

// GPU逆变换后的除法核函数
__global__ void inverse_divide(int *a, int n, int p) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    int inv_n = power_gpu(n, p - 2, p);
    a[idx] = (1LL * a[idx] * inv_n) % p;
}

// 点乘核函数
__global__ void pointwise_multiply(int *a, int *b, int *c, int n, int p) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    c[idx] = (1LL * a[idx] * b[idx]) % p;
}

// CUDA NTT函数
void NTT_cuda(std::vector<int> &a, int n, int p, int inv) {
    // 获取GPU信息
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    
    int max_threads = prop.maxThreadsPerBlock;
    int sm_count = prop.multiProcessorCount;
    
    // 设置块大小和网格大小
    int block_size = min(max_threads, 1024);
    int grid_size = (n + block_size - 1) / block_size;
    
    // 分配GPU内存
    int *d_a, *d_twiddle;
    CUDA_CHECK(cudaMalloc(&d_a, n * sizeof(int)));
    
    // 计算旋转因子总数
    int twiddle_count = 0;
    for (int len = 2; len <= n; len <<= 1) {
        twiddle_count += n / len;
    }
    CUDA_CHECK(cudaMalloc(&d_twiddle, twiddle_count * sizeof(int)));
    
    // 复制数据到GPU
    CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    
    // 预处理旋转因子
    int twiddle_grid = (twiddle_count + block_size - 1) / block_size;
    precompute_twiddle_factors<<<twiddle_grid, block_size>>>(d_twiddle, n, p, inv);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // 位反转
    bit_reverse_gpu<<<grid_size, block_size>>>(d_a, n);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // NTT主循环
    for (int len = 2; len <= n; len <<= 1) {
        int total_operations = n / 2;
        int ntt_grid = (total_operations + block_size - 1) / block_size;
        
        ntt_kernel<<<ntt_grid, block_size>>>(d_a, d_twiddle, n, len, p, inv);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    // 如果是逆变换，需要除以n
    if (inv == -1) {
        inverse_divide<<<grid_size, block_size>>>(d_a, n, p);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    // 复制结果回CPU
    CUDA_CHECK(cudaMemcpy(a.data(), d_a, n * sizeof(int), cudaMemcpyDeviceToHost));
    
    // 释放GPU内存
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_twiddle));
}

// 多项式乘法CUDA版本
void polynomial_multiply_cuda(std::vector<int> &a, std::vector<int> &b, std::vector<int> &c, int n, int p) {
    // 获取GPU信息
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    
    int block_size = min(prop.maxThreadsPerBlock, 1024);
    int grid_size = (n + block_size - 1) / block_size;
    
    // 分配GPU内存
    int *d_a, *d_b, *d_c;
    CUDA_CHECK(cudaMalloc(&d_a, n * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_b, n * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_c, n * sizeof(int)));
    
    // 复制数据到GPU
    CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    
    // 执行NTT
    NTT_cuda(a, n, p, 1);
    NTT_cuda(b, n, p, 1);
    
    // 复制变换后的数据到GPU
    CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    
    // 点乘
    pointwise_multiply<<<grid_size, block_size>>>(d_a, d_b, d_c, n, p);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // 复制结果回CPU
    CUDA_CHECK(cudaMemcpy(c.data(), d_c, n * sizeof(int), cudaMemcpyDeviceToHost));
    
    // 逆NTT
    NTT_cuda(c, n, p, -1);
    
    // 释放GPU内存
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_c));
}

// CPU版本的NTT (保留原有实现)
void bit_reverse(std::vector<int> &a, int n) {
    for(int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for(; j >= bit; bit >>= 1) {
            j -= bit;
        }
        j += bit;
        if(i < j) {
            std::swap(a[i], a[j]);
        }
    }
}

void NTT_iterative(std::vector<int> &a, int n, int p, int inv) {
    bit_reverse(a, n);

    for(int len = 2; len <= n; len <<= 1) {
        int wn = power(3, (p - 1) / len, p);
        if(inv == -1) {
            wn = power(wn, p - 2, p);
        }

        for(int i = 0; i < n; i += len) {
            int w = 1;
            for(int j = 0; j < len/2; ++j) {
                int u = a[i + j];
                int v = 1LL * w * a[i + j + len/2] % p;
                a[i + j] = (u + v) % p;
                a[i + j + len/2] = (u - v + p) % p;
                w = 1LL * w * wn % p;
            }
        }
    }

    if(inv == -1) {
        int inv_n = power(n, p - 2, p);
        for(int &x : a) {
            x = 1LL * x * inv_n % p;
        }
    }
}

int a[300000], b[300000], ab[300000];

int main(int argc, char *argv[]) {
    // 初始化CUDA
    CUDA_CHECK(cudaSetDevice(0));
    
    // 获取并打印GPU信息
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    std::cout << "GPU: " << prop.name << std::endl;
    std::cout << "SM数量: " << prop.multiProcessorCount << std::endl;
    std::cout << "最大线程数/块: " << prop.maxThreadsPerBlock << std::endl;
    std::cout << "最大块数/网格: " << prop.maxGridSize[0] << std::endl;
    
    int test_begin = 0;
    int test_end = 4;
    
    for(int i = test_begin; i <= test_end; ++i) {
        long double ans_cpu = 0, ans_gpu = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        
        int len = 1;
        while(len < 2 * n_) len <<= 1;
        
        int epochs = 50;
        
        std::cout << "\n测试用例 " << i << ": n = " << n_ << ", p = " << p_ << std::endl;
        
        // CPU版本测试
        for (int epoch = 0; epoch < epochs; ++epoch) {
            std::vector<int> a_1(len, 0), b_1(len, 0);
            for(int j = 0; j < n_; ++j) {
                a_1[j] = a[j];
                b_1[j] = b[j];
            }
            
            std::vector<int> c(len);
            auto Start = std::chrono::high_resolution_clock::now();

            NTT_iterative(a_1, len, p_, 1);
            NTT_iterative(b_1, len, p_, 1);

            for(int j = 0; j < len; ++j) {
                c[j] = 1LL * a_1[j] * b_1[j] % p_;
            }

            NTT_iterative(c, len, p_, -1);

            auto End = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
            ans_cpu += elapsed.count();
        }
        
        // GPU版本测试
        for (int epoch = 0; epoch < epochs; ++epoch) {
            std::vector<int> a_1(len, 0), b_1(len, 0);
            for(int j = 0; j < n_; ++j) {
                a_1[j] = a[j];
                b_1[j] = b[j];
            }
            
            std::vector<int> c(len);
            auto Start = std::chrono::high_resolution_clock::now();

            polynomial_multiply_cuda(a_1, b_1, c, len, p_);

            auto End = std::chrono::high_resolution_clock::now();

            for(int j = 0; j < 2 * n_ - 1; ++j) {
                ab[j] = c[j];
            }
            
            std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
            ans_gpu += elapsed.count();
        }

        fCheck(ab, n_, i);
        std::cout << "CPU平均延迟: " << double(ans_cpu / epochs) << " us" << std::endl;
        std::cout << "GPU平均延迟: " << double(ans_gpu / epochs) << " us" << std::endl;
        std::cout << "加速比: " << double(ans_cpu / ans_gpu) << "x" << std::endl;
        
        fWrite(ab, n_, i);
    }
    
    return 0;
}