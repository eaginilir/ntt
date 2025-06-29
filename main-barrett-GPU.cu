#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// CUDA错误检查宏
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " << cudaGetErrorString(err) << std::endl; \
            exit(1); \
        } \
    } while(0)

void fRead(int *a, int *b, int *n, int *p, int input_id){
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++){
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++){   
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id){
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++){
        int x;
        fin>>x;
        if(x != ab[i]){
            std::cout<<"多项式乘法结果错误"<<std::endl;
            return;
        }
    }
    std::cout<<"多项式乘法结果正确"<<std::endl;
    return;
}

void fWrite(int *ab, int n, int input_id){
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}

// CPU端快速幂实现
__host__ uint64_t power(int base, int exp, int mod) {
    uint64_t res = 1;
    while (exp > 0) {
        if (exp & 1) {
            res = 1LL * res * base % mod;
        }
        base = 1LL * base * base % mod;
        exp >>= 1;
    }
    return res;
}

// GPU端快速幂实现
__device__ uint64_t device_power(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t res = 1;
    while (exp > 0) {
        if (exp & 1) {
            res = (__uint128_t)res * base % mod;
        }
        base = (__uint128_t)base * base % mod;
        exp >>= 1;
    }
    return res;
}

// GPU端Barrett约简类
class __device__ Barrett {
public:
    uint64_t mod_;
    int shift_;
    uint64_t factor_;

    __device__ Barrett() {}
    
    __device__ Barrett(uint64_t mod): mod_(mod) {
        shift_ = 64;
        factor_ = (__uint128_t(1) << shift_) / mod_;
    }

    __device__ inline uint64_t reduce(uint64_t x) const {
        return x % mod_;
    }

    __device__ inline uint64_t reduce(__uint128_t x) const {
        uint64_t q = (uint64_t)((x * factor_) >> shift_);
        uint64_t r = (uint64_t)x - q * mod_;
        return r >= mod_ ? r - mod_ : r;
    }

    __device__ inline uint64_t multiply(uint64_t a, uint64_t b) const {
        return reduce((__uint128_t)a * b);
    }

    __device__ inline uint64_t add(uint64_t a, uint64_t b) const {
        uint64_t r = a + b;
        return r >= mod_ ? r - mod_ : r;
    }

    __device__ inline uint64_t sub(uint64_t a, uint64_t b) const {
        return a >= b ? a - b : mod_ - (b - a);
    }
};

// GPU端位反转核函数
__global__ void bit_reverse_kernel(uint64_t *a, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    int j = 0;
    int temp_idx = idx;
    int bits = 0;
    int temp_n = n;
    while (temp_n > 1) {
        bits++;
        temp_n >>= 1;
    }
    
    for (int i = 0; i < bits; i++) {
        j = (j << 1) | (temp_idx & 1);
        temp_idx >>= 1;
    }
    
    if (idx < j) {
        uint64_t temp = a[idx];
        a[idx] = a[j];
        a[j] = temp;
    }
}

// 预处理旋转因子的核函数
__global__ void precompute_twiddles_kernel(uint64_t *twiddles, int max_len, int p, Barrett barrett) {
    int len = blockIdx.x + 1;  // len从1开始到max_len
    if (len > max_len) return;
    
    int idx = threadIdx.x;
    int half_len = len / 2;
    
    if (idx >= half_len || len < 2) return;
    
    int g = 3;
    uint64_t wn = device_power(g, (p - 1) / len, p);
    uint64_t w = device_power(wn, idx, p);
    
    // 正变换的旋转因子
    twiddles[len * max_len + idx] = w;
    
    // 逆变换的旋转因子
    uint64_t inv_wn = device_power(wn, p - 2, p);
    uint64_t inv_w = device_power(inv_wn, idx, p);
    twiddles[len * max_len + half_len + idx] = inv_w;
}

// GPU NTT核函数
__global__ void ntt_kernel(uint64_t *a, uint64_t *twiddles, int n, int len, int inv, int max_len, Barrett barrett) {
    int i = blockIdx.x * len;  // 每个block处理一个长度为len的段
    int j = threadIdx.x;       // 每个thread处理段内的一个蝶形操作
    
    if (i >= n || j >= len/2) return;
    
    int half_len = len / 2;
    int twiddle_offset = len * max_len + (inv == -1 ? half_len : 0);
    
    uint64_t w = twiddles[twiddle_offset + j];
    
    uint64_t u = a[i + j];
    uint64_t v = barrett.multiply(w, a[i + j + half_len]);
    
    a[i + j] = barrett.add(u, v);
    a[i + j + half_len] = barrett.sub(u, v);
}

// 点乘核函数
__global__ void pointwise_multiply_kernel(uint64_t *a, uint64_t *b, uint64_t *c, int n, Barrett barrett) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    c[idx] = barrett.multiply(a[idx], b[idx]);
}

// 逆变换后的除法核函数
__global__ void divide_by_n_kernel(uint64_t *a, int n, int inv_n, Barrett barrett) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    a[idx] = barrett.multiply(a[idx], inv_n);
}

// GPU加速的NTT实现
void NTT_GPU(std::vector<uint64_t> &a, int n, int p, int inv) {
    // 获取GPU设备属性
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    
    int max_threads = prop.maxThreadsPerBlock;  // 通常为1024
    int sm_count = prop.multiProcessorCount;    // SM数量
    
    // 动态调整块大小
    int block_size = min(512, max_threads);  // 每个block的线程数
    int grid_size = (n + block_size - 1) / block_size;  // grid大小
    
    Barrett host_barrett(p);
    
    // 分配GPU内存
    uint64_t *d_a;
    uint64_t *d_twiddles;
    
    CUDA_CHECK(cudaMalloc(&d_a, n * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_twiddles, n * n * sizeof(uint64_t)));  // 为所有长度的旋转因子分配空间
    
    // 复制数据到GPU
    CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(uint64_t), cudaMemcpyHostToDevice));
    
    // 预处理所有旋转因子
    dim3 precomp_grid(n);
    dim3 precomp_block(n/2);
    precompute_twiddles_kernel<<<precomp_grid, precomp_block>>>(d_twiddles, n, p, host_barrett);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // 位反转
    bit_reverse_kernel<<<grid_size, block_size>>>(d_a, n);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // NTT主循环
    for (int len = 2; len <= n; len <<= 1) {
        int num_groups = n / len;
        int threads_per_group = len / 2;
        
        // 动态调整线程配置
        dim3 ntt_grid(num_groups);
        dim3 ntt_block(min(threads_per_group, max_threads));
        
        ntt_kernel<<<ntt_grid, ntt_block>>>(d_a, d_twiddles, n, len, inv, n, host_barrett);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    // 如果是逆变换，除以n
    if (inv == -1) {
        int inv_n = power(n, p - 2, p);
        divide_by_n_kernel<<<grid_size, block_size>>>(d_a, n, inv_n, host_barrett);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    // 复制结果回CPU
    CUDA_CHECK(cudaMemcpy(a.data(), d_a, n * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    
    // 释放GPU内存
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_twiddles));
}

// 完整的多项式乘法GPU实现
void poly_multiply_GPU(std::vector<uint64_t> &a, std::vector<uint64_t> &b, std::vector<uint64_t> &c, int len, int p) {
    // 获取GPU设备属性
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    
    int block_size = min(512, prop.maxThreadsPerBlock);
    int grid_size = (len + block_size - 1) / block_size;
    
    Barrett host_barrett(p);
    
    // 分配GPU内存
    uint64_t *d_a, *d_b, *d_c, *d_twiddles;
    
    CUDA_CHECK(cudaMalloc(&d_a, len * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_b, len * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_c, len * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_twiddles, len * len * sizeof(uint64_t)));
    
    // 复制数据到GPU
    CUDA_CHECK(cudaMemcpy(d_a, a.data(), len * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b.data(), len * sizeof(uint64_t), cudaMemcpyHostToDevice));
    
    // 预处理旋转因子
    dim3 precomp_grid(len);
    dim3 precomp_block(len/2);
    precompute_twiddles_kernel<<<precomp_grid, precomp_block>>>(d_twiddles, len, p, host_barrett);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // 前向NTT
    // a的NTT
    bit_reverse_kernel<<<grid_size, block_size>>>(d_a, len);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    for (int length = 2; length <= len; length <<= 1) {
        int num_groups = len / length;
        int threads_per_group = length / 2;
        
        dim3 ntt_grid(num_groups);
        dim3 ntt_block(min(threads_per_group, prop.maxThreadsPerBlock));
        
        ntt_kernel<<<ntt_grid, ntt_block>>>(d_a, d_twiddles, len, length, 1, len, host_barrett);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    // b的NTT
    bit_reverse_kernel<<<grid_size, block_size>>>(d_b, len);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    for (int length = 2; length <= len; length <<= 1) {
        int num_groups = len / length;
        int threads_per_group = length / 2;
        
        dim3 ntt_grid(num_groups);
        dim3 ntt_block(min(threads_per_group, prop.maxThreadsPerBlock));
        
        ntt_kernel<<<ntt_grid, ntt_block>>>(d_b, d_twiddles, len, length, 1, len, host_barrett);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    // 点乘
    pointwise_multiply_kernel<<<grid_size, block_size>>>(d_a, d_b, d_c, len, host_barrett);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // 逆NTT
    bit_reverse_kernel<<<grid_size, block_size>>>(d_c, len);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    for (int length = 2; length <= len; length <<= 1) {
        int num_groups = len / length;
        int threads_per_group = length / 2;
        
        dim3 ntt_grid(num_groups);
        dim3 ntt_block(min(threads_per_group, prop.maxThreadsPerBlock));
        
        ntt_kernel<<<ntt_grid, ntt_block>>>(d_c, d_twiddles, len, length, -1, len, host_barrett);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    // 除以n
    int inv_n = power(len, p - 2, p);
    divide_by_n_kernel<<<grid_size, block_size>>>(d_c, len, inv_n, host_barrett);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // 复制结果回CPU
    CUDA_CHECK(cudaMemcpy(c.data(), d_c, len * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    
    // 释放GPU内存
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_c));
    CUDA_CHECK(cudaFree(d_twiddles));
}

int a[300000], b[300000], ab[300000];

int main(int argc, char *argv[]) {
    // 初始化CUDA
    CUDA_CHECK(cudaSetDevice(0));
    
    // 打印GPU信息
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    std::cout << "GPU: " << prop.name << std::endl;
    std::cout << "SM Count: " << prop.multiProcessorCount << std::endl;
    std::cout << "Max Threads per Block: " << prop.maxThreadsPerBlock << std::endl;
    std::cout << "Max Threads per SM: " << prop.maxThreadsPerMultiProcessor << std::endl;

    int test_begin = 0;
    int test_end = 4;
    
    for(int i = test_begin; i <= test_end; ++i) {
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        
        int len = 1;
        while(len < 2*n_) {
            len <<= 1;
        }

        int epochs = 50;
        for (int epoch = 0; epoch < epochs; ++epoch) {
            std::vector<uint64_t> a_1(len, 0);
            std::vector<uint64_t> b_1(len, 0);
            std::vector<uint64_t> c(len, 0);
            
            for (int j = 0; j < n_; ++j) {
                a_1[j] = a[j];
                b_1[j] = b[j];
            }
            
            auto Start = std::chrono::high_resolution_clock::now();
            
            // 使用GPU加速的多项式乘法
            poly_multiply_GPU(a_1, b_1, c, len, p_);
            
            auto End = std::chrono::high_resolution_clock::now();
            
            for (int j = 0; j < 2 * n_ - 1; ++j) {
                ab[j] = c[j];
            }
            
            std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
            ans += elapsed.count();
        }
        
        fCheck(ab, n_, i);
        std::cout << "GPU average latency for n = " << n_ << " p = " << p_ << " : " << double(ans / epochs) << " (ms) " << std::endl;
        fWrite(ab, n_, i);
    }
    
    return 0;
}