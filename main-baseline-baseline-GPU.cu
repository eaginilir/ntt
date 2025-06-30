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

void fRead(int *a, int *b, int *n, int *p, int input_id)
{
    std::string strin = "./nttdata/" + std::to_string(input_id) + ".in";
    std::ifstream fin;
    fin.open(strin, std::ios::in);
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++)
    {
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++)
    {   
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id)
{
    std::string strout = "./nttdata/" + std::to_string(input_id) + ".out";
    std::ifstream fin;
    fin.open(strout, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++)
    {
        int x;
        fin>>x;
        if(x != ab[i])
        {
            std::cout<<"多项式乘法结果错误"<<std::endl;
            return;
        }
    }
    std::cout<<"多项式乘法结果正确"<<std::endl;
    return;
}

void fWrite(int *ab, int n, int input_id)
{
    std::string strout = "./files/" + std::to_string(input_id) + ".out";
    std::ofstream fout;
    fout.open(strout, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++)
    {
        fout<<ab[i]<<'\n';
    }
}

// GPU快速幂函数
__device__ uint64_t power_gpu(int base, int exp, int mod) 
{
    uint64_t res = 1;
    while (exp > 0) 
    {
        if (exp & 1) 
        {
            res = (1LL * res * base) % mod;
        }
        base = (1LL * base * base) % mod;
        exp >>= 1;
    }
    return res;
}

// CPU快速幂函数
uint64_t power(int base, int exp, int mod) 
{
    uint64_t res = 1;
    while (exp > 0) 
    {
        if (exp & 1) 
        {
            res = (1LL * res * base) % mod;
        }
        base = (1LL * base * base) % mod;
        exp >>= 1;
    }
    return res;
}

// 最简单的GPU位反转核函数 - 直接对应CPU版本
__global__ void bit_reverse_gpu(int *a, int n) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    
    // 对每个位置计算对应的位反转位置
    int j = 0;
    int temp_i = i;
    for (int bit = n >> 1; bit > 0; bit >>= 1) 
    {
        j = (j << 1) | (temp_i & 1);
        temp_i >>= 1;
    }
    
    // 只有当i < j时才交换，避免重复交换
    if (i < j) 
    {
        int temp = a[i];
        a[i] = a[j];
        a[j] = temp;
    }
}

// 最简单的GPU NTT核函数 - 直接对应CPU的嵌套循环
__global__ void ntt_kernel_simple(int *a, int n, int len, int p, int inv) 
{
    // 每个线程处理一个蝶形运算
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // 计算当前线程对应的i和j
    int total_butterflies = n / 2;
    if (idx >= total_butterflies) return;
    
    // 计算wn (len级别的单位根)
    int wn = power_gpu(3, (p - 1) / len, p);
    if (inv == -1) 
    {
        wn = power_gpu(wn, p - 2, p);
    }
    
    // 找到当前线程对应的段和段内位置
    int segment_size = len;
    int butterflies_per_segment = len / 2;
    int segment_idx = idx / butterflies_per_segment;
    int pos_in_segment = idx % butterflies_per_segment;
    
    int i = segment_idx * segment_size + pos_in_segment;
    int j = i + len / 2;
    
    if (j >= n) return;
    
    // 计算w = wn^pos_in_segment
    int w = power_gpu(wn, pos_in_segment, p);
    
    // 蝶形运算
    int u = a[i];
    int v = (1LL * w * a[j]) % p;
    a[i] = (u + v) % p;
    a[j] = (u - v + p) % p;
}

// 逆变换后的除法核函数
__global__ void inverse_divide(int *a, int n, int p) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    int inv_n = power_gpu(n, p - 2, p);
    a[idx] = (1LL * a[idx] * inv_n) % p;
}

// 点乘核函数
__global__ void pointwise_multiply(int *a, int *b, int *c, int n, int p) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    c[idx] = (1LL * a[idx] * b[idx]) % p;
}

// 最简单的CUDA NTT函数 - 直接对应CPU版本的结构
void NTT_cuda_simple(std::vector<int> &a, int n, int p, int inv) 
{
    // 分配GPU内存
    int *d_a;
    CUDA_CHECK(cudaMalloc(&d_a, n * sizeof(int)));
    
    // 复制数据到GPU
    CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    
    // 设置简单的块大小
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;
    
    // 位反转
    bit_reverse_gpu<<<grid_size, block_size>>>(d_a, n);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // NTT主循环 - 对应CPU版本的for(int len = 2; len <= n; len <<= 1)
    for (int len = 2; len <= n; len <<= 1) 
    {
        int butterflies = n / 2;  // 每一轮的蝶形运算总数
        int butterfly_grid = (butterflies + block_size - 1) / block_size;
        
        ntt_kernel_simple<<<butterfly_grid, block_size>>>(d_a, n, len, p, inv);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    // 如果是逆变换，需要除以n
    if (inv == -1) 
    {
        inverse_divide<<<grid_size, block_size>>>(d_a, n, p);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    
    // 复制结果回CPU
    CUDA_CHECK(cudaMemcpy(a.data(), d_a, n * sizeof(int), cudaMemcpyDeviceToHost));
    
    // 释放GPU内存
    CUDA_CHECK(cudaFree(d_a));
}

// 简单的多项式乘法CUDA版本
void polynomial_multiply_cuda_simple(std::vector<int> &a, std::vector<int> &b, std::vector<int> &c, int n, int p) 
{
    // 执行NTT
    NTT_cuda_simple(a, n, p, 1);
    NTT_cuda_simple(b, n, p, 1);
    
    // 点乘
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;
    
    int *d_a, *d_b, *d_c;
    CUDA_CHECK(cudaMalloc(&d_a, n * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_b, n * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_c, n * sizeof(int)));
    
    CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    
    pointwise_multiply<<<grid_size, block_size>>>(d_a, d_b, d_c, n, p);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    CUDA_CHECK(cudaMemcpy(c.data(), d_c, n * sizeof(int), cudaMemcpyDeviceToHost));
    
    // 逆NTT
    NTT_cuda_simple(c, n, p, -1);
    
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_c));
}

// CPU版本的NTT (保留原有实现)
void bit_reverse(std::vector<int> &a, int n) 
{
    for(int i = 1, j = 0; i < n; ++i) 
    {
        int bit = n >> 1;
        for(; j >= bit; bit >>= 1) 
        {
            j -= bit;
        }
        j += bit;
        if(i < j) 
        {
            std::swap(a[i], a[j]);
        }
    }
}

void NTT_iterative(std::vector<int> &a, int n, int p, int inv) 
{
    bit_reverse(a, n);

    for(int len = 2; len <= n; len <<= 1) 
    {
        int wn = power(3, (p - 1) / len, p);
        if(inv == -1) 
        {
            wn = power(wn, p - 2, p);
        }

        for(int i = 0; i < n; i += len) 
        {
            int w = 1;
            for(int j = 0; j < len/2; ++j) 
            {
                int u = a[i + j];
                int v = 1LL * w * a[i + j + len/2] % p;
                a[i + j] = (u + v) % p;
                a[i + j + len/2] = (u - v + p) % p;
                w = 1LL * w * wn % p;
            }
        }
    }

    if(inv == -1) 
    {
        int inv_n = power(n, p - 2, p);
        for(int &x : a) 
        {
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
        for (int epoch = 0; epoch < epochs; ++epoch) 
        {
            std::vector<int> a_1(len, 0), b_1(len, 0);
            for(int j = 0; j < n_; ++j) 
            {
                a_1[j] = a[j];
                b_1[j] = b[j];
            }
            
            std::vector<int> c(len);
            auto Start = std::chrono::high_resolution_clock::now();

            NTT_iterative(a_1, len, p_, 1);
            NTT_iterative(b_1, len, p_, 1);

            for(int j = 0; j < len; ++j) 
            {
                c[j] = 1LL * a_1[j] * b_1[j] % p_;
            }

            NTT_iterative(c, len, p_, -1);

            auto End = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
            ans_cpu += elapsed.count();
        }
        
        // GPU版本测试 - 使用简化版本
        for (int epoch = 0; epoch < epochs; ++epoch) 
        {
            std::vector<int> a_1(len, 0), b_1(len, 0);
            for(int j = 0; j < n_; ++j) 
            {
                a_1[j] = a[j];
                b_1[j] = b[j];
            }
            
            std::vector<int> c(len);
            auto Start = std::chrono::high_resolution_clock::now();

            polynomial_multiply_cuda_simple(a_1, b_1, c, len, p_);

            auto End = std::chrono::high_resolution_clock::now();

            for(int j = 0; j < 2 * n_ - 1; ++j) 
            {
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