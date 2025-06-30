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

#define THREADS_PER_BLOCK 512
#define MAX_POLY_LEN 262144

void fRead(int *a, int *b, int *n, int *p, int input_id){
    std::string strin = "D:\\study\\concurrent processing\\NTT\\ntt\\nttdata\\" + std::to_string(input_id) + ".in";
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
    std::string strout = "D:\\study\\concurrent processing\\NTT\\ntt\\nttdata\\" + std::to_string(input_id) + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++){
        int x;
        fin>>x;
        if(x != ab[i]){
            // std::cout<<"多项式乘法结果错误"<<std::endl;
            std::cout << "wrong" << std::endl;
            return;
        }
    }
    // std::cout<<"多项式乘法结果正确"<<std::endl;
    std::cout << "correct" << std::endl;
    return;
}

void fWrite(int *ab, int n, int input_id){
    std::string strout = "D:\\study\\concurrent processing\\NTT\\ntt\\files\\" + std::to_string(input_id) + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}

void poly_multiply(int *a, int *b, int *ab, int n, int p){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            ab[i+j]=(1LL * a[i] * b[j] % p + ab[i+j]) % p;
        }
    }
}

uint64_t power(int base, int exp, int mod) 
{
    uint64_t res = 1;
    while (exp > 0) 
    {
        if (exp & 1) 
        {
            res = 1LL * res * base % mod;
        }
        base = 1LL * base * base % mod;
        exp >>= 1;
    }
    return res;
}

__device__ inline int mod_mul(int a, int b, int p) {
    return (int)((1LL * a * b) % p);
}

__device__ inline int mod_add(int a, int b, int p) {
    int res = a + b;
    return (res >= p) ? res - p : res;
}

__device__ inline int mod_sub(int a, int b, int p) {
    int res = a - b;
    return (res < 0) ? res + p : res;
}

// CUDA kernel for one NTT stage (butterfly operation)
__global__ void ntt_stage(int *a, int len, int step, int *wns, int p) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int m = step << 1;
    if (idx >= len / 2) return;

    int j = idx % step;
    int i = (idx / step) * m + j;

    int u = a[i];
    int v = mod_mul(a[i + step], wns[j], p);
    a[i] = mod_add(u, v, p);
    a[i + step] = mod_sub(u, v, p);
}

// Bit reverse permutation (CPU side)
void bit_reverse(std::vector<int> &a, int n) {
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j >= bit; bit >>= 1)
            j -= bit;
        j += bit;
        if (i < j) std::swap(a[i], a[j]);
    }
}

// Modular exponentiation (CPU side)
int power(int base, int exp, int mod) {
    int64_t res = 1;
    while (exp) {
        if (exp & 1) res = res * base % mod;
        base = 1LL * base * base % mod;
        exp >>= 1;
    }
    return res;
}

// Host NTT logic
void NTT_cuda(std::vector<int> &a_host, int n, int p, int invert) {
    bit_reverse(a_host, n);

    int *a_dev, *wns_dev;
    cudaMalloc(&a_dev, sizeof(int) * n);
    cudaMemcpy(a_dev, a_host.data(), sizeof(int) * n, cudaMemcpyHostToDevice);

    for (int len = 2; len <= n; len <<= 1) {
        int wn = power(3, (p - 1) / len, p);
        if (invert == -1)
            wn = power(wn, p - 2, p);

        std::vector<int> wns(len / 2, 1);
        for (int i = 1; i < len / 2; ++i)
            wns[i] = (int)(1LL * wns[i - 1] * wn % p);

        cudaMalloc(&wns_dev, sizeof(int) * (len / 2));
        cudaMemcpy(wns_dev, wns.data(), sizeof(int) * (len / 2), cudaMemcpyHostToDevice);

        int numThreads = len / 2;
        int blocks = (numThreads + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        ntt_stage<<<blocks, THREADS_PER_BLOCK>>>(a_dev, n, len / 2, wns_dev, p);
        cudaDeviceSynchronize();

        cudaFree(wns_dev);
    }

    if (invert == -1) {
        int inv_n = power(n, p - 2, p);
        std::vector<int> a_tmp(n);
        cudaMemcpy(a_tmp.data(), a_dev, sizeof(int) * n, cudaMemcpyDeviceToHost);
        for (int &x : a_tmp) x = (int)(1LL * x * inv_n % p);
        a_host = a_tmp;
    } else {
        cudaMemcpy(a_host.data(), a_dev, sizeof(int) * n, cudaMemcpyDeviceToHost);
    }

    cudaFree(a_dev);
}

int a[300000], b[300000], ab[300000];
int main(int argc, char *argv[]) 
{
    int test_begin = 0;
    int test_end = 4;
    for(int i = test_begin; i <= test_end; ++i) 
    {
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        
        int len = 1;
        while(len < 2 * n_) len <<= 1;
        
        int epochs = 50;
        for (int i = 0; i < epochs; ++i)
        {
            std::vector<int> a_1(len, 0), b_1(len, 0);
            for(int i = 0; i < n_; ++i) 
            {
                a_1[i] = a[i];
                b_1[i] = b[i];
            }
            
            std::vector<int> c(len);
            auto Start = std::chrono::high_resolution_clock::now();

            NTT_cuda(a_1, len, p_, 1);
            NTT_cuda(b_1, len, p_, 1);

            for(int i = 0; i < len; ++i) 
            {
                c[i] = 1LL * a_1[i] * b_1[i] % p_;
            }

            NTT_cuda(c, len, p_, -1);

            auto End = std::chrono::high_resolution_clock::now();

            for(int i = 0; i < 2 * n_ - 1; ++i) 
            {
                ab[i] = c[i];
            }
            std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
            ans += elapsed.count();
        }

        fCheck(ab, n_, i);
        std::cout << "average latency for n = " << n_ << " p = " << p_ << " : " << double(ans / epochs) << " (us) " << std::endl;
        fWrite(ab, n_, i);
    }
    return 0;
}