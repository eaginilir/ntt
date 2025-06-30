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

// CPU快速幂
__host__ __device__ uint64_t power(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) {
            res = ((__uint128_t)res * base) % mod;
        }
        base = ((__uint128_t)base * base) % mod;
        exp >>= 1;
    }
    return res;
}

// Barrett结构体
struct Barrett {
    uint64_t mod;
    uint64_t factor;
    __host__ __device__ Barrett() {}
    __host__ __device__ Barrett(uint64_t m) : mod(m) {
        factor = ~0ULL / m;
    }
    __host__ __device__ inline uint64_t reduce(uint64_t x) const {
        uint64_t q = ((uint64_t)(((__uint128_t)x * factor) >> 64));
        uint64_t r = x - q * mod;
        return r >= mod ? r - mod : r;
    }
    __host__ __device__ inline uint64_t mul(uint64_t a, uint64_t b) const {
        return reduce((__uint128_t)a * b);
    }
    __host__ __device__ inline uint64_t add(uint64_t a, uint64_t b) const {
        uint64_t r = a + b;
        return r >= mod ? r - mod : r;
    }
    __host__ __device__ inline uint64_t sub(uint64_t a, uint64_t b) const {
        return a >= b ? a - b : mod - (b - a);
    }
};

// GPU内核 - 位反转
__global__ void bit_reverse_kernel(uint64_t *a, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    int j = 0, t = idx;
    for (int bit = n >> 1; bit > 0; bit >>= 1) {
        j = (j << 1) | (t & 1);
        t >>= 1;
    }
    if (idx < j) {
        uint64_t tmp = a[idx];
        a[idx] = a[j];
        a[j] = tmp;
    }
}

// GPU内核 - NTT主循环
__global__ void ntt_kernel(uint64_t *a, int n, int len, const uint64_t *twiddles, Barrett barrett, int twiddle_offset) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_butterflies = n / len;
    if (idx >= total_butterflies * (len / 2)) return;
    int butterfly_group = idx / (len / 2);
    int j = idx % (len / 2);
    int i = butterfly_group * len;
    uint64_t w = twiddles[twiddle_offset + j];
    uint64_t u = a[i + j];
    uint64_t v = barrett.mul(w, a[i + j + len / 2]);
    a[i + j] = barrett.add(u, v);
    a[i + j + len / 2] = barrett.sub(u, v);
}

// GPU内核 - 逆NTT的最终除法
__global__ void inv_ntt_final_kernel(uint64_t *a, int n, uint64_t inv_n, Barrett barrett) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    a[idx] = barrett.mul(a[idx], inv_n);
}

// GPU内核 - 向量乘法
__global__ void pointwise_mul_kernel(uint64_t *c, const uint64_t *a, const uint64_t *b, int n, Barrett barrett) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    c[idx] = barrett.mul(a[idx], b[idx]);
}

class CudaBarrettNTT {
private:
    Barrett barrett;
    uint64_t *d_twiddles;
    int max_n;
    uint64_t p;
public:
    CudaBarrettNTT(uint64_t p_, int max_n_) : p(p_), max_n(max_n_) {
        barrett = Barrett(p_);
        d_twiddles = nullptr;
        precompute_twiddles();
    }
    ~CudaBarrettNTT() {
        if (d_twiddles) {
            CUDA_CHECK(cudaFree(d_twiddles));
        }
    }
private:
    void precompute_twiddles() {
        int g = 3;
        std::vector<uint64_t> twiddles;
        // 正向
        for (int len = 2; len <= max_n; len <<= 1) {
            uint64_t wn = power(g, (p - 1) / len, p);
            uint64_t w = 1;
            for (int j = 0; j < len / 2; ++j) {
                twiddles.push_back(w);
                w = (__uint128_t(w) * wn) % p;
            }
        }
        // 逆向
        for (int len = 2; len <= max_n; len <<= 1) {
            uint64_t wn = power(g, (p - 1) / len, p);
            wn = power(wn, p - 2, p);
            uint64_t w = 1;
            for (int j = 0; j < len / 2; ++j) {
                twiddles.push_back(w);
                w = (__uint128_t(w) * wn) % p;
            }
        }
        CUDA_CHECK(cudaMalloc(&d_twiddles, twiddles.size() * sizeof(uint64_t)));
        CUDA_CHECK(cudaMemcpy(d_twiddles, twiddles.data(), twiddles.size() * sizeof(uint64_t), cudaMemcpyHostToDevice));
    }
    int get_twiddle_offset(int len, bool inverse) {
        int offset = 0;
        for (int l = 2; l < len; l <<= 1) {
            offset += l / 2;
        }
        if (inverse) {
            int forward_total = 0;
            for (int l = 2; l <= max_n; l <<= 1) {
                forward_total += l / 2;
            }
            offset += forward_total;
        }
        return offset;
    }
public:
    void ntt_gpu(uint64_t *d_a, int n, bool inverse) {
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
        int blockSize = std::min(1024, prop.maxThreadsPerBlock);
        int gridSize;

        // 位反转
        gridSize = (n + blockSize - 1) / blockSize;
        bit_reverse_kernel<<<gridSize, blockSize>>>(d_a, n);
        CUDA_CHECK(cudaDeviceSynchronize());

        // NTT主循环
        for (int len = 2; len <= n; len <<= 1) {
            int total_ops = n / len * (len / 2);
            gridSize = (total_ops + blockSize - 1) / blockSize;
            int twiddle_offset = get_twiddle_offset(len, inverse);
            ntt_kernel<<<gridSize, blockSize>>>(d_a, n, len, d_twiddles, barrett, twiddle_offset);
            CUDA_CHECK(cudaDeviceSynchronize());
        }

        // 逆NTT的最终除法
        if (inverse) {
            uint64_t inv_n = power(n, p - 2, p);
            gridSize = (n + blockSize - 1) / blockSize;
            inv_ntt_final_kernel<<<gridSize, blockSize>>>(d_a, n, inv_n, barrett);
            CUDA_CHECK(cudaDeviceSynchronize());
        }
    }

    void polynomial_multiply(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b,
                           std::vector<uint64_t> &result, int n) {
        int len = 1;
        while (len < 2 * n) len <<= 1;
        result.resize(len);
        uint64_t *d_a, *d_b, *d_c;
        CUDA_CHECK(cudaMalloc(&d_a, len * sizeof(uint64_t)));
        CUDA_CHECK(cudaMalloc(&d_b, len * sizeof(uint64_t)));
        CUDA_CHECK(cudaMalloc(&d_c, len * sizeof(uint64_t)));
        CUDA_CHECK(cudaMemset(d_a, 0, len * sizeof(uint64_t)));
        CUDA_CHECK(cudaMemset(d_b, 0, len * sizeof(uint64_t)));
        CUDA_CHECK(cudaMemcpy(d_a, a.data(), n * sizeof(uint64_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_b, b.data(), n * sizeof(uint64_t), cudaMemcpyHostToDevice));
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
        int blockSize = std::min(1024, prop.maxThreadsPerBlock);
        int gridSize = (len + blockSize - 1) / blockSize;

        ntt_gpu(d_a, len, false);
        ntt_gpu(d_b, len, false);
        pointwise_mul_kernel<<<gridSize, blockSize>>>(d_c, d_a, d_b, len, barrett);
        CUDA_CHECK(cudaDeviceSynchronize());
        ntt_gpu(d_c, len, true);
        CUDA_CHECK(cudaMemcpy(result.data(), d_c, len * sizeof(uint64_t), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaFree(d_a));
        CUDA_CHECK(cudaFree(d_b));
        CUDA_CHECK(cudaFree(d_c));
    }
};

int a[300000], b[300000], ab[300000];

int main(int argc, char *argv[]) {
    int deviceCount;
    CUDA_CHECK(cudaGetDeviceCount(&deviceCount));
    if (deviceCount == 0) {
        std::cerr << "No CUDA devices found!" << std::endl;
        return 1;
    }
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    std::cout << "Using GPU: " << prop.name << std::endl;
    std::cout << "SM count: " << prop.multiProcessorCount << std::endl;
    std::cout << "Max threads per block: " << prop.maxThreadsPerBlock << std::endl;
    std::cout << "Max threads per SM: " << prop.maxThreadsPerMultiProcessor << std::endl;

    int test_begin = 0;
    int test_end = 4;
    for(int i = test_begin; i <= test_end; ++i) {
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        std::cout << "Processing test case " << i << " with n=" << n_ << ", p=" << p_ << std::endl;
        int len = 1;
        while(len < 2 * n_) len <<= 1;
        CudaBarrettNTT* cuda_ntt = nullptr;
        try {
            cuda_ntt = new CudaBarrettNTT(p_, len);
        } catch (const std::exception& e) {
            std::cerr << "Failed to create CUDA Barrett NTT for p=" << p_ << ": " << e.what() << std::endl;
            continue;
        }
        int epochs = 50;
        for (int epoch = 0; epoch < epochs; ++epoch) {
            std::vector<uint64_t> a_vec(n_), b_vec(n_), result;
            for (int j = 0; j < n_; ++j) {
                a_vec[j] = a[j];
                b_vec[j] = b[j];
            }
            auto Start = std::chrono::high_resolution_clock::now();
            cuda_ntt->polynomial_multiply(a_vec, b_vec, result, n_);
            auto End = std::chrono::high_resolution_clock::now();
            for (int j = 0; j < 2 * n_ - 1; ++j) {
                ab[j] = result[j] % p_;
            }
            std::chrono::duration<double, std::ratio<1,1000>> elapsed = End - Start;
            ans += elapsed.count();
        }
        fCheck(ab, n_, i);
        std::cout << "GPU average latency for n = " << n_ << " p = " << p_ 
                  << " : " << double(ans / epochs) << " (us) " << std::endl;
        fWrite(ab, n_, i);
    }
    return 0;
}