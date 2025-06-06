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
#include <arm_neon.h>
#include <pthread.h>
#include <stdint.h>
#include <shared_mutex>
#include <cassert>
#include <thread>
#include <queue>
#include <functional>
#include <future>
#include <atomic>
#include <condition_variable>
#include <mpi.h>

// 可以自行添加需要的头文件

void fRead(uint64_t *a, uint64_t *b, uint64_t *n, uint64_t *p, int input_id){
    // 数据输入函数
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

void fCheck(uint64_t *ab, int n, int input_id){
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    std::string logout = "files/log.out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    std::ofstream fout(logout, std::ios::app);
    fin.open(data_path, std::ios::in);

    // if(input_id==4)
    // {
    //     bool correct = true;
    //     for (int i = 0; i < n * 2 - 1; i++) {
    //         uint64_t x;
    //         fin >> x;
    //         if (x != ab[i]) {
    //             fout << "错误位置: " << i 
    //                 << ", 期望值: " << x 
    //                 << ", 实际值: " << ab[i] 
    //                 << std::endl;
    //             correct = false;
    //         }
    //     }
    // }

    for (int i = 0; i < n * 2 - 1; i++){
        uint64_t x;
        fin>>x;
        if(x != ab[i]){
            std::cout<<"多项式乘法结果错误"<<std::endl;
            // std::cout << i << std::endl;
            // std::cout << x << std::endl;
            // std::cout << ab[i] << std::endl;
            return;
        }
    }
    std::cout<<"多项式乘法结果正确"<<std::endl;
    return;
}

void fWrite(uint64_t *ab, int n, int input_id){
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
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

std::string to_string(__uint128_t value) {
    if (value == 0) return "0";
    std::string result;
    while (value > 0) {
        result = char('0' + value % 10) + result;
        value /= 10;
    }
    return result;
}

void fWrite(__uint128_t *ab, int n, int input_id){
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<to_string(ab[i])<<'\n';
    }
}

void fWrite(uint64_t *ab, int n, std::string path,int test_id){
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string num = std::to_string(test_id);
    std::string strout = str1 + path + num +".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}

void fWrite(__uint128_t *ab, int n, std::string path,int test_id){
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string num = std::to_string(test_id);
    std::string strout = str1 + path + num + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<to_string(ab[i])<<'\n';
    }
}

void poly_multiply(int *a, int *b, int *ab, int n, int p){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            ab[i+j]=(1LL * a[i] * b[j] % p + ab[i+j]) % p;
        }
    }
}

//快速幂
uint64_t power(uint64_t base, uint64_t exp, uint64_t mod) 
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

//Barrett模乘
class barrett
{
public:
    barrett(uint64_t mod_) : mod(mod_) 
    {
        k = 64;
        mu = (__uint128_t(1) << k) / mod;

        mod_neg = (~mod) + 1;
    }

    inline uint64_t reduce(uint64_t x) const
    {
        if (x < mod) return x;

        __uint128_t q = ((__uint128_t)x * mu) >> k;

        uint64_t r = x - (uint64_t)q * mod;

        return r >= mod ? r - mod : r;
    }

    inline uint64_t mul(uint64_t a, uint64_t b) const
    {
        __uint128_t product = (__uint128_t)a * b;
        return reduce_128(product);
    }

    inline uint64_t reduce_128(__uint128_t x) const
    {
        if (x < mod) return (uint64_t)x;

        __uint128_t q = (x >> (k - 64)) * mu >> 64;

        __uint128_t qm = q * mod;
        uint64_t r = (uint64_t)(x - qm);

        if (r >= mod) r -= mod;
        if (r >= mod) r -= mod;
        
        return r;
    }

    // 模加法
    inline uint64_t add(uint64_t a, uint64_t b) const 
    {
        uint64_t sum = a + b;
        return sum >= mod ? sum - mod : sum;
    }

    // 模减法
    inline uint64_t sub(uint64_t a, uint64_t b) const 
    {
        return a >= b ? a - b : a + mod - b;
    }

    inline uint64_t fast_add(uint64_t a, uint64_t b) const 
    {
        uint64_t sum = a + b;
        uint64_t mask = -(sum >= mod);
        return sum - (mask & mod);
    }

    inline uint64_t fast_sub(uint64_t a, uint64_t b) const 
    {
        uint64_t diff = a - b;
        uint64_t mask = -(a < b);
        return diff + (mask & mod);
    }
public:
    uint64_t mod;
    __uint128_t mu;
    int k;
    uint64_t mod_neg;
};

// 线程池实现
class ThreadPool 
{
public:
    ThreadPool(size_t numThreads) : stop(false), activeThreads(0) 
    {
        for (size_t i = 0; i < numThreads; ++i) 
        {
            workers.emplace_back([this] 
                {
                while (true) 
                {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queueMutex);
                        condition.wait(lock, [this] { return stop || !tasks.empty(); });
                        if (stop && tasks.empty()) 
                        {
                            return;
                        }
                        task = std::move(tasks.front());
                        tasks.pop();
                        activeThreads++; // Increment active thread count
                    }
                    
                    task(); // Execute the task
                    
                    {
                        std::unique_lock<std::mutex> lock(queueMutex);
                        activeThreads--; // Decrement active thread count
                        lock.unlock();
                        if (tasks.empty() && activeThreads == 0)
                        {
                            completionCondition.notify_all(); // Notify when all tasks are done
                        }
                    }
                }
            });
        }
    }

    template<class F>
    void enqueue(F&& f) 
    {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            if (stop) 
            {
                throw std::runtime_error("enqueue on stopped ThreadPool");
            }
            tasks.emplace(std::forward<F>(f));
        }
        condition.notify_one();
    }

    // Improved waitForAll that ensures all tasks are completed
    void waitForAll() 
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        while (true) 
        {
            if (tasks.empty() && activeThreads == 0) break;

            // 立即抢占执行
            if (!tasks.empty()) 
            {
                auto task = std::move(tasks.front());
                tasks.pop();
                activeThreads++;
                lock.unlock();

                task();  // 直接执行，避免死锁

                lock.lock();
                activeThreads--;
                if (tasks.empty() && activeThreads == 0) 
                {
                    completionCondition.notify_all();
                }
            } 
            else 
            {
                completionCondition.wait(lock);
            }
        }
    }

    ~ThreadPool() 
    {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers) 
        {
            worker.join();
        }
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queueMutex;
    std::condition_variable condition;
    std::condition_variable completionCondition; // New condition variable for task completion
    bool stop;
    std::atomic<int> activeThreads; // Track number of active threads
};

//定义固定大小的线程池，与CPU核心数相匹配
const int MAX_THREADS = 8; // CPU为8核
ThreadPool* global_thread_pool = nullptr;

//初始化全局线程池
void init_thread_pool() 
{
    if (global_thread_pool == nullptr) 
    {
        global_thread_pool = new ThreadPool(MAX_THREADS);
    }
}

//清理全局线程池
void cleanup_thread_pool() 
{
    if (global_thread_pool != nullptr) 
    {
        delete global_thread_pool;
        global_thread_pool = nullptr;
    }
}

// 内存池管理
class MemoryPool 
{
private:
    std::vector<std::vector<uint64_t>> buffers;
    std::queue<int> available_buffers;
    std::mutex pool_mutex;
    int buffer_size;
    int max_buffers;

public:
    MemoryPool(int buf_size, int max_buf) : buffer_size(buf_size), max_buffers(max_buf) 
    {
        for (int i = 0; i < max_buffers; ++i) 
        {
            buffers.emplace_back(buffer_size, 0);
            available_buffers.push(i);
        }
    }

    std::vector<uint64_t>* get_buffer() 
    {
        std::lock_guard<std::mutex> lock(pool_mutex);
        if (available_buffers.empty()) 
        {
            return nullptr;
        }
        int idx = available_buffers.front();
        available_buffers.pop();
        return &buffers[idx];
    }

    void return_buffer(std::vector<uint64_t>* buffer) 
    {
        std::lock_guard<std::mutex> lock(pool_mutex);
        // 找到buffer对应的索引
        for (int i = 0; i < buffers.size(); ++i) 
        {
            if (&buffers[i] == buffer) 
            {
                available_buffers.push(i);
                break;
            }
        }
    }
};

// 全局内存池
MemoryPool* global_memory_pool = nullptr;

void init_memory_pool(int buffer_size, int max_buffers) 
{
    if (global_memory_pool == nullptr) 
    {
        global_memory_pool = new MemoryPool(buffer_size, max_buffers);
    }
}

void cleanup_memory_pool() 
{
    if (global_memory_pool != nullptr) 
    {
        delete global_memory_pool;
        global_memory_pool = nullptr;
    }
}

struct ModMulTaskArgs 
{
    const std::vector<uint64_t> *a;
    const std::vector<uint64_t> *b;
    std::vector<uint64_t> *c;
    int start;
    int end;
    barrett *m;
};

void *ModMul_worker(void *args) 
{
    ModMulTaskArgs *data = (ModMulTaskArgs*)args;
    const auto &a = *(data->a);
    const auto &b = *(data->b);
    auto &c = *(data->c);
    barrett *m = data->m;

    for (int i = data->start; i < data->end; ++i) 
    {
        c[i] = m->mul(a[i], b[i]);
    }
    return nullptr;
}

struct NTTTaskArgs 
{
    std::vector<uint64_t>* a;
    int n;
    uint64_t p;
    int len;
    uint64_t wn;
    int start_block, end_block;
    barrett* b_ctx;
};

void ntt_worker(NTTTaskArgs* args) 
{
    auto& a = *args->a;
    barrett& b_ctx = *args->b_ctx;
    uint64_t p = args->p;
    int len = args->len;
    uint64_t wn = args->wn;
    int start_block = args->start_block, end_block = args->end_block;

    // 预计算w的幂次
    std::vector<uint64_t> w_table(len / 2);
    uint64_t w_curr = 1;
    for (int j = 0; j < len / 2; ++j) 
    {
        w_table[j] = w_curr;
        w_curr = b_ctx.mul(w_curr, wn);
    }

    for (int b = start_block; b < end_block; ++b)
    {
        int i = b * len;
        uint64_t* A = &a[i];
        for (int jj = 0; jj < len / 2; ++jj) 
        {
            uint64_t w = w_table[jj];
            uint64_t u = A[jj];
            uint64_t v = b_ctx.mul(w, A[jj + len / 2]);
            A[jj] = b_ctx.fast_add(u, v);
            A[jj + len / 2] = b_ctx.fast_sub(u, v);
        }
    }
}

void bit_reverse(std::vector<uint64_t> &a,int n)
{
    for(int i = 1, j = 0; i < n;++i) 
    {
        int bit = n >> 1;
        for(; j >= bit;bit>>=1) 
        {
            j -= bit;
        }
        j += bit;
        if(i<j)
        {
            std::swap(a[i], a[j]);
        }
    }
}

void NTT_iterative(std::vector<uint64_t> &a, int n, uint64_t p, int inv, barrett &b_ctx) 
{
    int g = 3; // 原根
    bit_reverse(a, n);

    for(int len = 2; len <= n; len<<=1) 
    {
        uint64_t wn = power(g, (p - 1) / len, p);
        if(inv == -1) 
        {
            wn = power(wn, p - 2, p);
        }

        int total_blocks = n / len;
        int num_threads = std::min(MAX_THREADS, total_blocks);
        int chunk = (total_blocks + num_threads - 1) / num_threads;

        std::vector<NTTTaskArgs> args(num_threads);
        
        for (int t = 0; t < num_threads; ++t) 
        {
            int start = t * chunk;
            int end = std::min(start + chunk, total_blocks);
            args[t] = NTTTaskArgs{&a, n, p, len, wn, start, end, &b_ctx};
            
            global_thread_pool->enqueue([t, &args]() 
            {
                ntt_worker(&args[t]);
            });
        }
        
        global_thread_pool->waitForAll();
    }

    if(inv == -1) 
    {
        uint64_t inv_n = power(n, p - 2, p);
        
        int num_threads = std::min(MAX_THREADS, n);
        int chunk = (n + num_threads - 1) / num_threads;
        
        for (int t = 0; t < num_threads; ++t) 
        {
            int start = t * chunk;
            int end = std::min(start + chunk, n);
            
            global_thread_pool->enqueue([&a, &b_ctx, inv_n, start, end]() 
            {
                for (int i = start; i < end; ++i) 
                {
                    a[i] = b_ctx.mul(a[i], inv_n);
                }
            });
        }
        
        global_thread_pool->waitForAll();
    }
}

// MPI相关结构体和函数
struct MPICRTData 
{
    std::vector<uint64_t> mods;
    std::vector<std::vector<uint64_t>> a_mods;
    std::vector<std::vector<uint64_t>> b_mods;
    std::vector<std::vector<uint64_t>> res_mods;
    std::vector<barrett> barrett_instances;
    int len;
    int start_mod_idx;
    int end_mod_idx;
};

// 每个进程处理分配给它的模数
void process_mods_in_process(MPICRTData& data) 
{
    for (int k = data.start_mod_idx; k < data.end_mod_idx; ++k) 
    {
        // 执行NTT
        NTT_iterative(data.a_mods[k], data.len, data.mods[k], 1, data.barrett_instances[k]);
        NTT_iterative(data.b_mods[k], data.len, data.mods[k], 1, data.barrett_instances[k]);
        
        // 点乘
        int chunk_size = (data.len + MAX_THREADS - 1) / MAX_THREADS;
        std::vector<ModMulTaskArgs> mul_args(MAX_THREADS);
        
        for (int t = 0; t < MAX_THREADS; ++t) 
        {
            int start = t * chunk_size;
            int end = std::min(start + chunk_size, data.len);
            mul_args[t] = ModMulTaskArgs{&data.a_mods[k], &data.b_mods[k], &data.res_mods[k], start, end, &data.barrett_instances[k]};
            
            global_thread_pool->enqueue([t, &mul_args]() 
            {
                ModMul_worker(&mul_args[t]);
            });
        }
        
        global_thread_pool->waitForAll();
        
        // 逆NTT
        NTT_iterative(data.res_mods[k], data.len, data.mods[k], -1, data.barrett_instances[k]);
    }
}

struct CRTPrecomputed 
{
    std::vector<__uint128_t> Mi;
    std::vector<uint64_t> inv;
    __uint128_t M;
};

CRTPrecomputed crt_precompute(const std::vector<uint64_t>& mods) 
{
    CRTPrecomputed pre;
    pre.M = 1;
    for (uint64_t m : mods) 
    {
        pre.M *= m;
    }
    int k = mods.size();
    pre.Mi.resize(k);
    pre.inv.resize(k);
    for (int i = 0; i < k; ++i) 
    {
        pre.Mi[i] = pre.M / mods[i];
        pre.inv[i] = power(pre.Mi[i] % mods[i], mods[i] - 2, mods[i]);
    }
    return pre;
}

// CRT合并函数
void CRT_combine_parallel(std::vector<uint64_t>& result, const std::vector<std::vector<uint64_t>>& res_mods,const std::vector<uint64_t>& mods,const CRTPrecomputed& pre,uint64_t target_mod,int len) 
{
    int thread_count = MAX_THREADS;
    int block_size = (len + thread_count - 1) / thread_count;
    
    for (int t = 0; t < thread_count; ++t) 
    {
        int l = t * block_size;
        int r = std::min((t + 1) * block_size, len);
        
        global_thread_pool->enqueue([&result, &res_mods, &mods, &pre, target_mod, l, r]() 
        {
            int mod_count = mods.size();
            for (int i = l; i < r; ++i) 
            {
                __uint128_t res = 0;
                for (int k = 0; k < mod_count; ++k) 
                {
                    uint64_t t = res_mods[k][i] * pre.inv[k] % mods[k];
                    __uint128_t term = t * pre.Mi[k] % pre.M;
                    res = (res + term) % pre.M;
                }
                result[i] = (uint64_t)(res % target_mod);
            }
        });
    }
    
    global_thread_pool->waitForAll();
}

uint64_t a[300000],b[300000],ab[300000];
int main(int argc, char *argv[]) 
{
    // 初始化MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    init_thread_pool();
    
    // 定义模数
    std::vector<uint64_t> mods = {1004535809, 1224736769, 469762049, 998244353};
    int total_mods = mods.size();
    
    // 分配模数给各个进程
    int mods_per_process = (total_mods + size - 1) / size;
    int start_mod_idx = rank * mods_per_process;
    int end_mod_idx = std::min(start_mod_idx + mods_per_process, total_mods);
    
    // 只有rank 0进程负责输出结果
    if (rank == 0) 
    {
        std::cout << "使用 " << size << " 个MPI进程，每个进程使用 " << MAX_THREADS << " 个线程" << std::endl;
    }
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 1337006139375617
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 4;
    
    for (int i = test_begin; i <= test_end; ++i) 
    {
        uint64_t n_, p_;

        if (rank == 0) 
        {
            fRead(a, b, &n_, &p_, i);
        }

        MPI_Bcast(&n_, sizeof(uint64_t), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&p_, sizeof(uint64_t), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(a, n_ * sizeof(uint64_t), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(b, n_ * sizeof(uint64_t), MPI_BYTE, 0, MPI_COMM_WORLD);

        memset(ab, 0, sizeof(ab));
        barrett b_ctx(p_);

        int len = 1;
        while (len < 2 * n_) len <<= 1;

        std::vector<uint64_t> a_1(len, 0);
        std::vector<uint64_t> b_1(len, 0);
        std::vector<uint64_t> c(len, 0);
        for (int j = 0; j < n_; ++j) 
        {
            a_1[j] = a[j];
            b_1[j] = b[j];
        }

        double total_time = 0.0;
        for (int repeat = 0; repeat < 50; ++repeat) 
        {
            auto Start = std::chrono::high_resolution_clock::now();

            if (false) 
            {
                if (rank == 0) 
                {
                    std::vector<uint64_t> a_temp = a_1, b_temp = b_1;
                    std::fill(c.begin(), c.end(), 0);

                    NTT_iterative(a_temp, len, p_, 1, b_ctx);
                    NTT_iterative(b_temp, len, p_, 1, b_ctx);

                    int chunk_size = (len + MAX_THREADS - 1) / MAX_THREADS;
                    std::vector<ModMulTaskArgs> mul_args(MAX_THREADS);
                    
                    for (int t = 0; t < MAX_THREADS; ++t) 
                    {
                        int start = t * chunk_size;
                        int end = std::min(start + chunk_size, len);
                        mul_args[t] = ModMulTaskArgs{&a_temp, &b_temp, &c, start, end, &b_ctx};

                        global_thread_pool->enqueue([t, &mul_args]() {
                            ModMul_worker(&mul_args[t]);
                        });
                    }

                    global_thread_pool->waitForAll();
                    NTT_iterative(c, len, p_, -1, b_ctx);
                }
            } 
            else 
            {
                CRTPrecomputed pre = crt_precompute(mods);

                MPICRTData mpi_data;
                mpi_data.mods = mods;
                mpi_data.len = len;
                mpi_data.start_mod_idx = start_mod_idx;
                mpi_data.end_mod_idx = end_mod_idx;

                for (int k = 0; k < total_mods; ++k) 
                {
                    mpi_data.barrett_instances.emplace_back(mods[k]);
                }

                mpi_data.a_mods.resize(total_mods);
                mpi_data.b_mods.resize(total_mods);
                mpi_data.res_mods.resize(total_mods);
                
                for (int k = 0; k < total_mods; ++k) 
                {
                    mpi_data.a_mods[k].assign(len, 0);
                    mpi_data.b_mods[k].assign(len, 0);
                    mpi_data.res_mods[k].assign(len, 0);
                    for (int j = 0; j < n_; ++j) 
                    {
                        mpi_data.a_mods[k][j] = a_1[j] % mods[k];
                        mpi_data.b_mods[k][j] = b_1[j] % mods[k];
                    }
                }

                process_mods_in_process(mpi_data);

                if (rank == 0) 
                {
                    for (int src = 1; src < size; ++src) 
                    {
                        int src_start = src * mods_per_process;
                        int src_end = std::min(src_start + mods_per_process, total_mods);

                        for (int k = src_start; k < src_end; ++k) 
                        {
                            MPI_Recv(mpi_data.res_mods[k].data(), len * sizeof(uint64_t), MPI_BYTE, src, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
                    }
                } 
                else 
                {
                    for (int k = start_mod_idx; k < end_mod_idx; ++k) 
                    {
                        MPI_Send(mpi_data.res_mods[k].data(), len * sizeof(uint64_t), MPI_BYTE, 0, k, MPI_COMM_WORLD);
                    }
                }

                if (rank == 0) 
                {
                    CRT_combine_parallel(c, mpi_data.res_mods, mods, pre, p_, len);
                }
            }

            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::ratio<1,1000>> elapsed = End - Start;
            total_time += elapsed.count();
        }

        if (rank == 0) 
        {
            for (int j = 0; j < 2 * n_ - 1; ++j) 
            {
                ab[j] = c[j];
            }
            fCheck(ab, n_, i);
            std::cout << "average latency for n = " << n_ << " p = " << p_ << " : " << (total_time / 50.0) << " (us)" << std::endl;
            fWrite(ab, n_, i);
        }
    }
    
    cleanup_thread_pool();
    MPI_Finalize();
    return 0;
}