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

class montgomery
{
public:
    montgomery(uint64_t R_, uint64_t N_): N(N_), R(R_) 
    {
        logR = static_cast<int>(std::log2(R));
        uint64_t N_inv = modinv(N, R);
        N_inv_neg = R - N_inv;
        R2 = (__uint128_t(R) * R) % N;
        vec_modulus = vdupq_n_u64(N);
        vec_inv = vdupq_n_u64(N_inv);
        vec_R2 = vdupq_n_u64(R2);
    }

    static uint64_t gcd_extended(uint64_t a, uint64_t b, int64_t &x, int64_t &y) 
    {
        if(a==0) 
        {
            x = 0;
            y = 1;
            return b;
        }
        int64_t x1, y1;
        uint64_t gcd = gcd_extended(b % a, a, x1, y1);
        x = y1 - (b / a) * x1;
        y = x1;
        return gcd;
    }

    static uint64_t modinv(uint64_t a, uint64_t m) 
    {
        int64_t x, y;
        uint64_t g = gcd_extended(a, m, x, y);
        if(g!=1) 
        {
            throw std::invalid_argument("modular inverse does not exist");
        }
        return (x % int64_t(m) + m) % m;
    }

    inline uint64_t REDC(__uint128_t T) 
    {
        uint64_t mask = (1ULL << logR) - 1;
        uint64_t m = (uint64_t(T) & mask) * N_inv_neg & mask;
        __uint128_t t = T + __uint128_t(m) * N;
        t >>= logR;
        if(t>=N) 
        {
            return uint64_t(t - N);
        }
        return uint64_t(t);
    }

    //将普通数转换到Montgomery中
    inline uint64_t toMont(uint64_t x)
    {
        return REDC((__uint128_t)x * R2);
    }

    //将Montgomery数转换到普通数
    inline uint64_t fromMont(uint64_t x) 
    {
        return REDC((__uint128_t)x);
    }

    inline uint64_t ModMul(uint64_t a, uint64_t b) 
    {
        return REDC((__uint128_t)a * b);
    }

    inline uint64_t add(uint64_t a, uint64_t b) const 
    {
        uint64_t res = a + b;
        if (res >= N || res < a) 
        {
            res -= N;
        }
        return res;
    }

    inline uint64_t sub(uint64_t a, uint64_t b) const 
    {
        return (a >= b) ? (a - b) : (N - (b - a));
    }

    inline uint64_t montgomery_add(uint64_t a, uint64_t b, uint64_t p) 
    {
        uint64_t sum = a + b;
        uint64_t mask = -(sum >= p);
        return sum - (mask & p);
    }

    inline uint64_t montgomery_sub(uint64_t a, uint64_t b, uint64_t p) 
    {
        uint64_t diff = a - b;
        uint64_t mask = -(a < b);
        return diff + (mask & p);
    }

    void toMontgomery(std::vector<uint64_t> &vec)
    {
        for (auto &x : vec)
        {
            x = toMont(x);
        }
    }

    void fromMontgomery(std::vector<uint64_t> &vec)
    {
        for (auto &x : vec)
        {
            x = fromMont(x);
        }
    }

    //完整版乘法
    void ModMulSIMD(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b, std::vector<uint64_t> &res)
    {
        size_t n = a.size();
        res.resize(n);
        for (size_t i = 0; i < n; i += 2) 
        {
            //拆分为 32-bit 高低位
            uint32x2_t a_lo = {uint32_t(a[i]), uint32_t(a[i + 1])};
            uint32x2_t a_hi = {uint32_t(a[i] >> 32), uint32_t(a[i + 1] >> 32)};
            uint32x2_t b_lo = {uint32_t(b[i]), uint32_t(b[i + 1])};
            uint32x2_t b_hi = {uint32_t(b[i] >> 32), uint32_t(b[i + 1] >> 32)};

            //计算各部分乘积
            uint64x2_t lo = vmull_u32(a_lo, b_lo);
            uint64x2_t hi = vmull_u32(a_hi, b_hi); 
            uint64x2_t mid1 = vmull_u32(a_hi, b_lo); 
            uint64x2_t mid2 = vmull_u32(a_lo, b_hi);
            uint64x2_t mid = vaddq_u64(mid1, mid2);

            __uint128_t prod0 = (__uint128_t)vgetq_lane_u64(lo, 0);
            prod0 += (__uint128_t)vgetq_lane_u64(mid, 0) << 32;
            prod0 += (__uint128_t)vgetq_lane_u64(hi, 0) << 64;
            res[i] = REDC(prod0);

            __uint128_t prod1 = (__uint128_t)vgetq_lane_u64(lo, 1);
            prod1 += (__uint128_t)vgetq_lane_u64(mid, 1) << 32;
            prod1 += (__uint128_t)vgetq_lane_u64(hi, 1) << 64;
            res[i + 1] = REDC(prod1);
        }
    }

    inline uint64x2_t ModAddSIMD(uint64x2_t a, uint64x2_t b)
    {
        uint64x2_t sum = vaddq_u64(a, b);
        uint64x2_t cmp = vcgeq_u64(sum, vec_modulus);  // sum >= modulus
        return vbslq_u64(cmp, vsubq_u64(sum, vec_modulus), sum);
    }

    inline uint64x2_t ModSubSIMD(uint64x2_t a, uint64x2_t b)
    {
        uint64x2_t diff = vsubq_u64(a, b);
        uint64x2_t cmp = vcgtq_u64(b, a);  // b > a
        return vbslq_u64(cmp, vaddq_u64(diff, vec_modulus), diff);
    }

    inline uint64x2_t ModMulSIMD(uint64x2_t a, uint64x2_t b)
    {
        alignas(16) uint64_t a_raw[2], b_raw[2];
        vst1q_u64(a_raw, a);
        vst1q_u64(b_raw, b);

        uint64_t r0 = ModMul(a_raw[0], b_raw[0]);
        uint64_t r1 = ModMul(a_raw[1], b_raw[1]);

        return (uint64x2_t){r0, r1};
    }
public:
    uint64_t N, R, R2, N_inv_neg;
    int logR;
    uint64x2_t vec_modulus;
    uint64x2_t vec_inv;
    uint64x2_t vec_R2;
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

struct ModMulTaskArgs 
{
    const std::vector<uint64_t> *a;
    const std::vector<uint64_t> *b;
    std::vector<uint64_t> *c;
    int start;
    int end;
    montgomery *m;
};

void *ModMul_worker(void *args) 
{
    ModMulTaskArgs *data = (ModMulTaskArgs*)args;
    const auto &a = *(data->a);
    const auto &b = *(data->b);
    auto &c = *(data->c);
    montgomery *m = data->m;

    for (int i = data->start; i < data->end; ++i) 
    {
        c[i] = m->ModMul(a[i], b[i]);
    }
    return nullptr;
}

// 优化的ModMul_worker函数，使用SIMD
void ModMul_worker_simd(ModMulTaskArgs* args) 
{
    const auto& a = *(args->a);
    const auto& b = *(args->b);
    auto& c = *(args->c);
    montgomery* m = args->m;
    int start = args->start;
    int end = args->end;
    
    // SIMD处理，每次处理2个元素
    int i = start;
    for (; i + 1 < end; i += 2) 
    {
        uint64x2_t va = vld1q_u64(&a[i]);
        uint64x2_t vb = vld1q_u64(&b[i]);
        uint64x2_t result = m->ModMulSIMD(va, vb);
        vst1q_u64(&c[i], result);
    }
    
    // 处理剩余元素
    for (; i < end; ++i) 
    {
        c[i] = m->ModMul(a[i], b[i]);
    }
}

struct NTTTaskArgs 
{
    std::vector<uint64_t>* a;
    int n;
    uint64_t p;
    int len;
    uint64_t wnR;
    int start_block, end_block;
    montgomery* m;
};

void ntt_worker(NTTTaskArgs* args) 
{
    auto& a = *args->a;
    montgomery& m = *args->m;
    uint64_t p = args->p;
    int len = args->len;
    uint64_t wnR = args->wnR;
    int start_block = args->start_block, end_block = args->end_block;

    std::vector<uint64_t> w_table(len / 2);
    uint64_t w_mont = m.toMont(1);
    for (int j = 0; j < len / 2; ++j) 
    {
        w_table[j] = w_mont;
        w_mont = m.ModMul(w_mont, wnR);
    }

    for (int b = start_block; b < end_block; ++b)
    {
        int i = b * len;
        uint64_t* A = &a[i];
        for (int jj = 0; jj < len / 2; ++jj) 
        {
            uint64_t w = w_table[jj];
            uint64_t u = A[jj];
            uint64_t v = m.ModMul(w, A[jj + len / 2]);
            A[jj] = m.montgomery_add(u, v, p);
            A[jj + len / 2] = m.montgomery_sub(u, v, p);
        }
    }
}

// 添加SIMD优化的NTT工作函数
void ntt_worker_simd(NTTTaskArgs* args) 
{
    auto& a = *args->a;
    montgomery& m = *args->m;
    uint64_t p = args->p;
    int len = args->len;
    uint64_t wnR = args->wnR;
    int start_block = args->start_block, end_block = args->end_block;

    // 预计算w表，并转换为SIMD友好的格式
    std::vector<uint64_t> w_table(len / 2);
    uint64_t w_mont = m.toMont(1);
    for (int j = 0; j < len / 2; ++j) 
    {
        w_table[j] = w_mont;
        w_mont = m.ModMul(w_mont, wnR);
    }

    uint64x2_t vec_p = vdupq_n_u64(p);
    
    for (int b = start_block; b < end_block; ++b)
    {
        int i = b * len;
        uint64_t* A = &a[i];
        
        // SIMD优化的蝶形运算，每次处理2对数据
        for (int jj = 0; jj < len / 2; jj += 2) 
        {
            if (jj + 1 < len / 2) {
                // 加载权重
                uint64x2_t w = vld1q_u64(&w_table[jj]);
                
                // 加载数据
                uint64x2_t u = vld1q_u64(&A[jj]);
                uint64x2_t v_data = vld1q_u64(&A[jj + len / 2]);
                
                // v = w * A[jj + len/2] (Montgomery乘法)
                uint64x2_t v = m.ModMulSIMD(w, v_data);
                
                // 蝶形运算
                uint64x2_t add_result = m.ModAddSIMD(u, v);
                uint64x2_t sub_result = m.ModSubSIMD(u, v);
                
                // 存储结果
                vst1q_u64(&A[jj], add_result);
                vst1q_u64(&A[jj + len / 2], sub_result);
            } else {
                // 处理剩余的单个元素
                uint64_t w = w_table[jj];
                uint64_t u = A[jj];
                uint64_t v = m.ModMul(w, A[jj + len / 2]);
                A[jj] = m.montgomery_add(u, v, p);
                A[jj + len / 2] = m.montgomery_sub(u, v, p);
            }
        }
    }
}

//下面是迭代的写法，据说更快（想必更快）
//位转换置换，原来位转换置换有问题，索引交换错误
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

void NTT_iterative(std::vector<uint64_t> &a, int n, uint64_t p, int inv, montgomery &m) 
{
    int g = 3; // 原根
    bit_reverse(a, n); // 位反转置换

    // 确保线程池已初始化
    // init_thread_pool();

    for(int len = 2; len <= n; len<<=1) 
    {
        uint64_t wn = power(g, (p - 1) / len, p);
        if(inv == -1) 
        {
            wn = power(wn, p - 2, p);
        }
        uint64_t wnR = m.toMont(wn);

        int total_blocks = n / len;
        int num_threads = std::min(MAX_THREADS, total_blocks);
        int chunk = (total_blocks + num_threads - 1) / num_threads;

        std::vector<NTTTaskArgs> args(num_threads);
        
        // 使用线程池分发任务
        for (int t = 0; t < num_threads; ++t) 
        {
            int start = t * chunk;
            int end = std::min(start + chunk, total_blocks);
            args[t] = NTTTaskArgs{&a, n, p, len, wnR, start, end, &m};
            
            global_thread_pool->enqueue([t, &args]() 
            {
                ntt_worker(&args[t]);
            });
        }
        
        // 等待所有NTT任务完成
        global_thread_pool->waitForAll();
    }

    if(inv == -1) 
    {
        uint64_t inv_n = power(n, p - 2, p);
        uint64_t invR = m.toMont(inv_n);
        
        // 并行处理除法
        int num_threads = std::min(MAX_THREADS, n);
        int chunk = (n + num_threads - 1) / num_threads;
        
        std::vector<std::pair<int, int>> ranges(num_threads);
        for (int t = 0; t < num_threads; ++t) 
        {
            int start = t * chunk;
            int end = std::min(start + chunk, n);
            ranges[t] = {start, end};
            
            global_thread_pool->enqueue([&a, &m, invR, start, end]() 
            {
                for (int i = start; i < end; ++i) 
                {
                    a[i] = m.ModMul(a[i], invR);
                }
            });
        }
        
        // 等待所有除法任务完成
        global_thread_pool->waitForAll();
    }
}

// 优化的NTT函数，使用SIMD版本的工作函数
void NTT_iterative_simd(std::vector<uint64_t> &a, int n, uint64_t p, int inv, montgomery &m) 
{
    int g = 3; // 原根
    bit_reverse(a, n); // 位反转置换

    for(int len = 2; len <= n; len<<=1) 
    {
        uint64_t wn = power(g, (p - 1) / len, p);
        if(inv == -1) 
        {
            wn = power(wn, p - 2, p);
        }
        uint64_t wnR = m.toMont(wn);

        int total_blocks = n / len;
        int num_threads = std::min(MAX_THREADS, total_blocks);
        int chunk = (total_blocks + num_threads - 1) / num_threads;

        std::vector<NTTTaskArgs> args(num_threads);

        for (int t = 0; t < num_threads; ++t) 
        {
            int start = t * chunk;
            int end = std::min(start + chunk, total_blocks);
            args[t] = NTTTaskArgs{&a, n, p, len, wnR, start, end, &m};
            
            global_thread_pool->enqueue([t, &args]() 
            {
                ntt_worker_simd(&args[t]); // 使用SIMD版本
            });
        }

        global_thread_pool->waitForAll();
    }

    if(inv == -1) 
    {
        uint64_t inv_n = power(n, p - 2, p);
        uint64_t invR = m.toMont(inv_n);

        int num_threads = std::min(MAX_THREADS, n);
        int chunk = (n + num_threads - 1) / num_threads;
        
        std::vector<std::pair<int, int>> ranges(num_threads);
        for (int t = 0; t < num_threads; ++t) 
        {
            int start = t * chunk;
            int end = std::min(start + chunk, n);
            ranges[t] = {start, end};
            
            global_thread_pool->enqueue([&a, &m, invR, start, end]() 
            {
                // SIMD优化的逆变换后处理
                int i = start;
                uint64x2_t vinv = vdupq_n_u64(invR);
                for (; i + 1 < end; i += 2) 
                {
                    uint64x2_t va = vld1q_u64(&a[i]);
                    uint64x2_t result = m.ModMulSIMD(va, vinv);
                    vst1q_u64(&a[i], result);
                }
                // 处理剩余元素
                for (; i < end; ++i) 
                {
                    a[i] = m.ModMul(a[i], invR);
                }
            });
        }

        global_thread_pool->waitForAll();
    }
}

struct CRTTaskArgs
{
    std::vector<uint64_t> *a;
    std::vector<uint64_t> *b;
    std::vector<uint64_t> *result;
    uint64_t p;
    int len;
    montgomery *m;
};

void CRT_worker(CRTTaskArgs* args) 
{
    auto& a = *args->a;
    auto& b = *args->b;
    auto& result = *args->result;
    montgomery& m = *args->m;
    int len = args->len;
    uint64_t p = args->p;

    NTT_iterative(a, len, p, 1, m);
    NTT_iterative(b, len, p, 1, m);

    // 用线程池进行点乘优化
    int num_chunks = MAX_THREADS;
    int chunk_size = (len + num_chunks - 1) / num_chunks;
    
    std::vector<ModMulTaskArgs> mul_args(num_chunks);
    
    for (int t = 0; t < num_chunks; ++t) 
    {
        int start = t * chunk_size;
        int end = std::min(start + chunk_size, len);
        mul_args[t] = ModMulTaskArgs{&a, &b, &result, start, end, &m};
        
        global_thread_pool->enqueue([t, &mul_args]() 
        {
            ModMul_worker(&mul_args[t]);
        });
    }
    
    // 等待所有点乘任务完成
    global_thread_pool->waitForAll();
    
    NTT_iterative(result, len, p, -1, m);
    m.fromMontgomery(result);
}

void CRT_worker_simd(CRTTaskArgs* args) 
{
    auto& a = *args->a;
    auto& b = *args->b;
    auto& result = *args->result;
    montgomery& m = *args->m;
    int len = args->len;
    uint64_t p = args->p;

    NTT_iterative_simd(a, len, p, 1, m);  // 使用SIMD版本
    NTT_iterative_simd(b, len, p, 1, m);  // 使用SIMD版本

    int num_chunks = MAX_THREADS;
    int chunk_size = (len + num_chunks - 1) / num_chunks;
    
    std::vector<ModMulTaskArgs> mul_args(num_chunks);
    
    for (int t = 0; t < num_chunks; ++t) 
    {
        int start = t * chunk_size;
        int end = std::min(start + chunk_size, len);
        mul_args[t] = ModMulTaskArgs{&a, &b, &result, start, end, &m};
        
        global_thread_pool->enqueue([t, &mul_args]() 
        {
            ModMul_worker_simd(&mul_args[t]); // 使用SIMD版本
        });
    }

    global_thread_pool->waitForAll();
    
    NTT_iterative_simd(result, len, p, -1, m); // 使用SIMD版本
    m.fromMontgomery(result);
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

struct CRTCombineArgs 
{
    std::vector<uint64_t>* result;
    const std::vector<std::vector<uint64_t>>* res_mods;
    const std::vector<uint64_t>* mods;
    const CRTPrecomputed* pre;
    uint64_t target_mod;
    int l, r;
};

void CRT_combine_worker(CRTCombineArgs* args) 
{
    int l = args->l, r = args->r;
    const auto& res_mods = *(args->res_mods);
    const auto& mods = *(args->mods);
    auto& result = *(args->result);
    const auto& pre = *(args->pre);
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
        result[i] = (uint64_t)(res % args->target_mod);
    }
}

//基4的位逆序置换函数
void bit_reverse_radix4(std::vector<uint64_t> &a, int n) 
{
    int log4n = 0;
    int temp = n;
    while (temp > 1) 
    {
        temp >>= 2;
        log4n++;
    }
    for (int i = 0; i < n; ++i) 
    {
        int reversed = 0;
        int num = i;
        for (int j = 0; j < log4n; ++j) 
        {
            reversed = (reversed << 2) | (num & 3);
            num >>= 2;
        }
        if (i < reversed) 
        {
            std::swap(a[i], a[reversed]);
        }
    }
}

void NTT_radix4(std::vector<uint64_t> &a, int n, int p, int inv,montgomery &m) 
{
    int g = 3;
    bit_reverse_radix4(a, n);

    for (int len = 4; len <= n; len <<= 2) 
    {
        int step = len >> 2;
        uint64_t w = power(g, (p - 1) / len, p);
        if (inv == -1) 
        {
            w = power(w, p - 2, p);
        }
        uint64_t imag = power(w, step, p);

        uint64_t wR = m.toMont(w);
        uint64_t w2R = m.ModMul(wR, wR);
        uint64_t w3R = m.ModMul(w2R, wR);
        uint64_t imagR = m.toMont(imag);

        for (int i = 0; i < n; i += len) 
        {
            uint64_t w1R_current = m.toMont(1);
            uint64_t w2R_lane_current = m.toMont(1);
            uint64_t w3R_lane_current = m.toMont(1);

            for (int j = 0; j < step; j += 2) 
            {
                if (j + 1 >= step) break; //处理奇数step时的边界

                //生成旋转因子向量
                uint64_t w1R_j1 = m.ModMul(w1R_current, wR);
                uint64x2_t w1R_vec = vcombine_u64(vcreate_u64(w1R_current), vcreate_u64(w1R_j1));

                uint64_t w2R_j1 = m.ModMul(w2R_lane_current, w2R);
                uint64x2_t w2R_lane_vec = vcombine_u64(vcreate_u64(w2R_lane_current), vcreate_u64(w2R_j1));

                uint64_t w3R_j1 = m.ModMul(w3R_lane_current, w3R);
                uint64x2_t w3R_lane_vec = vcombine_u64(vcreate_u64(w3R_lane_current), vcreate_u64(w3R_j1));

                uint64x2_t imagR_vec = vdupq_n_u64(imagR);

                //加载数据
                uint64x2_t a0 = vld1q_u64(&a[i + j]);
                uint64x2_t a1 = vld1q_u64(&a[i + j + step]);
                uint64x2_t a2 = vld1q_u64(&a[i + j + 2 * step]);
                uint64x2_t a3 = vld1q_u64(&a[i + j + 3 * step]);

                //计算t1, t2, t3
                uint64x2_t t1 = m.ModMulSIMD(a1, w1R_vec);
                uint64x2_t t2 = m.ModMulSIMD(a2, w2R_lane_vec);
                uint64x2_t t3 = m.ModMulSIMD(a3, w3R_lane_vec);

                //计算中间项
                uint64x2_t t1i = m.ModMulSIMD(t1, imagR_vec);
                uint64x2_t t3i = m.ModMulSIMD(t3, imagR_vec);

                //计算y0-y3
                uint64x2_t y0 = m.ModAddSIMD(m.ModAddSIMD(m.ModAddSIMD(a0, t1), t2), t3);
                uint64x2_t y1 = m.ModSubSIMD(m.ModSubSIMD(m.ModAddSIMD(a0, t1i), t2), t3i);
                uint64x2_t y2 = m.ModSubSIMD(m.ModAddSIMD(m.ModSubSIMD(a0, t1), t2), t3);
                uint64x2_t y3 = m.ModAddSIMD(m.ModSubSIMD(m.ModSubSIMD(a0, t1i), t2), t3i);

                //存储结果
                vst1q_u64(&a[i + j], y0);
                vst1q_u64(&a[i + j + step], y1);
                vst1q_u64(&a[i + j + 2 * step], y2);
                vst1q_u64(&a[i + j + 3 * step], y3);

                //更新旋转因子
                w1R_current = m.ModMul(w1R_current, m.ModMul(wR, wR));
                w2R_lane_current = m.ModMul(w2R_lane_current, m.ModMul(w2R, w2R));
                w3R_lane_current = m.ModMul(w3R_lane_current, m.ModMul(w3R, w3R));
            }
            //处理剩余的j（step为奇数）
            if (step % 2 != 0) 
            {
                int j = step - 1;
                uint64_t a0R = a[i + j];
                uint64_t a1R = a[i + j + step];
                uint64_t a2R = a[i + j + 2 * step];
                uint64_t a3R = a[i + j + 3 * step];

                uint64_t t1R = m.ModMul(a1R, w1R_current);
                uint64_t t2R = m.ModMul(a2R, w2R_lane_current);
                uint64_t t3R = m.ModMul(a3R, w3R_lane_current);

                uint64_t t1iR = m.ModMul(t1R, imagR);
                uint64_t t3iR = m.ModMul(t3R, imagR);

                uint64_t y0R = m.add(m.add(m.add(a0R, t1R), t2R), t3R);
                uint64_t y1R = m.sub(m.sub(m.add(a0R, t1iR), t2R), t3iR);
                uint64_t y2R = m.sub(m.add(m.sub(a0R, t1R), t2R), t3R);
                uint64_t y3R = m.add(m.sub(m.sub(a0R, t1iR), t2R), t3iR);

                a[i + j] = y0R;
                a[i + j + step] = y1R;
                a[i + j + 2 * step] = y2R;
                a[i + j + 3 * step] = y3R;
            }
        }
    }
    if (inv == -1) 
    {
        // int inv_n = power(n, p - 2, p);
        uint64_t inv_n = power(n, p - 2, p);
        uint64_t invR = m.toMont(inv_n);
        uint64x2_t invR_vec = vdupq_n_u64(invR);
        for (size_t i = 0; i < a.size(); i += 2) 
        {
            uint64x2_t x = vld1q_u64(&a[i]);
            x = m.ModMulSIMD(x, invR_vec);
            vst1q_u64(&a[i], x);
        }
    }
}

// int a[300000], b[300000];
uint64_t a[300000],b[300000],ab[300000];
int main(int argc, char *argv[])
{
    init_thread_pool();
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 1337006139375617
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 4;
    for(int i = test_begin; i <= test_end; ++i){
        long double ans = 0;
        uint64_t n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        // uint64_t R = 1ULL << 32;
        //换大数有利于减少REDC的次数
        uint64_t R = 1ULL << 63;
        montgomery m(R, p_);
        //多模数分解
        //基本能确定是四个模数的问题了，这四个模数还是太小了
        std::vector<uint64_t> mods = {1004535809, 1224736769, 469762049, 998244353};
        CRTPrecomputed pre = crt_precompute(mods);
        std::vector<montgomery> montgomery_instances = {
            montgomery(R, mods[0]),
            montgomery(R, mods[1]),
            montgomery(R, mods[2]),
            montgomery(R, mods[3])
        };
        int len = 1;
        while(len<2*n_)
        {
            len <<= 1;
        }
        std::vector<uint64_t> a_1(len, 0);
        std::vector<uint64_t> b_1(len, 0);
        for (int i = 0; i < n_; ++i)
        {
            a_1[i] = a[i];
            b_1[i] = b[i];
        }
        std::vector<uint64_t> c(len, 0);
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
        //小模数直接用朴素，等大模数再用CRT多线程
        if(p_<(1ULL<<50))
        {
            m.toMontgomery(a_1);
            m.toMontgomery(b_1);
            // NTT_iterative_simd(a_1, len, p_, 1, m);
            // NTT_iterative_simd(b_1, len, p_, 1, m);
            NTT_iterative(a_1, len, p_, 1, m);
            NTT_iterative(b_1, len, p_, 1, m);
            // 使用线程池进行点乘
            int chunk_size = (len + MAX_THREADS - 1) / MAX_THREADS;
            std::vector<ModMulTaskArgs> mul_args(MAX_THREADS);
            
            for (int t = 0; t < MAX_THREADS; ++t) 
            {
                int start = t * chunk_size;
                int end = std::min(start + chunk_size, len);
                mul_args[t] = ModMulTaskArgs{&a_1, &b_1, &c, start, end, &m};
                
                global_thread_pool->enqueue([t, &mul_args]() 
                {
                    ModMul_worker_simd(&mul_args[t]);
                });
            }
            
            // 等待所有点乘任务完成
            global_thread_pool->waitForAll();
            // NTT_iterative_simd(c, len, p_, -1, m);
            NTT_iterative(c, len, p_, -1, m);
            m.fromMontgomery(c);
        }
        else
        {
            std::vector<std::vector<uint64_t>> a_mods(4, std::vector<uint64_t>(len, 0));
            std::vector<std::vector<uint64_t>> b_mods(4, std::vector<uint64_t>(len, 0));
            std::vector<std::vector<uint64_t>> res_mods(4);
            
            // 准备数据
            for (int k = 0; k < 4; ++k) 
            {
                auto& m_k = montgomery_instances[k];
                for (int j = 0; j < n_; ++j) 
                {
                    a_mods[k][j] = m_k.toMont(a_1[j]);
                    b_mods[k][j] = m_k.toMont(b_1[j]);
                }
                res_mods[k].resize(len, 0);
            }
            
            // 使用线程池并行处理4个模数的NTT
            std::vector<CRTTaskArgs> crt_args(4);
            for (int k = 0; k < 4; ++k) 
            {
                crt_args[k] = CRTTaskArgs{&a_mods[k], &b_mods[k], &res_mods[k], mods[k], len, &montgomery_instances[k]};
                // CRT_worker_simd(&crt_args[k]);
                CRT_worker(&crt_args[k]);
            }
            
            // 等待所有CRT任务完成
            global_thread_pool->waitForAll();
            
            // CRT合并
            int thread_count = MAX_THREADS;
            std::vector<CRTCombineArgs> crt_thread_data(thread_count);
            
            int block_size = (len + thread_count - 1) / thread_count;
            for (int t = 0; t < thread_count; ++t) 
            {
                int l = t * block_size;
                int r = std::min((t + 1) * block_size, len);
                crt_thread_data[t] = CRTCombineArgs{&c, &res_mods, &mods, &pre, p_, l, r};
                
                global_thread_pool->enqueue([t, &crt_thread_data]() 
                {
                    CRT_combine_worker(&crt_thread_data[t]);
                });
            }
            
            // 等待所有CRT合并任务完成
            global_thread_pool->waitForAll();
        }
        auto End = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 2 * n_ - 1; ++i)
        {
            ab[i] = c[i];
        }
        // fWrite(ab, len / 2, "true_result_mods4_", i);
        std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout<<"average latency for n = "<<n_<<" p = "<<p_<<" : "<<ans<<" (us) "<<std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
    }
    cleanup_thread_pool();
    return 0;
}
