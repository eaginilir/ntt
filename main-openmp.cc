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
// #include <arm_neon.h>
#include <pthread.h>
#include <stdint.h>
#include <shared_mutex>
#include <cassert>


// 可以自行添加需要的头文件

void fRead(uint64_t *a, uint64_t *b, uint64_t *n, uint64_t *p, int input_id){
    std::string strin = "D:\\study\\concurrent processing\\NTT\\ntt\\nttdata\\4.in";
    std::ifstream fin(strin);

    if (!fin.is_open()) {
        std::cerr << "❌ 文件打不开，路径是: " << strin << std::endl;
        return;
    }

    fin >> *n >> *p;
    if (fin.fail()) {
        std::cerr << "❌ 无法读取 n 和 p，请检查文件格式。" << std::endl;
        return;
    }

    std::cout << "✅ 成功读取 n = " << *n << ", p = " << *p << std::endl;

    for (int i = 0; i < *n; i++) {
        fin >> a[i];
    }
    for (int i = 0; i < *n; i++) {
        fin >> b[i];
    }
}


void fCheck(uint64_t *ab, int n, int input_id){
    // 判断多项式乘法结果是否正确
    std::string logout = "files/log.out";
    std::string strin = "D:\\study\\concurrent processing\\NTT\\ntt\\nttdata\\4.out";
    std::ifstream fin;
    std::ofstream fout(logout, std::ios::app);
    fin.open(strin, std::ios::in);

    if (!fin.is_open()) {
        std::cerr << "❌ 无法打开输出文件: " << strin << std::endl;
        return;
    }

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
    std::string strout = "D:\\study\\concurrent processing\\NTT\\ntt\\files\\4.out";
    std::ofstream fout;
    fout.open(strout, std::ios::out);
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
        // vec_modulus = vdupq_n_u64(N);
        // vec_inv = vdupq_n_u64(N_inv);
        // vec_R2 = vdupq_n_u64(R2);
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
        #pragma omp parallel for num_threads(8)
        for (auto &x : vec)
        {
            x = toMont(x);
        }
    }

    void fromMontgomery(std::vector<uint64_t> &vec)
    {
        #pragma omp parallel for num_threads(8)
        for (auto &x : vec)
        {
            x = fromMont(x);
        }
    }

    // //完整版乘法
    // void ModMulSIMD(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b, std::vector<uint64_t> &res)
    // {
    //     size_t n = a.size();
    //     res.resize(n);
    //     for (size_t i = 0; i < n; i += 2) 
    //     {
    //         //拆分为 32-bit 高低位
    //         uint32x2_t a_lo = {uint32_t(a[i]), uint32_t(a[i + 1])};
    //         uint32x2_t a_hi = {uint32_t(a[i] >> 32), uint32_t(a[i + 1] >> 32)};
    //         uint32x2_t b_lo = {uint32_t(b[i]), uint32_t(b[i + 1])};
    //         uint32x2_t b_hi = {uint32_t(b[i] >> 32), uint32_t(b[i + 1] >> 32)};

    //         //计算各部分乘积
    //         uint64x2_t lo = vmull_u32(a_lo, b_lo);
    //         uint64x2_t hi = vmull_u32(a_hi, b_hi); 
    //         uint64x2_t mid1 = vmull_u32(a_hi, b_lo); 
    //         uint64x2_t mid2 = vmull_u32(a_lo, b_hi);
    //         uint64x2_t mid = vaddq_u64(mid1, mid2);

    //         __uint128_t prod0 = (__uint128_t)vgetq_lane_u64(lo, 0);
    //         prod0 += (__uint128_t)vgetq_lane_u64(mid, 0) << 32;
    //         prod0 += (__uint128_t)vgetq_lane_u64(hi, 0) << 64;
    //         res[i] = REDC(prod0);

    //         __uint128_t prod1 = (__uint128_t)vgetq_lane_u64(lo, 1);
    //         prod1 += (__uint128_t)vgetq_lane_u64(mid, 1) << 32;
    //         prod1 += (__uint128_t)vgetq_lane_u64(hi, 1) << 64;
    //         res[i + 1] = REDC(prod1);
    //     }
    // }

    // inline uint64x2_t ModAddSIMD(uint64x2_t a, uint64x2_t b)
    // {
    //     uint64x2_t sum = vaddq_u64(a, b);
    //     uint64x2_t cmp = vcgeq_u64(sum, vec_modulus);  // sum >= modulus
    //     return vbslq_u64(cmp, vsubq_u64(sum, vec_modulus), sum);
    // }

    // inline uint64x2_t ModSubSIMD(uint64x2_t a, uint64x2_t b)
    // {
    //     uint64x2_t diff = vsubq_u64(a, b);
    //     uint64x2_t cmp = vcgtq_u64(b, a);  // b > a
    //     return vbslq_u64(cmp, vaddq_u64(diff, vec_modulus), diff);
    // }

    // inline uint64x2_t ModMulSIMD(uint64x2_t a, uint64x2_t b)
    // {
    //     alignas(16) uint64_t a_raw[2], b_raw[2];
    //     vst1q_u64(a_raw, a);
    //     vst1q_u64(b_raw, b);

    //     uint64_t r0 = ModMul(a_raw[0], b_raw[0]);
    //     uint64_t r1 = ModMul(a_raw[1], b_raw[1]);

    //     return (uint64x2_t){r0, r1};
    // }
public:
    uint64_t N, R, R2, N_inv_neg;
    int logR;
    // uint64x2_t vec_modulus;
    // uint64x2_t vec_inv;
    // uint64x2_t vec_R2;
};

const int max_thread = 8;//CPU为8核

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

void NTT_iterative(std::vector<uint64_t> &a, int n, uint64_t p, int inv,montgomery &m)
{
    int g = 3; //原根
    bit_reverse(a, n); //位反转置换

    for(int len = 2; len <= n;len<<=1) 
    {
        uint64_t wn = power(g, (p - 1) / len, p);
        if(inv==-1)
        {
            wn = power(wn, p - 2, p);
        }
        uint64_t wnR= m.toMont(wn);

        int total_blocks = n / len;

        #pragma omp parallel for num_threads(max_thread) schedule(static)
        for (int b = 0; b < total_blocks; ++b)
        {
            int i = b * len;
            uint64_t w_mont = m.toMont(1);
            for (int j = 0; j < len / 2; ++j)
            {
                uint64_t u = a[i + j];
                uint64_t v = m.ModMul(w_mont, a[i + j + len / 2]);
                a[i + j] = (u + v) % p;
                a[i + j + len / 2] = (u + p - v) % p;
                w_mont = m.ModMul(w_mont, wnR);
            }
        }
    }

    if(inv == -1) 
    {
        uint64_t inv_n = power(n, p - 2, p);
        uint64_t invR = m.toMont(inv_n);
        for (uint64_t &x : a)
        {
            x = m.ModMul(x, invR);
        }
    }
}

void CRT_parallel(std::vector<uint64_t>& a, std::vector<uint64_t>& b,std::vector<uint64_t>& result, int len, uint64_t p, montgomery& m)
{
    m.toMontgomery(a);
    m.toMontgomery(b);

    NTT_iterative(a, len, p, 1, m);
    NTT_iterative(b, len, p, 1, m);

    #pragma omp parallel for num_threads(max_thread)
    for (int i = 0; i < len; ++i)
    {
        result[i] = m.ModMul(a[i], b[i]);
    }

    NTT_iterative(result, len, p, -1, m);
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

void CRT_combine_parallel(std::vector<uint64_t>& result,const std::vector<std::vector<uint64_t>>& res_mods,const std::vector<uint64_t>& mods,const CRTPrecomputed& pre,uint64_t target_mod)
{
    int n = result.size();
    int mod_count = mods.size();

    #pragma omp parallel for num_threads(max_thread)
    for (int i = 0; i < n; ++i)
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

// int a[300000], b[300000];
uint64_t a[300000],b[300000],ab[300000];
int main(int argc, char *argv[])
{
    
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 1337006139375617
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 4;
    int test_end = 4;
    for(int i = test_begin; i <= test_end; ++i){
        long double ans = 0;
        uint64_t n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        // uint64_t R = 1ULL << 32;
        //换大数有利于减少REDC的次数
        uint64_t R = 1ULL << 60;
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
            NTT_iterative(a_1, len, p_, 1, m);
            NTT_iterative(b_1, len, p_, 1, m);
            // m.ModMulSIMD(a_1, b_1, c);
            #pragma omp parallel for num_threads(8)
            for (int i = 0; i < len; ++i) 
            {
                c[i] = m.ModMul(a_1[i], b_1[i]);
            }
            NTT_iterative(c, len, p_, -1, m);
            m.fromMontgomery(c);
        }
        else
        {
            std::vector<std::vector<uint64_t>> a_mods(4, std::vector<uint64_t>(len, 0));
            std::vector<std::vector<uint64_t>> b_mods(4, std::vector<uint64_t>(len, 0));
            std::vector<std::vector<uint64_t>> res_mods(4, std::vector<uint64_t>(len, 0));

            #pragma omp parallel for num_threads(4)
            for (int k = 0; k < 4; ++k)
            {
                for (int j = 0; j < n_; ++j) 
                {
                    a_mods[k][j] = a[j];
                    b_mods[k][j] = b[j];
                }
                CRT_parallel(a_mods[k], b_mods[k], res_mods[k], len, mods[k], montgomery_instances[k]);
            }

            CRT_combine_parallel(c, res_mods, mods, pre, p_);
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
    return 0;
}
