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
// 可以自行添加需要的头文件

void fRead(int *a, int *b, int *n, int *p, int input_id){
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

void fCheck(int *ab, int n, int input_id){
    // 判断多项式乘法结果是否正确
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

void poly_multiply(int *a, int *b, int *ab, int n, int p){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            ab[i+j]=(1LL * a[i] * b[j] % p + ab[i+j]) % p;
        }
    }
}

//这里我们分开尝试递归和迭代写法，先从递归开始吧（这个好写）

//快速幂
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

//不用看下面这个递归版本，这个版本还没测试过
std::vector<int> NTT_recursive(std::vector<int> &a,int n,int p,int inv)
{
    if(n==1)//长度为0，直接返回
    {
        return a;
    }
    int g = 3; //原根
    int wn = power(g, (p - 1) / n, p);
    if (inv == -1)
    {
        wn = power(wn, p - 2, p); //模逆元
    }
    int half = n / 2;
    std::vector<int> a_e(half);
    std::vector<int> a_o(half);
    for (int i = 0; i < n; i += 2)
    {
        a_e[i / 2] = a[i];
        a_o[i / 2] = a[i + 1];
    }
    a_e = NTT_recursive(a_e, half, p, inv);
    a_o = NTT_recursive(a_o, half, p, inv);
    std::vector<int> result(n);
    int w = 1;
    for (int i = 0; i < half; ++i)
    {
        int u = a_e[i];
        int v = 1LL * w * a_o[i] % p;
        result[i] = (u + v) % p;
        result[i + half] = (u - v + p) % p;
        w = 1LL * w * wn % p;
    }
    return result;
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

    uint64_t REDC(__uint128_t T) 
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

    uint64_t ModMul(uint64_t a, uint64_t b) 
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

    //粗略乘法也能过四个测试样例，挺神奇的，可能因为数据都是
    // void ModMulSIMD(const std::vector<uint64_t> &a,const std::vector<uint64_t> &b,std::vector<uint64_t> &res)
    // {
    //     size_t n = a.size();
    //     res.resize(n);
    //     for (size_t i = 0; i < n; i += 2) 
    //     {
    //     //拆成32-bit低高位
    //     uint32x2_t a_lo = {uint32_t(a[i]), uint32_t(a[i + 1])};
    //     uint32x2_t a_hi = {uint32_t(a[i] >> 32), uint32_t(a[i + 1] >> 32)};
    //     uint32x2_t b_lo = {uint32_t(b[i]), uint32_t(b[i + 1])};
    //     uint32x2_t b_hi = {uint32_t(b[i] >> 32), uint32_t(b[i + 1] >> 32)};

    //     uint64x2_t res_lo = vmull_u32(a_lo, b_lo); //a_lo * b_lo
    //     uint64x2_t res_hi = vmull_u32(a_hi, b_hi); //a_hi * b_hi
    //     uint64x2_t result = vaddq_u64(res_lo, res_hi); //简化估算，不完整乘法（适用于测试）

    //     res[i] = REDC((__uint128_t)vgetq_lane_u64(result, 0));
    //     res[i + 1] = REDC((__uint128_t)vgetq_lane_u64(result, 1));
    //     }
    // }

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

void NTT_iterative(std::vector<uint64_t> &a, int n, int p, int inv)
{
    uint64_t R = 1ULL << 32;
    montgomery m(R,p);

    int g = 3; //原根
    bit_reverse(a, n); //位反转置换

    for(int len = 2; len <= n;len<<=1) 
    {
        int wn = power(g, (p - 1) / len, p);
        if(inv==-1)
        {
            wn = power(wn, p - 2, p);
        }

        for(int i = 0; i < n;i+=len) 
        {
            int w_mont = 1;
            for (int j = 0; j < len / 2;++j) 
            {
                uint64_t u = a[i + j];
                // int v = 1LL * w * a[i + j + len / 2] % p;
                uint64_t v = m.ModMul(w_mont, a[i + j + len/2]);
                a[i + j] = (u + v) % p;
                a[i + j + len / 2] = (u - v + p) % p;
                w_mont = m.ModMul(w_mont, wn);
            }
        }
    }

    if(inv == -1) 
    {
        int inv_n = power(n, p - 2, p);
        for (uint64_t &x : a)
        {
            // x = 1LL * x * inv_n % p;
            x = m.ModMul(x, inv_n);
        }
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
        int inv_n = power(n, p - 2, p);
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

int a[300000],b[300000], ab[300000];
int main(int argc, char *argv[])
{
    
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 3;
    for(int i = test_begin; i <= test_end; ++i){
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        uint64_t R = 1ULL << 32;
        montgomery m(R,p_);
        int len = 1;
        while(len<2*n_)
        {
            len <<= 2;
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
        m.toMontgomery(a_1);
        m.toMontgomery(b_1);
        NTT_radix4(a_1, len, p_, 1, m);
        NTT_radix4(b_1, len, p_, 1, m);
        m.ModMulSIMD(a_1, b_1, c);
        NTT_radix4(c, len, p_, -1, m);
        m.fromMontgomery(c);
        auto End = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 2 * n_ - 1; ++i)
        {
            ab[i] = c[i];
        }
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
