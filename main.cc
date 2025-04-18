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

//如果是逆运算整个结果要除n
class montgomery
{
public:
    montgomery(uint64_t R_, uint64_t N_): N(N_), R(R_) 
    {
        logR = static_cast<int>(std::log2(R));
        uint64_t N_inv = modinv(N, R);
        N_inv_neg = R - N_inv;
        R2 = (__uint128_t(R) * R) % N;
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

    uint64_t ModMul(uint64_t a, uint64_t b) 
    {
        if(a>=N||b>=N)
        {
            throw std::invalid_argument("input integer must be smaller than the modulus N");
        }
        uint64_t aR = REDC(__uint128_t(a) * R2);
        uint64_t bR = REDC(__uint128_t(b) * R2);
        uint64_t abR = REDC(__uint128_t(aR) * bR);
        return REDC(abR);
    }
public:
    uint64_t N, R, R2, N_inv_neg;
    int logR;
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

// 基4的位逆序置换函数
void bit_reverse_radix4(std::vector<uint64_t> &a, int n) 
{
    int log4n = 0;
    int temp = n;
    while (temp > 1) 
    {
        temp >>= 2;
        log4n++;
    }
    std::vector<uint64_t> rev(n, 0);
    for (int i = 0; i < n; ++i) 
    {
        int reversed = 0;
        int num = i;
        for (int j = 0; j < log4n; ++j) 
        {
            reversed = (reversed << 2) | (num & 3);
            num >>= 2;
        }
        rev[i] = reversed;
    }
    for (int i = 0; i < n; ++i) 
    {
        if (i < rev[i]) 
        {
            std::swap(a[i], a[rev[i]]);
        }
    }
}

void NTT_radix4(std::vector<uint64_t> &a,int n, int p, int inv)
{
    uint64_t R = 1ULL << 32;
    montgomery m(R, p);

    int g = 3;
    // bit_reverse(a, n);
    bit_reverse_radix4(a, n);

    for (int len = 4; len <= n;len<<= 2) 
    {
        int step = len >> 2;
        uint64_t w1 = power(g, (p - 1) / len, p);
        if(inv == -1)
        {
            w1 = power(w1, p - 2, p);
        }
        uint64_t w2 = (1ULL * w1 * w1) % p;
        uint64_t w3 = (1ULL * w2 * w1) % p;

        for (int i = 0; i < n;i+=len)
        {
            std::vector<uint64_t> w1_powers(step), w2_powers(step), w3_powers(step);
            w1_powers[0] = w2_powers[0] = w3_powers[0] = 1;
            for (int j = 1; j < step; ++j)
            {
                w1_powers[j] = m.ModMul(w1_powers[j - 1], w1);
                w2_powers[j] = m.ModMul(w2_powers[j - 1], w2);
                w3_powers[j] = m.ModMul(w3_powers[j - 1], w3);

                uint64_t x0 = a[i + j];
                uint64_t x1 = a[i + j + step];
                uint64_t x2 = a[i + j + 2 * step];
                uint64_t x3 = a[i + j + 3 * step];

                // uint64_t t1 = m.ModMul(x1, power(w1, j, p));
                // uint64_t t2 = m.ModMul(x2, power(w2, j, p));
                // uint64_t t3 = m.ModMul(x3, power(w3, j, p));

                uint64_t t1 = m.ModMul(x1, w1_powers[j]);
                uint64_t t2 = m.ModMul(x2, w2_powers[j]);
                uint64_t t3 = m.ModMul(x3, w3_powers[j]);

                uint64_t s0 = (x0 + t2) % p;
                uint64_t s1 = (x0 + p - t2) % p;
                uint64_t s2 = (t1 + p - t3) % p;
                uint64_t s3 = (t1 + t3) % p;

                a[i + j] = (s0 + s3) % p;
                a[i + j + step] = (s1 + s2) % p;
                a[i + j + 2 * step] = (s0 + p - s3) % p;
                a[i + j + 3 * step] = (s1 + p - s2) % p;
            }
        }
    }

    if(inv==-1) 
    {
        int inv_n = power(n, p - 2, p);
        for(auto &x : a) 
        {
            x = m.ModMul(x, inv_n);
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
    int test_end = 1;
    for(int i = test_begin; i <= test_end; ++i){
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        uint64_t R = 1ULL << 32;
        montgomery m(R,p_);
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
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
        // poly_multiply(a, b, ab, n_, p_);
        // NTT_iterative(a_1, len, p_, 1);
        // NTT_iterative(b_1, len, p_, 1);
        NTT_radix4(a_1, len, p_, 1);
        NTT_radix4(b_1, len, p_, 1);
        std::vector<uint64_t> c(len, 0);
        for (int i = 0; i < len; ++i)
        {
            c[i] = m.ModMul(a_1[i], b_1[i]);
        }
        // NTT_iterative(c, len, p_, -1);
        NTT_radix4(c, len, p_, -1);
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
