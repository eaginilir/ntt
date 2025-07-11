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

class Barrett 
{
public:
    Barrett(uint64_t mod): mod_(mod) 
    {
        // 计算 2^64 / mod
        shift_ = 64;
        factor_ = (__uint128_t(1) << shift_) / mod_;
    }

    inline uint64_t reduce(uint64_t x) const 
    {
        return x % mod_;
    }

    inline uint64_t reduce(__uint128_t x) const 
    {
        uint64_t q = (uint64_t)((x * factor_) >> shift_);
        uint64_t r = (uint64_t)x - q * mod_;
        return r >= mod_ ? r - mod_ : r;
    }

    inline uint64_t multiply(uint64_t a, uint64_t b) const 
    {
        return reduce((__uint128_t)a * b);
    }

    inline uint64_t add(uint64_t a, uint64_t b) const 
    {
        uint64_t r = a + b;
        return r >= mod_ ? r - mod_ : r;
    }

    inline uint64_t sub(uint64_t a, uint64_t b) const 
    {
        return a >= b ? a - b : mod_ - (b - a);
    }

private:
    uint64_t mod_;
    int shift_;
    uint64_t factor_;
};

// 更新bit_reverse函数保持不变
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

// 修改NTT实现使用Barrett规约
void NTT_iterative(std::vector<uint64_t> &a, int n, int p, int inv, Barrett &barrett) 
{
    int g = 3; //原根
    bit_reverse(a, n);

    for(int len = 2; len <= n; len <<= 1) 
    {
        int wn = power(g, (p - 1) / len, p);
        if(inv == -1) 
        {
            wn = power(wn, p - 2, p);
        }

        for(int i = 0; i < n; i += len) 
        {
            uint64_t w = 1;
            for(int j = 0; j < len/2; ++j) 
            {
                uint64_t u = a[i + j];
                uint64_t v = barrett.multiply(w, a[i + j + len/2]);
                a[i + j] = barrett.add(u, v);
                a[i + j + len/2] = barrett.sub(u, v);
                w = barrett.multiply(w, wn);
            }
        }
    }

    if(inv == -1) 
    {
        int inv_n = power(n, p - 2, p);
        for(uint64_t &x : a) 
        {
            x = barrett.multiply(x, inv_n);
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
    int test_end = 4;
    for(int i = test_begin; i <= test_end; ++i)
    {
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        Barrett barrett(p_);
        int len = 1;
        while(len < 2*n_) 
        {
            len <<= 1;
        }

        int epochs = 50;
        for (int i = 0; i < epochs; ++i) 
        {
            std::vector<uint64_t> a_1(len, 0);
            std::vector<uint64_t> b_1(len, 0);
            for (int i = 0; i < n_; ++i) 
            {
                a_1[i] = a[i];
                b_1[i] = b[i];
            }
            
            auto Start = std::chrono::high_resolution_clock::now();
            NTT_iterative(a_1, len, p_, 1, barrett);
            NTT_iterative(b_1, len, p_, 1, barrett);
            
            std::vector<uint64_t> c(len, 0);
            for (int i = 0; i < len; ++i) 
            {
                c[i] = barrett.multiply(a_1[i], b_1[i]);
            }
            
            NTT_iterative(c, len, p_, -1, barrett);
            auto End = std::chrono::high_resolution_clock::now();
            
            for (int i = 0; i < 2 * n_ - 1; ++i) 
            {
                ab[i] = c[i];
            }
            std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
            ans += elapsed.count();
        }
        fCheck(ab, n_, i);
        std::cout << "average latency for n = " << n_ << " p = " << p_ << " : " << double(ans / epochs) << " (us) " << std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
    }
    return 0;
}