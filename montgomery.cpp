#include "montgomery.h"

// CPU版本快速幂
uint64_t power(uint64_t base, uint64_t exp, uint64_t mod) {
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

// 构造函数
MontgomeryModArith::MontgomeryModArith() : p(0), r(0), r2(0), p_inv(0) {}

MontgomeryModArith::MontgomeryModArith(uint64_t p_) : p(p_) {
    r = (1ULL << 63) % p; 
    r = (r << 1) % p; // r = 2^64 mod p
    r2 = ((__uint128_t)r * r) % p;
    p_inv = mont_get_pinv(p);
}

uint64_t MontgomeryModArith::mont_get_pinv(uint64_t p) {
    uint64_t ret = 1;
    for (int i = 0; i < 6; ++i) ret *= 2 - p * ret;
    return ~ret + 1;
}

uint64_t MontgomeryModArith::to_mont(uint64_t a) const {
    return mont_mul(a, r2);
}

uint64_t MontgomeryModArith::from_mont(uint64_t a) const {
    return mont_reduce(a);
}

uint64_t MontgomeryModArith::mont_mul(uint64_t a, uint64_t b) const {
    return mont_reduce((__uint128_t)a * b);
}

uint64_t MontgomeryModArith::mont_reduce(__uint128_t t) const {
    uint64_t m = uint64_t(t) * p_inv;
    __uint128_t u = t + (__uint128_t)m * p;
    u >>= 64;
    if (u >= p) u -= p;
    return (uint64_t)u;
}

uint64_t MontgomeryModArith::add(uint64_t a, uint64_t b) const {
    uint64_t res = a + b;
    return res >= p ? res - p : res;
}

uint64_t MontgomeryModArith::sub(uint64_t a, uint64_t b) const {
    return a >= b ? a - b : a + p - b;
}

uint64_t MontgomeryModArith::mul(uint64_t a, uint64_t b) const {
    return mont_mul(a, b);
}
