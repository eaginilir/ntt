#ifndef CPU_UINT128_H
#define CPU_UINT128_H

#include <cstdint>

// CPU端的快速幂函数，使用__uint128_t
uint64_t power(uint64_t base, uint64_t exp, uint64_t mod);

// Montgomery模运算结构体（CPU端实现）
struct MontgomeryModArith {
    uint64_t p, r, r2, p_inv;

    MontgomeryModArith();
    explicit MontgomeryModArith(uint64_t p_);

    static uint64_t mont_get_pinv(uint64_t p);

    uint64_t to_mont(uint64_t a) const;
    uint64_t from_mont(uint64_t a) const;
    uint64_t mont_mul(uint64_t a, uint64_t b) const;
    uint64_t mont_reduce(__uint128_t t) const;

    uint64_t add(uint64_t a, uint64_t b) const;
    uint64_t sub(uint64_t a, uint64_t b) const;
    uint64_t mul(uint64_t a, uint64_t b) const;
};

#endif // CPU_UINT128_H
