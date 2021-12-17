#pragma once

#include "../ctm.hpp"

namespace ctm {

namespace detail {

template <typename T>
constexpr T sqrt_newton_iter(T x, T xn) {
    return abs(xn - x / xn) / x < limits<T>::epsilon()
               ? xn
               : sqrt_newton_iter(x, T{0.5} * (xn + x / xn));
}

template <typename T>
constexpr T sqrt_impl(T x, T f) {
    return x < T{0} ? limits<T>::signaling_NaN()
           : detail::isnan(x) || detail::isinf(x) || x == T{0} ? x
           : limits<T>::epsilon() > abs(T{1} - x)              ? f
           : x < T{1} ? sqrt_impl(T{4} * x, f * T{0.5})
           : x > T{4} ? sqrt_impl(x / T{4}, f * T{2})
                      : sqrt_newton_iter(x, T{0x1.4dbf86a314dc0p-2} * x +
                                                T{0x1.8e38e38e38e39p-1}) *
                            f;
}

/**
 * @brief This ieee754 conforming sqrt rounds correctly.
 *
 * @param x
 * @return constexpr double
 */
constexpr double __ieee754_sqrt(double x) {
    if (isnan(x) || isinf(x)) { return x * x + x; }

    if (x == 0) {
        return x; /* $sqrt(\pm 0) = \pm 0$ */
    } else if (x < 0) {
        return (x - x) / (x - x);
    }

    constexpr int sign = (int)0x80000000;
    constexpr double ONE = 1.;
    constexpr double TINY = 1.e-300;

    int ix0 = f64_high(x);
    unsigned ix1 = f64_low(x);
    int m = ix0 >> 20;
    if (m == 0) { /* subnormal x */
        while (ix0 == 0) {
            m -= 21;
            ix0 |= (ix1 >> 11);
            ix1 <<= 21;
        }
        int i = 0;
        for (i = 0; (ix0 & 0x00100000) == 0; i++) ix0 <<= 1;
        m -= i - 1;
        ix0 |= (ix1 >> (32 - i));
        ix1 <<= i;
    }
    m -= 1023;
    ix0 = (ix0 & 0x000fffff) | 0x00100000;
    if (m & 1) { /* odd m, double x to make it even */
        ix0 += ix0 + ((ix1 & sign) >> 31);
        ix1 += ix1;
    }
    m >>= 1;

    /* generate sqrt(x) bit by bit */
    ix0 += ix0 + ((ix1 & sign) >> 31);
    ix1 += ix1;

    int q{}, s0{}, t{};
    unsigned q1{}, s1{}, t1{};
    q = q1 = s0 = s1 = 0;    /* [q,q1] = sqrt(x) */
    unsigned r = 0x00200000; /* r = moving bit from right to left */

    while (r != 0) {
        t = s0 + r;
        if (t <= ix0) {
            s0 = t + r;
            ix0 -= t;
            q += r;
        }
        ix0 += ix0 + ((ix1 & sign) >> 31);
        ix1 += ix1;
        r >>= 1;
    }

    r = sign;
    while (r != 0) {
        t1 = s1 + r;
        t = s0;
        if ((t < ix0) || ((t == ix0) && (t1 <= ix1))) {
            s1 = t1 + r;
            if (((t1 & sign) == sign) && (s1 & sign) == 0) s0 += 1;
            ix0 -= t;
            if (ix1 < t1) ix0 -= 1;
            ix1 -= t1;
            q1 += r;
        }
        ix0 += ix0 + ((ix1 & sign) >> 31);
        ix1 += ix1;
        r >>= 1;
    }

    /* use floating add to find out rounding direction */
    if ((ix0 | ix1) != 0) {
        double z = ONE - TINY; /* trigger inexact flag */
        if (z >= ONE) {
            z = ONE + TINY;
            if (q1 == (unsigned)0xffffffff) {
                q1 = 0;
                q += 1;
            } else if (z > ONE) {
                if (q1 == (unsigned)0xfffffffe) q += 1;
                q1 += 2;
            } else
                q1 += (q1 & 1);
        }
    }
    ix0 = (q >> 1) + 0x3fe00000;
    ix1 = q1 >> 1;
    if ((q & 1) == 1) ix1 |= sign;
    ix0 += (m << 20);

    return f64_from_bits(ix0, ix1);
}

}  // namespace detail

/**
 * @brief Square root function, rounding to zero.
 *
 * @tparam T
 * @param x
 * @return constexpr numeric_t<T>
 */
template <typename T>
constexpr numeric_t<T> sqrt(T x) {
    // return detail::sqrt_impl(static_cast<numeric_t<T>>(x), numeric_t<T>{1});
    return detail::__ieee754_sqrt(static_cast<double>(x));
}

}  // namespace ctm