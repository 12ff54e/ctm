#pragma once

#include <cstdint>

#include "../ctm.hpp"

namespace ctm {

namespace detail {

template <typename T>
constexpr T scalbn_impl_u(T x, int n) {
    return n == 0 ? x : scalbn_impl_u(x, n - 1) * T{2};
}

template <typename T>
constexpr T scalbn_impl_d(T x, int n) {
    return n == 0 ? x : scalbn_impl_d(x, n + 1) / T{2};
}

template <typename T>
constexpr int ilogb_impl(T x) {
    return x >= T{2} ? (ilogb_impl(x / T{2}) + 1)
                     : (x < T{1} ? (ilogb_impl(x * T{2}) - 1) : 0);
}

}  // namespace detail

template <typename T>
constexpr numeric_t<T> scalbn(T x, int n) {
    return n > 0 ? detail::scalbn_impl_u(static_cast<numeric_t<T>>(x), n)
                 : detail::scalbn_impl_d(static_cast<numeric_t<T>>(x), n);
}

template <typename T>
constexpr int ilogb(T x) {
    return detail::ilogb_impl(abs(x));
}

namespace detail {

/**
 * @brief Clears the lower 32 bits of a double.
 *
 * @param x
 * @return constexpr double
 */
constexpr double clear_f64_low_bits(double x) {
    return ctm::scalbn(static_cast<double>(static_cast<std::int_least32_t>(
                           ctm::scalbn(x, 20 - ctm::ilogb(x)))),
                       ctm::ilogb(x) - 20);
}

/**
 * @brief Construct a 64 bit float point number from its IEEE754 binary
 * representation
 *
 * @param hi hight 32 bit
 * @param lo low 32 bit
 * @return constexpr double
 */
constexpr double f64_from_bits(std::uint_least32_t hi, std::uint_least32_t lo) {
    // u_int32_t sign = hi >> 31;
    // int32_t exp = static_cast<int32_t>((hi & 0x7ff00000) >> 20) - 1075;
    // int64_t mantissa = (hi & 0xfffff) | 0x100000;
    // mantissa = (sign == 0 ? 1 : -1) * ((mantissa << 32) + lo);
    // return ctm::scalbn(mantissa, exp);

    return ctm::scalbn(
        ((hi >> 31) == 0 ? 1 : -1) *
            ((static_cast<std::int_least64_t>((hi & 0xfffff) | 0x100000)
              << 32) +
             lo),
        static_cast<std::int_least32_t>((hi & 0x7ff00000) >> 20) - 1075);
}

constexpr std::int_least32_t f64_high(double x) {
    // int32_t signbit = x > 0 ? 0 : 1;
    // int32_t exp = ctm::ilogb(x);
    // int32_t mantissa_h = static_cast<int32_t>(ctm::scalbn(abs(x), 20 - exp));
    // return (signbit << 31) | ((exp + 1023) << 20) | (mantissa_h &
    // 0xffefffff);

    return ((x > 0 ? 0 : 1) << 31) | ((ctm::ilogb(x) + 1023) << 20) |
           (static_cast<std::int_least32_t>(
                ctm::scalbn(abs(x), 20 - ctm::ilogb(x))) &
            0xffefffff);
}
constexpr std::uint_least32_t f64_low(double x) {
    return static_cast<std::uint_least32_t>(
        static_cast<std::uint_least64_t>(
            ctm::scalbn(ctm::abs(x), 52 - ctm::ilogb(x))) &
        0xffffffff);
}

}  // namespace detail

}  // namespace ctm
