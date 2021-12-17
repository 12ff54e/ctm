#include "../ctm.hpp"

namespace ctm {

namespace detail {

template <typename T>
constexpr T ceiling_do(T x, T x_int) {
    return x_int + (x < T{0} && x != x_int ? T{0} : T{1});
}

template <typename T>
constexpr T ceiling_impl(T x) {
    return isnan(x) || isinf(x) ? x + x
           : abs(x) > 0x1p+52 || x == T{0}
               ? x
               : ceiling_do(x, static_cast<T>(static_cast<long long int>(x)));
}

}  // namespace detail

template <typename T>
constexpr numeric_t<T> ceiling(T x) {
    return std::is_integral_v<T>
               ? x
               : detail::ceiling_impl(static_cast<numeric_t<T>>(x));
}

}  // namespace ctm
