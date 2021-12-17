#include "../ctm.hpp"

namespace ctm {

namespace detail {

template <typename T>
constexpr T floor_do(T x, T x_int) {
    return x_int - (x < T{0} && x != x_int ? T{1} : T{0});
}

template <typename T>
constexpr T floor_impl(T x) {
    return isnan(x) || isinf(x) ? x + x
           : abs(x) > 0x1p+52 || x == T{0}
               ? x
               : floor_do(x, static_cast<T>(static_cast<long long int>(x)));
}

}  // namespace detail

template <typename T>
constexpr numeric_t<T> floor(T x) {
    return std::is_integral_v<T>
               ? x
               : detail::floor_impl(static_cast<numeric_t<T>>(x));
}

}  // namespace ctm
