#pragma once

#include "../ctm.hpp"

namespace ctm {

namespace detail {

template <typename T>
constexpr T abs_impl(T x) {
    return x < T{0} ? -x : x;
}
}  // namespace detail

template <typename T>
constexpr numeric_t<T> abs(T x) {
    return x == T{0} ? numeric_t<T>{0}
                     : detail::abs_impl(static_cast<numeric_t<T>>(x));
}

}  // namespace ctm
