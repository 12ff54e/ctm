#pragma once

namespace ctm {

namespace detail {

template <typename T>
constexpr bool isnan(T x) {
    return x != x;
}

template <typename T>
constexpr bool isinf(T x) {
    return x == limits<T>::infinity() || x == -limits<T>::infinity();
}

}  // namespace detail

}  // namespace ctm
