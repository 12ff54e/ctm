#pragma once

#include <type_traits>

#include "../ctm.hpp"

namespace ctm {

template <typename T, typename U>
constexpr std::common_type_t<T, U> min(T a, U b) {
    return detail::isnan(a) || detail::isnan(b)
               ? limits<std::common_type_t<T, U>>::quiet_NaN()
           : a > b ? b
                   : a;
}

template <typename T, typename U>
constexpr std::common_type_t<T, U> max(T a, U b) {
    return detail::isnan(a) || detail::isnan(b)
               ? limits<std::common_type_t<T, U>>::quiet_NaN()
           : a < b ? b
                   : a;
}
}  // namespace ctm
