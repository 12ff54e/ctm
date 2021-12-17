#pragma once

#include <limits>
#include <type_traits>

namespace ctm {

// Math consants

constexpr double PI = 0x1.921fb54442d18p+1;
constexpr double SQRT2 = 0x1.6a09e667f3bcdp+0;
constexpr double SQRT3 = 0x1.bb67ae8584caap+0;

/**
 * @brief Casting integral type to double, keeping float-point type unchanged
 *
 * @tparam T
 */
template <typename T>
using numeric_t = std::conditional_t<std::is_integral_v<T>, double, T>;

template <typename T>
using limits = std::numeric_limits<T>;

}  // namespace ctm
