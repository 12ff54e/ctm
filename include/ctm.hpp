#pragma once

#include "impl/ctm_base.hpp"

// forward declerations

namespace ctm {
namespace detail {
constexpr double cos_impl(double x, double y);
constexpr double sin_impl(double x, double y, int iy);
}  // namespace detail
}  // namespace ctm

// classification functions

#include "impl/Classification.hpp"

// basic functions

#include "impl/abs.hpp"
#include "impl/min-max.hpp"

// float-point manipulation

#include "impl/Float-Point-Manipulation.hpp"

// rounding functions

#include "impl/ceiling.hpp"
#include "impl/floor.hpp"

// power functions

#include "impl/sqrt.hpp"

// ieee

#include "impl/rem_pio2.hpp"

// trigonometry functions

#include "impl/cos.hpp"
#include "impl/sin.hpp"
#include "impl/tan.hpp"
