#pragma once

#include "../ctm.hpp"

namespace ctm {

namespace detail {

constexpr double sin_impl(double x, double y, int iy) {
    constexpr double
        S1 = -1.66666666666666324348e-01, /* 0xBFC55555, 0x55555549 */
        S2 = 8.33333333332248946124e-03,  /* 0x3F811111, 0x1110F8A6 */
        S3 = -1.98412698298579493134e-04, /* 0xBF2A01A0, 0x19C161D5 */
        S4 = 2.75573137070700676789e-06,  /* 0x3EC71DE3, 0x57B1FE7D */
        S5 = -2.50507602534068634195e-08, /* 0xBE5AE5E6, 0x8A2B9CEB */
        S6 = 1.58969099521155010221e-10;  /* 0x3DE5D93A, 0x5ACFD57C */

    if (abs(x) < 0x1p-27) { return x; }

    double z = x * x;
    double v = z * x;
    double r = S2 + z * (S3 + z * (S4 + z * (S5 + z * S6)));

    if (iy == 0) {
        return x + v * (S1 + z * r);
    } else {
        return x - ((z * (.5 * y - v * r) - y) - v * S1);
    }
}

}  // namespace detail

constexpr double sin(double x) {
    if (detail::isnan(x) || detail::isinf(x)) { return x - x; }

    if (abs(x) < 0x1.921fcp-1) { return detail::sin_impl(x, 0, 0); }

    double y[2] = {0};
    int n = rem_pio2(x, y);
    switch (n & 3) {
        case 0:
            return detail::sin_impl(y[0], y[1], 1);
        case 1:
            return detail::cos_impl(y[0], y[1]);
        case 2:
            return -detail::sin_impl(y[0], y[1], 1);
        default:
            return -detail::cos_impl(y[0], y[1]);
    }
}

}  // namespace ctm
