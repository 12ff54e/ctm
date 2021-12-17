#pragma once

#include "../ctm.hpp"

namespace ctm {

namespace detail {

constexpr double cos_impl(double x, double y) {
    constexpr double
        ONE = 1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
        C1 = 4.16666666666666019037e-02,  /* 0x3FA55555, 0x5555554C */
        C2 = -1.38888888888741095749e-03, /* 0xBF56C16C, 0x16C15177 */
        C3 = 2.48015872894767294178e-05,  /* 0x3EFA01A0, 0x19CB1590 */
        C4 = -2.75573143513906633035e-07, /* 0xBE927E4F, 0x809C52AD */
        C5 = 2.08757232129817482790e-09,  /* 0x3E21EE9E, 0xBDB4B1C4 */
        C6 = -1.13596475577881948265e-11; /* 0xBDA8FAE9, 0xBE8838D4 */

    if (abs(x) < 0x1p-27) { return ONE; }

    double z = x * x;
    double r = z * (C1 + z * (C2 + z * (C3 + z * (C4 + z * (C5 + z * C6)))));
    if (abs(x) < 0x1.33333p-2) {
        return ONE - (.5 * z - (z * r - x * y));
    } else {
        double qx = abs(x) > 0x1.9p-1 ? 0x1.2p-2 : .25 * x;
        qx = clear_f64_low_bits(qx);
        return ONE - qx - (.5 * z - qx - (z * r - x * y));
    }
}
}  // namespace detail

constexpr double cos(double x) {
    if (detail::isnan(x) || detail::isinf(x)) { return x - x; }

    if (abs(x) < 0x1.921fcp-1) { return detail::cos_impl(x, 0); }

    double y[2] = {0};
    int n = rem_pio2(x, y);
    switch (n & 3) {
        case 0:
            return detail::cos_impl(y[0], y[1]);
        case 1:
            return -detail::sin_impl(y[0], y[1], 1);
        case 2:
            return -detail::cos_impl(y[0], y[1]);
        default:
            return detail::sin_impl(y[0], y[1], 1);
    }
}

}  // namespace ctm
