#include "../ctm.hpp"

namespace ctm {

namespace detail {

constexpr double tan_impl(double x, double y, int iy) {
    constexpr double COEF[] = {
        3.33333333333334091986e-01,  /* 3FD55555, 55555563 */
        1.33333333333201242699e-01,  /* 3FC11111, 1110FE7A */
        5.39682539762260521377e-02,  /* 3FABA1BA, 1BB341FE */
        2.18694882948595424599e-02,  /* 3F9664F4, 8406D637 */
        8.86323982359930005737e-03,  /* 3F8226E3, E96E8493 */
        3.59207910759131235356e-03,  /* 3F6D6D22, C9560328 */
        1.45620945432529025516e-03,  /* 3F57DBC8, FEE08315 */
        5.88041240820264096874e-04,  /* 3F4344D8, F2F26501 */
        2.46463134818469906812e-04,  /* 3F3026F7, 1A8D1068 */
        7.81794442939557092300e-05,  /* 3F147E88, A03792A6 */
        7.14072491382608190305e-05,  /* 3F12B80F, 32F0A7E9 */
        -1.85586374855275456654e-05, /* BEF375CB, DB605373 */
        2.59073051863633712884e-05,  /* 3EFB2A70, 74BF7AD4 */
    };

    constexpr double ONE = 0x1p0;
    constexpr double PI_OVER_4 = 0x1.921fb54442d18p-1;
    constexpr double PI_OVER_4_L = 0x1.1a62633145c07p-55;

    const double x_abs = abs(x);
    if (x_abs < 0x1p-28) {
        if (x_abs == 0 && iy + 1 == 0) {
            return ONE / x_abs;
        } else {
            if (iy == 1) {
                return x;
            } else {
                double z = x + y;
                double w = z;
                z = clear_f64_low_bits(z);
                double t = -ONE / (x + y);
                double a = t;
                t = clear_f64_low_bits(t);
                return t + a * (ONE + t * z + t * (y - (z - x)));
            }
        }
    }

    double z = PI_OVER_4 - x;
    double w = PI_OVER_4_L - y;
    if (x_abs >= 0x1.59428p-1) {
        if (x < 0) {
            x = -x;
            y = -y;
        }
        x = z + w;
        y = 0.;
    }
    z = x * x;
    w = z * z;

    double r =
        COEF[1] +
        w * (COEF[3] +
             w * (COEF[5] + w * (COEF[7] + w * (COEF[9] + w * COEF[11]))));
    double v =
        z *
        (COEF[2] +
         w * (COEF[4] +
              w * (COEF[6] + w * (COEF[8] + w * (COEF[10] + w * COEF[12])))));
    double s = z * x;
    r = y + z * (s * (r + v) + y);
    r += COEF[0] * s;
    w = x + r;

    if (x_abs >= 0x1.59428p-1) {
        v = static_cast<double>(iy);
        return static_cast<double>(x > 0 ? 1 : -1) *
               (v - 2.0 * (x - (w * w / (w + v) - r)));
    }

    if (iy == 1) {
        return w;
    } else {
        z = w;
        z = clear_f64_low_bits(z);
        double t = -1. / w;
        double a = t;
        t = clear_f64_low_bits(t);
        return t + a * (ONE + t * z + t * (r - (z - x)));
    }
}

}  // namespace detail

constexpr double tan(double x) {
    if (detail::isnan(x) || detail::isinf(x)) { return x - x; }

    if (abs(x) < 0x1.921fcp-1) {
        return detail::tan_impl(x, 0., 1);
    } else {
        double y[2] = {0.};
        int n = rem_pio2(x, y);
        return detail::tan_impl(y[0], y[1], 1 - ((n & 1) << 1));
    }
}

}  // namespace ctm
