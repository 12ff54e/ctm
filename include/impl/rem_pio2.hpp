#pragma once

#include "../ctm.hpp"

namespace ctm {

namespace detail {

constexpr int rem_pio2_large_impl(double* x,
                                  double* y,
                                  int e0,
                                  int nx,
                                  int prec) {
    constexpr int TWO_OVER_PI_BITS[] = {
        0xA2F983, 0x6E4E44, 0x1529FC, 0x2757D1, 0xF534DD, 0xC0DB62, 0x95993C,
        0x439041, 0xFE5163, 0xABDEBB, 0xC561B7, 0x246E3A, 0x424DD2, 0xE00649,
        0x2EEA09, 0xD1921C, 0xFE1DEB, 0x1CB129, 0xA73EE8, 0x8235F5, 0x2EBB44,
        0x84E99C, 0x7026B4, 0x5F7E41, 0x3991D6, 0x398353, 0x39F49C, 0x845F8B,
        0xBDF928, 0x3B1FF8, 0x97FFDE, 0x05980F, 0xEF2F11, 0x8B5A0A, 0x6D1F6D,
        0x367ECF, 0x27CB09, 0xB74F46, 0x3F669E, 0x5FEA2D, 0x7527BA, 0xC7EBE5,
        0xF17B3D, 0x0739F7, 0x8A5292, 0xEA6BFB, 0x5FB11F, 0x8D5D08, 0x560330,
        0x46FC7B, 0x6BABF0, 0xCFBC20, 0x9AF436, 0x1DA9E3, 0x91615E, 0xE61B08,
        0x659985, 0x5F14A0, 0x68408D, 0xFFD880, 0x4D7327, 0x310606, 0x1556CA,
        0x73A8C9, 0x60E27B, 0xC08C6B};

    constexpr double PI_OVER_TWO[] = {
        1.57079625129699707031e+00, 7.54978941586159635335e-08,
        5.39030252995776476554e-15, 3.28200341580791294123e-22,
        1.27065575308067607349e-29, 1.22933308981111328932e-36,
        2.73370053816464559624e-44, 2.16741683877804819444e-51,
    };

    constexpr double ZERO = 0.;
    constexpr double ONE = 1.;
    constexpr double TWO_24 = 0x1p+24;
    constexpr double TWO_NEG_24 = 0x1p-24;

    constexpr int init_jk[] = {2, 3, 4, 6};

    int jk = init_jk[prec];
    int jp = jk;

    int jx = nx - 1;
    int jv = (e0 - 3) / 24;
    if (jv < 0) { jv = 0; }
    int q0 = e0 - 24 * (jv + 1);

    int j = jv - jx;
    int m = jx + jk;
    double f[20] = {0.};
    for (int i = 0; i <= m; i++, j++) {
        f[i] = (j < 0) ? ZERO : static_cast<double>(TWO_OVER_PI_BITS[j]);
    }

    double q[20] = {0.};
    for (int i = 0; i <= jk; i++) {
        for (int j = 0, fw = 0.0; j <= jx; j++) {
            fw += x[j] * f[jx + i - j];
            q[i] = fw;
        }
    }

    int jz = jk;

    double z = q[jz];
    double fw = 0.;
    int iq[20] = {0};
    int n = 0;
    int ih = 0;

    do {
        z = q[jz];

        for (int i = 0, j = jz; j > 0; i++, j--) {
            fw = static_cast<double>(static_cast<int>(0x1p-24 * z));
            iq[i] = static_cast<int>(z - 0x1p24 * fw);
            z = q[j - 1] + fw;
        }

        z = -8. * floor(.125 * z * (2 << q0));
        n = static_cast<int>(z);
        z -= static_cast<double>(n);
        ih = 0;
        if (q0 > 0) {
            int i = (iq[jz - 1] >> (24 - q0));
            n += i;
            iq[jz - 1] -= i << (24 - q0);
            ih = iq[jz - 1] >> (23 - q0);
        } else if (q0 == 0) {
            ih = iq[jz - 1] >> 23;
        } else if (z >= 0.5) {
            ih = 2;
        }

        if (ih > 0) { /* q > 0.5 */
            n += 1;
            int carry = 0;
            for (int i = 0; i < jz; i++) { /* compute 1-q */
                j = iq[i];
                if (carry == 0) {
                    if (j != 0) {
                        carry = 1;
                        iq[i] = 0x1000000 - j;
                    }
                } else
                    iq[i] = 0xffffff - j;
            }
            if (q0 > 0) { /* rare case: chance is 1 in 12 */
                switch (q0) {
                    case 1:
                        iq[jz - 1] &= 0x7fffff;
                        break;
                    case 2:
                        iq[jz - 1] &= 0x3fffff;
                        break;
                }
            }
            if (ih == 2) {
                z = ONE - z;
                if (carry != 0) z -= scalbn(ONE, q0);
            }
        }

        /* check if recomputation is needed */
        if (z == ZERO) {
            j = 0;
            for (int i = jz - 1; i >= jk; i--) { j |= iq[i]; }

            if (j == 0) { /* need recomputation */
                int k = 1;
                for (; iq[jk - k] == 0; k++)
                    ;
                /* k = no. of terms needed */

                for (int i = jz + 1; i <= jz + k; i++) {
                    /* add q[jz+1] to q[jz+k] */
                    f[jx + i] = static_cast<double>(TWO_OVER_PI_BITS[jv + i]);
                    for (j = 0, fw = 0.0; j <= jx; j++) {
                        fw += x[j] * f[jx + i - j];
                    }
                    q[i] = fw;
                }
                jz += k;
                continue;
            }
        }
        break;
    } while (true);

    /* chop off zero terms */
    if (z == 0.0) {
        jz -= 1;
        q0 -= 24;
        while (iq[jz] == 0) {
            jz--;
            q0 -= 24;
        }
    } else { /* break z into 24-bit if necessary */
        z = scalbn(z, -q0);
        if (z >= TWO_24) {
            fw = static_cast<double>(static_cast<int>(TWO_NEG_24 * z));
            iq[jz] = static_cast<int>(z - TWO_24 * fw);
            jz += 1;
            q0 += 24;
            iq[jz] = static_cast<int>(fw);
        } else
            iq[jz] = static_cast<int>(z);
    }

    /* convert integer "bit" chunk to floating-point value */
    fw = scalbn(ONE, q0);
    for (int i = jz; i >= 0; i--) {
        q[i] = fw * static_cast<double>(iq[i]);
        fw *= TWO_NEG_24;
    }

    /* compute PIo2[0,...,jp]*q[jz,...,0] */
    double fq[20] = {ZERO};
    for (int i = jz; i >= 0; i--) {
        int k = 0;
        for (fw = 0.0; k <= jp && k <= jz - i; k++)
            fw += PI_OVER_TWO[k] * q[i + k];
        fq[jz - i] = fw;
    }

    /* compress fq[] into y[] */
    switch (prec) {
        case 0:
            fw = 0.0;
            for (int i = jz; i >= 0; i--) fw += fq[i];
            y[0] = (ih == 0) ? fw : -fw;
            break;
        case 1:
        case 2:
            fw = 0.0;
            for (int i = jz; i >= 0; i--) fw += fq[i];
            y[0] = (ih == 0) ? fw : -fw;
            fw = fq[0] - fw;
            for (int i = 1; i <= jz; i++) fw += fq[i];
            y[1] = (ih == 0) ? fw : -fw;
            break;
        case 3: /* painful */
            for (int i = jz; i > 0; i--) {
                fw = fq[i - 1] + fq[i];
                fq[i] += fq[i - 1] - fw;
                fq[i - 1] = fw;
            }
            for (int i = jz; i > 1; i--) {
                fw = fq[i - 1] + fq[i];
                fq[i] += fq[i - 1] - fw;
                fq[i - 1] = fw;
            }
            fw = 0.0;
            for (int i = jz; i >= 2; i--) fw += fq[i];
            if (ih == 0) {
                y[0] = fq[0];
                y[1] = fq[1];
                y[2] = fw;
            } else {
                y[0] = -fq[0];
                y[1] = -fq[1];
                y[2] = -fw;
            }
    }
    return n & 7;
}

}  // namespace detail

constexpr int rem_pio2(double x, double* y) {
    constexpr double N_PI_OVER_2_D[] = {
        0x1.921FBp0,  0x1.921FBp+1, 0x1.2D97Cp+2, 0x1.921FBp+2, 0x1.F6A7Ap+2,
        0x1.2D97Cp+3, 0x1.5FDBBp+3, 0x1.921FBp+3, 0x1.C463Ap+3, 0x1.F6A7Ap+3,
        0x1.1475Cp+4, 0x1.2D97Cp+4, 0x1.46B9Cp+4, 0x1.5FDBBp+4, 0x1.78FDBp+4,
        0x1.921FBp+4, 0x1.AB41Bp+4, 0x1.C463Ap+4, 0x1.DD85Ap+4, 0x1.F6A7Ap+4,
        0x1.07E4Cp+5, 0x1.1475Cp+5, 0x1.2106Cp+5, 0x1.2D97Cp+5, 0x1.3A28Cp+5,
        0x1.46B9Cp+5, 0x1.534ACp+5, 0x1.5FDBBp+5, 0x1.6C6CBp+5, 0x1.78FDBp+5,
        0x1.858EBp+5, 0x1.921FBp+5};

    constexpr double N_PI_OVER_2_U[] = {
        0x1.921FCp0,  0x1.921FCp+1, 0x1.2D97Dp+2, 0x1.921FCp+2, 0x1.F6A7Bp+2,
        0x1.2D97Dp+3, 0x1.5FDBCp+3, 0x1.921FCp+3, 0x1.C463Bp+3, 0x1.F6A7Bp+3,
        0x1.1475Dp+4, 0x1.2D97Dp+4, 0x1.46B9Dp+4, 0x1.5FDBCp+4, 0x1.78FDCp+4,
        0x1.921FCp+4, 0x1.AB41Cp+4, 0x1.C463Bp+4, 0x1.DD85Bp+4, 0x1.F6A7Bp+4,
        0x1.07E4Dp+5, 0x1.1475Dp+5, 0x1.2106Dp+5, 0x1.2D97Dp+5, 0x1.3A28Dp+5,
        0x1.46B9Dp+5, 0x1.534ADp+5, 0x1.5FDBCp+5, 0x1.6C6CCp+5, 0x1.78FDCp+5,
        0x1.858ECp+5, 0x1.921FCp+5};

    constexpr double ZERO = 0., HALF = 0x1p-1, TWO_24 = 0x1p+24,
                     TWO_OVER_PI = 6.36619772367581382433e-01,
                     PI_OVER_TWO_1 = 1.57079632673412561417e+00,
                     PI_OVER_TWO_1t = 6.07710050650619224932e-11,
                     PI_OVER_TWO_2 = 6.07710050630396597660e-11,
                     PI_OVER_TWO_2t = 2.02226624879595063154e-21,
                     PI_OVER_TWO_3 = 2.02226624871116645580e-21,
                     PI_OVER_TWO_3t = 8.47842766036889956997e-32;

    if (detail::isnan(x) || detail::isinf(x)) {
        y[0] = y[1] = x - x;
        return 0;
    }

    const double x_abs = abs(x);

    if (x_abs <= 0x1.921fbp-1) {
        /* |x| ~<= pi/4 , no need for reduction */
        y[0] = x;
        y[1] = 0;
        return 0;
    }

    if (x_abs < 0x1.2d97cp+1) {
        /* |x| < 3pi/4, special case with n=+-1 */
        if (x > 0) {
            double z = x - PI_OVER_TWO_1;
            if (x_abs < 0x1.921fbp0 || x_abs >= 0x1.921fcp0) {
                /* 33+53 bit pi is good enough */
                y[0] = z - PI_OVER_TWO_1t;
                y[1] = (z - y[0]) - PI_OVER_TWO_1t;
            } else { /* near pi/2, use 33+33+53 bit pi */
                z -= PI_OVER_TWO_2;
                y[0] = z - PI_OVER_TWO_2t;
                y[1] = (z - y[0]) - PI_OVER_TWO_2t;
            }
            return 1;
        } else { /* negative x */
            double z = x + PI_OVER_TWO_1;
            if (x_abs < 0x1.921fbp0 || x_abs >= 0x1.921fcp0) {
                /* 33+53 bit pi is good enough */
                y[0] = z + PI_OVER_TWO_1t;
                y[1] = (z - y[0]) + PI_OVER_TWO_1t;
            } else { /* near pi/2, use 33+33+53 bit pi */
                z += PI_OVER_TWO_2;
                y[0] = z + PI_OVER_TWO_2t;
                y[1] = (z - y[0]) + PI_OVER_TWO_2t;
            }
            return -1;
        }
    }

    if (x_abs <= 0x1.921fbp+20) { /* |x| ~<= 2^19*(pi/2), medium size */
        int n = static_cast<int>(x_abs * TWO_OVER_PI + HALF);
        double fn = static_cast<double>(n);
        double r = x_abs - fn * PI_OVER_TWO_1;
        double w = fn * PI_OVER_TWO_1t; /* 1st round good to 85 bit */
        if (n < 32 &&
            (x_abs < N_PI_OVER_2_D[n - 1] || x_abs >= N_PI_OVER_2_U[n - 1])) {
            y[0] = r - w; /* quick check no cancellation */
        } else {
            int j = ilogb(x);
            y[0] = r - w;
            int i = j - (ilogb(y[0]) & 0x7ff);

            if (i > 16) {
                double t = r;
                double w = fn * PI_OVER_TWO_2;
                r = t - w;
                w = fn * PI_OVER_TWO_2t - ((t - r) - w);
                y[0] = r - w;
                i = j - (ilogb(y[0]) & 0x7ff);
                if (i > 49) { /* 3rd iteration need, 151 bits acc */
                    t = r;    /* will cover all possible cases */
                    w = fn * PI_OVER_TWO_3;
                    r = t - w;
                    w = fn * PI_OVER_TWO_3t - ((t - r) - w);
                    y[0] = r - w;
                }
            }
        }
        y[1] = (r - y[0]) - w;
        if (x < 0) {
            y[0] = -y[0];
            y[1] = -y[1];
            return -n;
        } else {
            return n;
        }
    }

    /*
     * all other (large) arguments
     */

    int e0 = ilogb(x) - 23;
    double z = scalbn(x_abs, -e0);
    double tx[3] = {0.};
    for (int i = 0; i < 2; i++) {
        tx[i] = static_cast<double>(static_cast<int>(z));
        z = (z - tx[i]) * TWO_24;
    }
    tx[2] = z;

    int nx = 3;
    while (tx[nx - 1] == ZERO) { nx--; } /* skip zero term */
    int n = detail::rem_pio2_large_impl(tx, y, e0, nx, 2);
    if (x < 0) {
        y[0] = -y[0];
        y[1] = -y[1];
        return -n;
    }

    return n;
}

}  // namespace ctm
