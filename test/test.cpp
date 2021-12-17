// #include <iostream>
#include "../ctm.hpp"

constexpr bool approx_eq(double a, double b) {
    return ctm::abs(a - b) <=
           ctm::limits<double>::epsilon() * ctm::min(ctm::abs(a), ctm::abs(b));
}

int main(int argc, char const* argv[]) {
    // basic

    static_assert(ctm::abs(-2.3) == 2.3);
    static_assert(ctm::abs(4.9) == 4.9);

    static_assert(ctm::min(1.4, 2.3) == 1.4);

    // float point

    static_assert(ctm::ilogb(12.5) == 3);
    static_assert(ctm::scalbn(2.5, 3) == 20.);

    static_assert(ctm::detail::f64_from_bits(ctm::detail::f64_high(ctm::PI),
                                             ctm::detail::f64_low(ctm::PI)) ==
                  ctm::PI);
    static_assert(ctm::detail::f64_from_bits(
                      ctm::detail::f64_high(-ctm::SQRT2),
                      ctm::detail::f64_low(-ctm::SQRT2)) == -ctm::SQRT2);

    // floor test

    static_assert(ctm::floor(573.2) == 573.);
    static_assert(ctm::floor(-156.98) == -157.);
    static_assert(ctm::floor(std::numeric_limits<double>::infinity()) ==
                  std::numeric_limits<double>::infinity());

    // ceiling test

    static_assert(ctm::ceiling(573.2) == 574.);
    static_assert(ctm::ceiling(-156.98) == -156.);
    static_assert(ctm::ceiling(std::numeric_limits<double>::infinity()) ==
                  std::numeric_limits<double>::infinity());

    // sqrt

    static_assert(ctm::sqrt(2.) == ctm::SQRT2);
    static_assert(ctm::sqrt(3.) == ctm::SQRT3);
    static_assert(ctm::sqrt(9801) == 99.);
    static_assert(ctm::sqrt(-0.) == -0.);
    static_assert(ctm::sqrt(ctm::limits<double>::infinity()) ==
                  ctm::limits<double>::infinity());

    // trig

    static_assert(approx_eq(ctm::sin(ctm::PI / 6), .5));
    static_assert(approx_eq(ctm::sin(ctm::PI / 4), .5 * ctm::SQRT2));
    static_assert(approx_eq(ctm::sin(ctm::PI / 3), .5 * ctm::SQRT3));

    static_assert(approx_eq(ctm::cos(ctm::PI / 6), .5 * ctm::SQRT3));
    static_assert(approx_eq(ctm::cos(ctm::PI / 4), .5 * ctm::SQRT2));
    static_assert(approx_eq(ctm::cos(ctm::PI / 3), .5));

    static_assert(approx_eq(ctm::tan(ctm::PI / 6), ctm::SQRT3 / 3));
    static_assert(approx_eq(ctm::tan(ctm::PI / 4), 1.));
    static_assert(
        approx_eq(ctm::tan(ctm::PI / 3 + ctm::limits<double>::epsilon()),
                  ctm::SQRT3 + ctm::limits<double>::epsilon()));

    return 0;
}
