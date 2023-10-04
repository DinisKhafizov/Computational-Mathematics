#include "Quadrature.hpp"

const double trapezoid(const double a, const double b, const double (*y)(const double)) {
    return (b - a) * (y(a) + y(b)) / 2;
}

const double simpson(const double a, const double b, const double (*y)(const double)) {
    return (b - a) * (y(a) + 4 * y((a + b)/2) + y(b)) / 6;
}

const double three_eights(const double a, const double b, const double (*y)(const double)) {
    return (b - a) * (y(a) + 3 * y((2 * a + b)/3) + 3 * y((a + 2 * b)/3) + y(b)) / 8;
}
const double gauss_quadrature(const double a, const double b, const size_t deg, const double (*y)(const double)) {
    const double k1 = (a + b)/2., k2 = (b - a)/2.;
    double res;
    if (deg == 2) {
        res = (y(k1 + k2 / sqrt(3)) + y(k1 - k2 / sqrt(3)));
    }
    if (deg == 3) {
        res = (5./9 * y(k1 + k2 * sqrt(3./5)) + 8./9 * y(k1) + 5./9 * y(k1 - k2 * sqrt(3./5)));
    }
    if (deg == 4) {
        res = (0.347855 * y(k1 + k2 * 0.861136) + 0.652145 * y(k1 + k2 * 0.339981) + 0.652145 * y(k1 - k2 * 0.339981) + 0.347855 * y(k1 - k2 * 0.861136));
    }
    return res * k2;
}