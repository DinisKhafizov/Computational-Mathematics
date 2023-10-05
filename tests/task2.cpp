#include <gtest/gtest.h>
#include <cmath>

const double f(const double x) {
    return x/exp(x * x);
}
const double f_der1(const double x) {
    return (1 - 2 * x * x)/exp(x * x);
}
const double f_der2(const double x) {
    return (4 * x * x * x - 6 * x)/exp(x * x);
}
const double iter_func0(const double x) {
    return (x - f_der1(x)/f_der2(x));
}
const double iter_func1(const double x, const double y_max) {
    return y_max * exp(x * x)/2.;
}
const double iter_func2(const double x, const double y_max) {
    return sqrt(log(2. * x / y_max));
}
const double simple_iteration1(const double (*op)(const double), const double x_0, const double tolerance) {
    double x = x_0, temp = op(x_0);
    while (fabs(x - temp) >= tolerance) {
        x = temp;
        temp = op(x);
    }
    return x;
}
const double simple_iteration2(const double (*op)(const double, const double),const double x0, const double y_max, const double tolerance) {
    double x = x0, temp = op(x0, y_max);
    while (fabs(x - temp) >= tolerance) {
        x = temp;
        temp = op(x, y_max);
    }
    return x;
}

TEST(homework, num2) {
    const double x0 = 0.5, tolerance = 1e-3;
    const double x_max = simple_iteration1(iter_func0, x0, tolerance * 1e-3);
    const double y_max = f(x_max);
    const double x1 = simple_iteration2(iter_func1, x_max - 0.2, y_max, tolerance/2);
    const double x2 = simple_iteration2(iter_func2, x_max + 0.2, y_max, tolerance/2);
    std::cout << std::endl << "x1 = " << x1 << std::endl;
    std::cout << "x2 = " << x2 << std::endl;
    std::cout << "delta = x2 - x1 = " << x2 - x1 << std::endl << std::endl;
}
