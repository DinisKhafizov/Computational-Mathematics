#include <gtest/gtest.h>
#include <cmath>

const double y_1(const double x) {
    return x/exp(x * x);
}

const double y_der(const double x) {
    return 1/sqrt(2);
}

const double y_2(const double y_max, const double x) {
    return exp(x * x) * y_max / 2;
}

const double y_3(const double y_max, const double x) {
    return sqrt(log(2 * x / y_max));
}

const double simple_iteration1(const double x_0, const double tolerance) {
    double x = x_0, temp = y_der(x_0);
    while (fabs(x - temp) >= tolerance) {
        x = temp;
        temp = y_der(x);
    }
    return x;
}

const double simple_iteration2(const double (*op)(const double, const double),const double x_max, const double y_max, const double tolerance, const double shift) {
    double x = x_max + shift, temp = op(y_max, x_max + shift);
    while (fabs(x - temp) >= tolerance) {
        x = temp;
        temp = op(y_max, x);
    }
    return x;
}

TEST(homework, num2) {
    double x = 0, tolerance = 0.001;
    x = simple_iteration1(x, tolerance);
    double f_max = y_1(x);
    double x1 = simple_iteration2(y_2, x, f_max, tolerance/2, -0.2);
    double x2 = simple_iteration2(y_3, x, f_max, tolerance/2, 0.2);
    std::cout << std::endl << "x1 = " << x1 << std::endl;
    std::cout << "x2 = " << x2 << std::endl;
    std::cout << "delta = x2 - x1 = " << x2 - x1 << std::endl << std::endl;
}