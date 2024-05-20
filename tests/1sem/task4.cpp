#include <gtest/gtest.h>

#include <cmath>

const double u1(const double x) { return sqrt(1 - x * x); }
const double u2(const double y) { return atan(y); }
const std::pair<double, double> simple_iteration_sys(const double x0,
                                                     const double y0,
                                                     const double tolerance) {
    double x = x0, y = y0, temp_x = u2(y0), temp_y = u1(x0);
    while (fabs(x - temp_x) >= tolerance && fabs(y - temp_y) >= tolerance) {
        x = temp_x;
        y = temp_y;
        temp_x = u2(temp_y);
        temp_y = u1(temp_x);
    }
    std::pair<double, double> res = {x, y};
    return res;
}

TEST(homework, num4) {
    double x0 = 0.5, y0 = 0.5;
    std::pair<double, double> res = simple_iteration_sys(x0, y0, 0.000001);
    std::cout << std::endl << "x = " << res.first << std::endl;
    std::cout << "y = " << res.second << std::endl;
    std::cout << std::endl << "x = " << -res.first << std::endl;
    std::cout << "y = " << -res.second << std::endl << std::endl;
}