#include <gtest/gtest.h>
#include "../src/Interpolation/Cubic_spline.hpp"
#include "../src/Interpolation/Newton.hpp"
#include "../src/Interpolation/Lagrange.hpp"
TEST(newton, x) { //x
    std::vector<double> x = {0, 1, 2, 3, 4, 5}, y = {0, 1, 2, 3, 4, 5};
    std::vector<double> res(6);
    res = newton_interpolation(x, y);
    get_and_write_Yvalues("Newton_x.txt", x, y, res, 100);
}
TEST(newton, x2) { // x**2
    std::vector<double> x = {0, 1, 2, 3, 4, 5}, y = {0, 1, 4, 9, 16, 25};
    std::vector<double> res(6);
    res = newton_interpolation(x, y);
    get_and_write_Yvalues("Newton_x2.txt", x, y, res, 100);
}
TEST(newton, x2_x_sqx) { // x**2 - x - x^(1/2)
    std::vector<double> x = {3, 7, 11, 15, 19, 20}, y = {4.26794919, 39.35424869, 106.68337521, 206.12701665, 337.64110106, 375.52786405};
    std::vector<double> res(6);
    res = newton_interpolation(x, y);
    get_and_write_Yvalues("Newton_x2_x_sqx.txt", x, y, res, 100);
}
TEST(lagrange, x) { //x
    std::vector<double> x = {0, 1, 2, 3, 4, 5}, y = {0, 1, 2, 3, 4, 5};
    Lagrange_interpolation graph1(x, 6);
    graph1.calc_and_write_mesh("Lagrange_x.txt", y, 100);
}
TEST(lagrange, x2) { // x**2
    std::vector<double> x = {0, 1, 2, 3, 4, 5}, y = {0, 1, 4, 9, 16, 25};
    Lagrange_interpolation graph1(x, 6);
    graph1.calc_and_write_mesh("Lagrange_x2.txt", y, 100);
}
TEST(lagrange, x2_x_sqx) { // x**2 - x - x^(1/2)
    std::vector<double> x = {3, 7, 11, 15, 19, 20}, y = {4.26794919, 39.35424869, 106.68337521, 206.12701665, 337.64110106, 375.52786405};
    Lagrange_interpolation graph1(x, 6);
    graph1.calc_and_write_mesh("Lagrange_x2_x_sqx.txt", y, 100);
}
TEST(cubic_spline, x) { //x
    std::vector<double> x = {0, 1, 2, 3, 4, 5, 6}, y = {0, 1, 2, 3, 4, 5, 6};
    std::vector<double> res;
    res = cubic_spline_interopolation(x, y);
}
TEST(cubic_spline, x2) { //x**2
    std::vector<double> x = {0, 1, 2, 3, 4, 5, 6}, y = {0, 1, 4, 9, 16, 25, 36};
    std::vector<double> res;
    res = cubic_spline_interopolation(x, y);
}