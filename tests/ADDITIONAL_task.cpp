#include <gtest/gtest.h>
#include "../src/ODE/Runge-kutta.hpp"

std::vector<double> f_add(const double t, const std::vector<double> &x) {
    std::vector<double> res(1);
    res[0] = sin(t);
    return res;
}

TEST(homework_3, ADDITIONAL) {
    std::vector<double> t = linspace(0, 100, 101), u_0 = {1};
    std::vector<double> result = runge_kutta_6(t, f_add, u_0, 1);
    std::vector<double> result2 = runge_kutta_2(t, f_add, u_0, 1);
    WritingToFile_sys(result, "ODE_add.txt", 2);
    WritingToFile_sys(result2, "ODE_add2.txt", 2);
}