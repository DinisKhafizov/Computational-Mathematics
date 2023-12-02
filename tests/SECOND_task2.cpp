#include <gtest/gtest.h>
#include "../src/ODE/Runge-kutta.hpp"

//Classical Runge-Kutta 4 orders with big splitting.

std::vector<double> f1(const double t, const std::vector<double> &x) {
    std::vector<double> result(2);
    const double mu = 1000;
    result[0] = x[1];
    result[1] = mu * (1 - x[1] * x[1]) * x[1] - x[0];
    return result;
}


TEST(homework_2, num2) {
    std::vector<double> t = linspace(0, 1e3, 1e6 + 1);
    std::vector<double> u_0 = {0, 1e-3};
    std::vector<double> result = runge_kutta_4_sys(t, f1, u_0, 1e-3);
    WritingToFile_sys(result, "ODE_task2.txt", 1);
}