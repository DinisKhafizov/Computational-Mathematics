#include <gtest/gtest.h>
#include "../src/ODE/Runge-kutta.hpp"

// Shooting method using dichotomy method

std::vector<double> f_1(const double t, const std::vector<double> &x) {
    std::vector<double> result(2);
    result[0] = x[1];
    result[1] = t * sqrt(x[0]);
    return result;
}

TEST(homework_2, num_1) {
    std::vector<double> t = linspace(0, 1, 1e4 + 1), u_0 = {0, 0}, result;
    double left = 0, right = 100, temp_y = 5, epsilon = 1e-3;
    while (fabs(temp_y) > epsilon) {
        u_0[1] = (left + right)/2;
        result = runge_kutta_4_sys(t, f_1, u_0, 1e-4);
        temp_y = result[1e4 - 2] - 2;
        if (temp_y > 0) {
            right = (left + right)/2;
        }
        else {
            left = (left + right)/2;
        }
    }
    WritingToFile_sys(result, "ODE_task1.txt", 1);
}