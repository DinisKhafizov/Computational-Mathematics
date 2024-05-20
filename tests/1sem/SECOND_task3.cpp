#include <gtest/gtest.h>


#include "../src/ODE/Runge-kutta.hpp"

// Shooting method, using dichotomy and claasical R-K 4

std::vector<double> g(const double t, const std::vector<double> &x) {
    std::vector<double> result(2);
    const double cos_ = cos(t), sin_ = sin(t), t2 = t * t, exp_ = exp(t);
    result[0] = x[1];
    result[1] = -(t2 - 3) * (x[1] + x[0] * cos_) + 2 - 6 * t + 2 * t2 * t +
                (t2 - 3) * exp_ * sin_ * (1 + cos_) +
                cos_ * (exp_ + (t2 - 1) + t2 * t2 - 3 * t2);
    return result;
}


TEST(homework_2, num_3) {
    std::vector<double> t = linspace(0, 4, 1e3 + 1), u_0 = {0, -1e-3}, result;
    double temp_y = 5, epsilon = 1e-3, left = -100, right = 100, pi2 = M_PI * M_PI;
    while (fabs(temp_y) > epsilon) {
        u_0[1] = (left + right) / 2;
        result = runge_kutta_4_sys(t, g, u_0, 4e-3);
        temp_y = result[785 * 2] - pi2;
        if (temp_y < 0) {
            left = (left + right) / 2;
        } else if (temp_y > 0) {
            right = (right + left) / 2;
        }
    }
    std::cout << "y(0.5) = " << result[125 * 2]
              << "\ny(1) = " << result[250 * 2]
              << "\ny(1.5) = " << result[375 * 2]
              << "\ny(2) = " << result[500 * 2]
              << "\ny(2.5) = " << result[625 * 2]
              << "\ny(3) = " << result[750 * 2] << std::endl;
    WritingToFile_sys(result, "ODE_task3.txt", 1);
}

int s21_abs(int x) {
    return x > 0? x : -x;
}

int s21_floor(double x) {

}