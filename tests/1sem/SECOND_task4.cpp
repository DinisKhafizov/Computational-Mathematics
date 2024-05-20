#include <gtest/gtest.h>

#include "../src/Interpolation/Cubic_spline.hpp" // for Sweep only..
#include "../src/ODE/Runge-kutta.hpp"

std::vector<double> g4(const double t, const std::vector<double> &x) {
    std::vector<double> result(2);
    result[0] = x[1];
    result[1] = 1 / (t + M_E) / (t + M_E) + M_E / log(t + M_E) * x[0] * x[0] -
                exp(x[1]) * x[0];
    return result;
}

const double p(const double t) {
    return -(t + M_E) * (t + M_E) * M_E * log(t + M_E);
}

const double q(const double t) { return M_E * (t + M_E); }

const double r(const double t) {
    return 1 / (t + M_E) / (t + M_E) - 1 / (t + M_E);
}

std::vector<double> g4_2(const double t, const std::vector<double> &x) {
    std::vector<double> result(2);
    result[0] = x[1];
    result[1] = p(t) * x[1] + q(t) * x[0] + r(t);
    return result;
}

TEST(homework_2, num_4) {
    const size_t N = 1e3 + 1;
    const double h = 8e-3;
    std::vector<double> t = linspace(0, 8, N);
    double temp_y = 5, epsilon = 1e-3, left = -3, right = 3;

    // Shooting method (everywhere using classical R-K 4)

    std::vector<double> u_0 = {M_E, 0}, result;
    while (fabs(temp_y) > epsilon) {
        u_0[1] = (left + right) / 2;
        result = runge_kutta_4_sys(t, g4, u_0, h);
        temp_y = result[584 * 2] - 2 * M_E * M_E;
        if (temp_y < 0) {
            left = (left + right) / 2;
        } else if (temp_y > 0) {
            right = (right + left) / 2;
        }
    }
    WritingToFile_sys(result, "ODE_task4_1.txt", 1);
    std::cout << "\n\033[32mShooting method\n\033[0m"
              << "y(e + 0.5) = " << result[125 * 2]
              << "\ny(e + 1) = " << result[188 * 2]
              << "\ny(e + 1.5) = " << result[250 * 2]
              << "\ny(e + 2) = " << result[313 * 2]
              << "\ny(e + 2.5) = " << result[375 * 2] << std::endl << std::endl;

    // Newton method (zero approx) using shooting method 
    
    u_0[0] = 0;
    left = -3;
    right = 3;
    temp_y = 5;
    while (fabs(temp_y) > epsilon) {
        u_0[1] = (left + right) / 2;
        result = runge_kutta_4_sys(t, g4_2, u_0, h);
        temp_y = result[584 * 2];
        if (temp_y < 0) {
            left = (left + right) / 2;
        } else if (temp_y > 0) {
            right = (right + left) / 2;
        }
    }
    WritingToFile_sys(result, "ODE_task4_2.txt", 1);
    std::cout << "\033[32mNewton method with shooting\n\033[0m"
              << "y(e + 0.5) = "
              << result[125 * 2] + (t[125] + M_E) * log(t[125] + M_E)
              << "\ny(e + 1) = "
              << result[188 * 2] + (t[188] + M_E) * log(t[188] + M_E)
              << "\ny(e + 1.5) = "
              << result[250 * 2] + (t[250] + M_E) * log(t[250] + M_E)
              << "\ny(e + 2) = "
              << result[313 * 2] + (t[313] + M_E) * log(t[313] + M_E)
              << "\ny(e + 2.5) = "
              << result[375 * 2] + (t[375] + M_E) * log(t[375] + M_E)
              << std::endl << std::endl;

    //Newton method using sweep

    std::vector<double> c(N - 1, 1), a(N - 1, 0), b(N, 1), f(N);
    c.push_back(0);
    for (size_t i = 0, end = N - 2; i < end; ++i) {
        a[i + 1] = 1 - p(t[i]) * h;
        b[i + 1] = -2 + p(t[i]) * h - h * h * q(t[i]);
        f[i] = h * h * r(t[i]);
    }
    TridiagonalMatrix A(c, b, a);
    result = Sweep(A, f);
    WritingToFile_sys(result, "ODE_task4_2.txt", 2);
    std::cout << "\033[32mNewton method with sweep\n\033[0m"
              << "y(e + 0.5) = "
              << result[125] + (t[125] + M_E) * log(t[125] + M_E)
              << "\ny(e + 1) = "
              << result[188] + (t[188] + M_E) * log(t[188] + M_E)
              << "\ny(e + 1.5) = "
              << result[250] + (t[250] + M_E) * log(t[250] + M_E)
              << "\ny(e + 2) = "
              << result[313] + (t[313] + M_E) * log(t[313] + M_E)
              << "\ny(e + 2.5) = "
              << result[375] + (t[375] + M_E) * log(t[375] + M_E)
              << std::endl << std::endl;
}