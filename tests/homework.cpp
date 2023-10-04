#include <gtest/gtest.h>
#include <cmath>
#include "Interpolation/Cubic_spline.hpp"
#include "Interpolation/Newton.hpp"
#include "Integration/Quadrature.hpp"
#include "Integration/Rectangle.hpp"

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
    double x = x_0;
    while (fabs((x - y_der(x))) >= tolerance) {
        x = y_der(x);
    }
    return x;
}

const double simple_iteration2(const double (*op)(const double, const double),const double x_max, const double y_max, const double tolerance, const double shift) {
    double x = x_max + shift;
    while (fabs((x - op(y_max, x))) >= tolerance) {
        x = op(y_max, x);
    }
    return x;
}

TEST(homework, num2) {
    double x = 0, tolerance = 0.001;
    x = simple_iteration1(x, tolerance);
    double f_max = y_1(x);
    double x1 = simple_iteration2(y_2, x, f_max, tolerance/2, -0.2);
    double x2 = simple_iteration2(y_3, x, f_max, tolerance/2, 0.2);
    std::cout << "x1 = " << x1 << std::endl;
    std::cout << "x2 = " << x2 << std::endl;
    std::cout << "delta = x2 - x1 = " << x2 - x1 << std::endl;
}

TEST(homework, num3) {
    std::vector<double> x = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};
    std::vector<double> y = {92228496, 106021537, 123202624, 132164569, 151325798, 179323175, 203211926, 226545805, 248709873, 281421906};
    double x_predicting = 2023, y_true = 336e6;
    double delta = x_predicting - x.back();
    std::vector<double> res_cubic = cubic_spline_interopolation(x, y);
    std::vector<double> res_newton = newton_interpolation(x, y);
    double y_pred_newton = get_polynoms_value(x_predicting, res_newton, x);
    double y_pred_cubic = res_cubic[32] + res_cubic[33] * delta + res_cubic[34] * delta  * delta / 2 + res_cubic[35] * delta * delta * delta / 6;
    std::cout << std::endl;
    std::cout << "USA 2023 population: " << y_true << std::endl << std::endl;
    std::cout << "Cubic spline method prediction: " << y_pred_cubic << " with error = " << fabs(y_true - y_pred_cubic)/y_true * 100 << "%" << std::endl << std::endl;
    std::cout << "Newton method prediciton: " << y_pred_newton << " with error = " << fabs(y_true - y_pred_newton)/y_true * 100 << "%" << std::endl;
    std::cout << std::endl;  
}

const double u1(const double x) {
    return sqrt(1 - x * x);
}
const double u2(const double y) {
    return atan(y);
}
const std::pair<double, double> simple_iteration_sys(const double x0, const double y0, const double tolerance){
    double x = x0, y = y0, temp;
    while(fabs(x - u2(y)) >= tolerance && fabs(y - u1(x)) >= tolerance) {
        temp = x;
        x = u2(y);
        y = u1(temp); 
    }
    std::pair<double, double> res = {x, y};
    return res;
}

TEST(homework, num4) {
    double x0 = 0.5, y0 = 0.5;
    std::pair<double, double> res = simple_iteration_sys(x0, y0, 0.000001);
    std::cout << res.first << " " << res.second;
}

const double y_5(const double x) {
    return sin(100 * x) * exp(-x * x) * cos(2 * x); 
}

TEST(homework, num5) {
    const size_t number_of_dots = 2000;
    const int a = 0, b = 3;
    std::vector<double> x_lin = linspace(a, b, number_of_dots);
    double I_left_rect = left_rectangle_method(x_lin, y_5), 
           I_right_rect = right_rectangle_method(x_lin, y_5), 
           I_cent_rect = central_rectangle_method(x_lin, y_5);
    double I_trap = 0, 
           I_simp = 0, 
           I_th_eight = 0, 
           I_gauss = 0;
    for (size_t i = 0, end = number_of_dots - 1; i < end; ++i) {
        I_trap += trapezoid(x_lin[i], x_lin[i + 1], y_5);
        I_simp += simpson(x_lin[i], x_lin[i + 1], y_5);
        I_th_eight += three_eights(x_lin[i], x_lin[i + 1], y_5);
        I_gauss += gauss_quadrature(x_lin[i], x_lin[i + 1], 3, y_5);
    }
    std::cout << std::endl;
}