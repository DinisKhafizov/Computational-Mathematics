#include <gtest/gtest.h>
#include "Interpolation/Cubic_spline.hpp"
#include "Interpolation/Newton.hpp"

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
    std::cout << "USA 2023 population: " << y_true << std::endl;
    std::cout << "Cubic spline method prediction: " << y_pred_cubic << " with error = " << fabs(y_true - y_pred_cubic)/y_true * 100 << "%" << std::endl;
    std::cout << "Newton method prediciton: " << y_pred_newton << " with error = " << fabs(y_true - y_pred_newton)/y_true * 100 << "%" << std::endl << std::endl;
}