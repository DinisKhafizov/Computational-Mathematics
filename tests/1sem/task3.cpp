#include <gtest/gtest.h>

#include "Interpolation/Cubic_spline.hpp"
#include "Interpolation/Newton.hpp"

TEST(homework, num3) {
    /*
    std::vector<double> x = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980,
    1990, 2000}; std::vector<double> y = {92228496, 106021537, 123202624,
    132164569, 151325798, 179323175, 203211926, 226545805, 248709873,
    281421906}; double x_predicting = 2023, y_true = 336e6; double delta =
    x_predicting - x.back(); std::vector<double> res_cubic =
    cubic_spline_interopolation(x, y); std::vector<double> res_newton =
    newton_interpolation(x, y); double y_pred_newton =
    get_polynoms_value(x_predicting, res_newton, x); double y_pred_cubic =
    res_cubic[32] + res_cubic[33] * delta + res_cubic[34] * delta  * delta / 2 +
    res_cubic[35] * delta * delta * delta / 6; std::cout << std::endl; std::cout
    << "USA 2023 population: " << y_true << std::endl; std::cout << "Cubic
    spline method prediction: " << y_pred_cubic << " with error = " <<
    fabs(y_true - y_pred_cubic)/y_true * 100 << "%" << std::endl; std::cout <<
    "Newton method prediciton: " << y_pred_newton << " with error = " <<
    fabs(y_true - y_pred_newton)/y_true * 100 << "%" << std::endl << std::endl;
    */
    std::vector<double> x = {0,   0.2,  0.4, 0.6,  0.8, 1,    1.2, 1.4,
                             1.6, 1.8,  2.1, 2.3,  2.5, 2.7,  2.8, 3,
                             3.2, 3.4,  3.6, 4,    4.1, 4.15, 4.2, 4.23,
                             4.3, 4.35, 4.4, 4.45, 4.6, 4.8,  5};
    std::vector<double> y = {-7.278, -7.189, -7.156, -7.056, -7.178,
                             -6.823, -5.834, -3.902, -2.08,  0.0750000000000002,
                             3.252,  3.974,  4.085,  2.674,  1.774,
                             0.253,  -1.491, -3.457, -5.279, -6.09,
                             -5.046, 0.208,  3.23,   5.218,  6.773,
                             5.584,  3.341,  1.197,  -1.18,  -7.223,
                             -7.278};
    std::vector<double> res = cubic_spline_interopolation(x, y);
    WritingToFile(res);
}