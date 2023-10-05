#include <gtest/gtest.h>
#include "Integration/Quadrature.hpp"
#include "Integration/Rectangle.hpp"
#include "UsefulOperations/VectorOperations.hpp"

const double y_5(const double x) {
    return sin(100 * x) * exp(-x * x) * cos(2 * x); 
}

TEST(homework, num5) {
    const size_t number_of_dots = 5000;
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
    std::cout << std::endl <<"For " << number_of_dots << " number of dots:" << std::endl;
    std::cout <<  "Left rectangle: I = " << I_left_rect << std::endl;
    std::cout << "Central rectangle: I = " << I_cent_rect << std::endl;
    std::cout << "Right rectangle: I = " << I_right_rect << std::endl;
    std::cout << "Trapezoid: I = " << I_trap << std::endl;
    std::cout << "Simpson: I = " << I_simp << std::endl;
    std::cout << "Three eights: I = " << I_th_eight << std::endl;
    std::cout << "Gauss (third degree): I = " << I_gauss << std::endl << std::endl;
}