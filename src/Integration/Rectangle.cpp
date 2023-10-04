#include "Rectangle.hpp"

const double right_rectangle_method(const std::vector<double> &x, const double (*y)(const double)) {
    double integral = 0;
    for (size_t i = 0, end = size(x) - 1; i < end; ++i) {
        integral += y(x[i + 1]) * (x[i + 1] - x[i]);
    }
    return integral;
}

const double central_rectangle_method(const std::vector<double> &x, const double (*y)(const double)) {
    double integral = 0;
    for (size_t i = 0, end = size(x) - 1; i < end; ++i) {
        integral += y((x[i + 1] + x[i])/2.) * (x[i + 1] - x[i]);
    }  
    return integral;
}
const double left_rectangle_method(const std::vector<double> &x, const double (*y)(const double)) {
    double integral = 0;
    for (size_t i = 0, end = size(x) - 1; i < end; ++i) {
        integral += y(x[i]) * (x[i + 1] - x[i]);
    }
    return integral;
}