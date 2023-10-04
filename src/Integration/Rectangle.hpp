#ifndef RECTANGLE_INTEGRATION_METH
#define RECTANGLE_INTGRATION_METH
#include <vector>
#include <iostream>
const double right_rectangle_method(const std::vector<double> &x, const double (*y)(const double));
const double central_rectangle_method(const std::vector<double> &x, const double (*y)(const double));
const double left_rectangle_method(const std::vector<double> &x, const double (*y)(const double));

#endif