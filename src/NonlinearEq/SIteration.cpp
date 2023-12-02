#include "SIteration.hpp"

const double simple_iteration(const double x_0, const double tolerance, const double (*y)(const double)) {
    double x = x_0, temp = y(x_0);
    while (fabs(x - temp) >= tolerance) {
        x = temp;
        temp = y(temp);
    }
    return x;
}