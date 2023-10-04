#include "SIteration.hpp"

const double y(const double x) {
    return 0;  // changable function!
}

double simple_iteration(const double x_0, const double tolerance) {
    double x = x_0;
    while ((x - y(x)) >= tolerance) {
        x = y(x);
    }
    return x;
}