#include "Newton_eq.hpp"

const double newton_eq(const double (*y)(const double), const double (*y_der)(const double), const double x_0, const double tolerance) {
    double x = x_0, temp = x_0 - y(x_0)/y_der(x_0);
    while(fabs(x - temp) >= tolerance) {
        x = temp;
        temp = temp - y(temp)/y_der(temp);
    }
    return x;
}