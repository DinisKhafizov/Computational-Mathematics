#include "CG.hpp"

std::vector<double> CGD(const Matrix &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance) {
    std::vector<double> x = x_0, r = A * x_0 - b;
    std::vector<double> d = r;
    double alpha, koef1 = r * r, koef2;
    while(first_norm(r) > tolerance) {
        alpha = (d * r)/(d * (A * d));
        x = x - alpha * d;
        r = A * x - b;
        koef2 = r * r;
        d = r + d * koef2/koef1;
        koef1 = koef2; 
    }
    return x;
}
std::vector<double> CGD_precond(const Matrix &A, const std::vector<double> &x_0, const std::vector<double> &b, double tolerance) {
    std::vector<double> x = x_0, r = A * x_0 - b;
    std::vector<double> w = CGD(A, x_0, r, tolerance);
    std::vector<double> d = w, r_last, w_last;
    double alpha, beta, norm = first_norm(r);

    while(norm > tolerance) {
        r_last = r;
        w_last = w;
        alpha = (r * w)/(d * (A * d));
        x = x - alpha * d;
        r = A * x - b;
        norm = first_norm(r);
        if (norm > tolerance) {
            w = CGD(A, x_0, r, tolerance);
            beta = (r * w)/(r_last * w_last);
            d = w + beta * d;
        }
    }
    return x;
}