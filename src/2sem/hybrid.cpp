#include "funcs.hpp"

std::vector<double> hybrid(const double (*y)(const double), double h, double tau, double L, double c, double time) {
    int N = static_cast<int>(L/h);
    int M = static_cast<int>(time/tau);
    long double beta, alpha, laks, ugolok;
    std::vector<double> u_1(N), u_2(N);
    for (int i = 0; i < N; ++i) {
        u_1[i] = y(i * h);
    }
    for (int i = 1; i < M; ++i) {
        u_2[0] = 10;
        for (int j = 1, end = N - 1; j < end; ++j) {
            beta = fabs(u_1[j + 1] - 2 * u_1[j] + u_1[j - 1]) / (1e-9 + fabs((u_1[j + 1] - u_1[j]) * (u_1[j] - u_1[j - 1])));
            alpha = beta / (beta + 1);
            laks = c * c * tau * tau / h / h / 2 * (u_1[j + 1] - 2 * u_1[j] + u_1[j - 1]) - c * tau / h / 2 * (u_1[j + 1] - u_1[j - 1]) + u_1[j];
            ugolok = u_1[j] - c * tau / h * (u_1[j] - u_1[j - 1]);
            u_2[j] = alpha * ugolok + (1 - alpha) * laks;
        }
        u_2[N - 1] = u_2[N - 2];
        u_1 = u_2;
    }
    return u_2;
}