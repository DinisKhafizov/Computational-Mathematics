#include "funcs.hpp"

std::vector<double> laks_vendroff(const double (*y)(const double), double h, double tau, double L, double c, double time) {
    int N = static_cast<int>(L/h);
    int M = static_cast<int>(time/tau);
    std::vector<double> u_1(N), u_2(N);
    for (int i = 0; i < N; ++i) {
        u_1[i] = y(i * h);
    }

    for (size_t i = 1; i < M; ++i) {
        u_2[0] = u_1[0];
        for (size_t j = 1; j < N - 1; ++j) {
            u_2[j] =  c * c * tau * tau / h / h / 2 * (u_1[j + 1] - 2 * u_1[j] + u_1[j - 1]) - c * tau / h / 2 * (u_1[j + 1] - u_1[j - 1]) + u_1[j];
        }
        u_2[N - 1] = u_2[N - 2];
        u_1 = u_2;
    }
    return u_2;
}

std::vector<double> laks_vendroff_task(const double (*y)(const double), double h, double tau, double L, double c, double time, std::vector<double> &TIMES) {
    int N = static_cast<int>(L/h);
    int M = static_cast<int>(time/tau);
    int counter = 0;
    double c1 = (c * c * tau * tau) / (2 * h * h);
    double c2 = (c * tau) / (2 * h);
    std::vector<double> u_1(N), u_2(N), result(N * TIMES.size());
    for (int i = 0; i < N; ++i) {
        u_1[i] = y(i * h);
    }
    if (TIMES[0] == 0) {
        ++counter;
        for (int i = 0; i < N; ++i) {
            result[i] = u_1[i];
        }
    }
    for (size_t i = 1; i < M; ++i) {
        u_2[0] = c1 * (u_1[1] - 2 * u_1[0] + u_1[N - 1]) - c2 * (u_1[1] - u_1[N - 1]) + u_1[0];
        for (size_t j = 1, end = N - 1; j < end; ++j) {
            u_2[j] =  c1 * (u_1[j + 1] - 2 * u_1[j] + u_1[j - 1]) - c2 * (u_1[j + 1] - u_1[j - 1]) + u_1[j];
        }
        u_2[N - 1] = c1 * (u_1[0] - 2 * u_1[N - 1] + u_1[N - 2]) - c2 * (u_1[0] - u_1[N - 2]) + u_1[N - 1];
        u_1 = u_2;
        if ((TIMES[counter] >= (i - 1) * tau) && (TIMES[counter] <= i * tau)) {
            for (int i = 0; i < N; ++i) {
                result[i + counter * N] = u_1[i];
            }
            ++counter;
        }
    }
    return result;
}