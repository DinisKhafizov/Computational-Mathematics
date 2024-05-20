#include "funcs.hpp"


// Граничные условия периодичны
std::vector<double> ugolok_left(const double (*y)(const double), double h, double tau, double L, double c, double time) {
    int N = static_cast<int>(L/h);
    int M = static_cast<int>(time/tau);
    std::vector<double> u_1(N), u_2(N);
    for (int i = 0; i < N; ++i) {
        u_1[i] = y(i * h);
    }
    for (size_t i = 1; i < M; ++i) {
        u_2[0] = u_1[0];
        for (size_t j = 1; j < N; ++j) {
            u_2[j] =  u_1[j] - c * tau / h * (u_1[j] - u_1[j - 1]);
        }
        u_1 = u_2;
    }
    return u_2;
}

std::vector<double> ugolok_left_task(const double (*y)(const double), double h, double tau, double L, double c, double time, std::vector<double> &TIMES) {
    int N = static_cast<int>(L/h);
    int M = static_cast<int>(time/tau);
    int counter = 0;
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
        for (size_t j = 1; j < N; ++j) {
            u_2[j] =  u_1[j] - c * tau / h * (u_1[j] - u_1[j - 1]);
        }
        u_2[0] = u_2[N - 1];
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