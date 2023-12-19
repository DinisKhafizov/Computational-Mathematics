#include "Runge-kutta.hpp"

std::vector<double> runge_kutta_4(const std::vector<double> &t,
                                  double (*f)(const double, const double),
                                  const double u0, const double h) {
    std::vector<double> result(size(t));
    double k1, k2, k3, k4;
    result[0] = u0;
    for (size_t i = 0, end = size(t) - 1; i < end; ++i) {
        k1 = f(t[i], result[i]);
        k2 = f(t[i] + h / 2, result[i] + h * k1 / 2);
        k3 = f(t[i] + h / 2, result[i] + h * k2 / 2);
        k4 = f(t[i] + h, result[i] + h * k3);
        result[i + 1] = result[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    return result;
}

std::vector<double> runge_kutta_4_sys(
    const std::vector<double> &t,
    std::vector<double> (*f)(const double, const std::vector<double> &),
    const std::vector<double> &u_0, const double h) {
    size_t N = size(u_0);
    std::vector<double> result(size(t) * N);
    std::vector<double> k1, k2, k3, k4;
    for (size_t i = 0; i < N; ++i) {
        result[i] = u_0[i];
    }
    for (size_t i = 0, end = size(t) - 1; i < end; ++i) {
        std::vector<double> temp(N);
        for (size_t j = 0; j < N; ++j) {
            temp[j] = result[i * N + j];
        }
        k1 = f(t[i], temp);
        k2 = f(t[i] + h / 2, temp + h / 2 * k1);
        k3 = f(t[i] + h / 2, temp + h / 2 * k2);
        k4 = f(t[i] + h, temp + h * k3);
        for (size_t j = 0; j < N; ++j) {
            result[(i + 1) * N + j] =
                result[i * N + j] +
                h / 6 * (k1[j] + 2 * (k2[j] + k3[j]) + k4[j]);
        }
    }
    return result;
}

std::vector<double> runge_kutta_6(
    const std::vector<double> &t,
    std::vector<double> (*f)(const double, const std::vector<double> &),
    const std::vector<double> &u_0, const double h) {
    size_t N = size(u_0);
    std::vector<double> result(size(t) * N);
    std::vector<double> k1, k2, k3, k4, k5, k6;
    for (size_t i = 0; i < N; ++i) {
        result[i] = u_0[i];
    }
    for (size_t i = 0, end = size(t) - 1; i < end; ++i) {
        std::vector<double> temp(N);
        for (size_t j = 0; j < N; ++j) {
            temp[j] = result[i * N + j];
        }
        k1 = f(t[i], temp);
        k2 = f(t[i] + h / 10, temp + h * k1 / 10);
        k3 = f(t[i] + 2 * h / 10, temp + 2. / 10 * h * (k1 + k2));
        k4 = f(t[i] + 3 * h / 10, temp + h / 40 * (3 * k1 + 9 * k2));
        k5 = f(t[i] + 8 * h / 10, temp + h / 10 * (3 * k1 - 9 * k2 + 6 * k3));
        k6 = f(t[i] + 9 * h / 10, temp + h * (-11. / 18 * k1 + 5. / 2 * k2 - 70. / 27 * k3 + 35. / 27 * k4 - 1631. / 55296 * k5));
        for (size_t j = 0; j < N; ++j) {
            result[(i + 1) * N + j] =
                result[i * N + j] +
                h * (37. / 378 * k1[j] + 250. / 621 * k3[j] + 125. / 594 * k4[j] + 512. / 1771 * k6[j]);
        }
    }
    return result;
}

std::vector<double> runge_kutta_2(
    const std::vector<double> &t,
    std::vector<double> (*f)(const double, const std::vector<double> &),
    const std::vector<double> &u_0, const double h) {
    size_t N = size(u_0);
    std::vector<double> result(size(t) * N);
    std::vector<double> k1, k2;
    for (size_t i = 0; i < N; ++i) {
        result[i] = u_0[i];
    }
    for (size_t i = 0, end = size(t) - 1; i < end; ++i) {
        std::vector<double> temp(N);
        for (size_t j = 0; j < N; ++j) {
            temp[j] = result[i * N + j];
        }
        k1 = f(t[i], temp);
        k2 = f(t[i] + h, temp + h * k1);
        for (size_t j = 0; j < N; ++j) {
            result[(i + 1) * N + j] =
                result[i * N + j] + h * k2[j];
        }
    }
    return result;
    }

void WritingToFile_sys(const std::vector<double> &res, std::string filename,
                       int type) {
    std::ofstream fout;
    fout.open(filename);
    if (type == 1) {
        for (size_t i = 0, end = static_cast<size_t>(size(res) / 2); i < end;
             ++i) {
            fout << res[i * 2] << ";";
        }
    } else if (type == 2) {
        for (size_t i = 0, end = size(res); i < end; ++i) {
            fout << res[i] << ";";
        }
    }
}