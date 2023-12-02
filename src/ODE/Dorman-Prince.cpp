#include "Dorman-Prince.hpp"

std::pair<std::vector<double>, std::vector<double>> dorman_prince(
    const double a, const double b, const std::vector<double> &u_0,
    std::vector<double> (*f)(const double, const std::vector<double> &),
    const double epsilon, const double h0) {
    const std::vector<double> a_ij = {1. / 5, 0, 0, 0, 0, 0,
                                      3. / 40, 9. / 40, 0, 0, 0, 0,
                                      44. / 45, -56. / 15, 32. / 9, 0, 0, 0, 
                                      19372. / 6561, -25360. / 2187, 64448. / 6561, -212. / 729, 0, 0,
                                      9017. / 3168, -355. / 33, 46732. / 5247, 49. / 176, -5103. / 18656, 0,
                                      35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84};
    const std::vector<double> c = {1. / 5, 3. / 10, 4. / 5, 8. / 9, 1, 1};
    const std::vector<double> b1 = {35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0};
    const std::vector<double> b2 = {5179. / 57600, 0, 7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40};
    const int N = size(u_0);
    double t = a, h = h0, error, norm, h_opt;
    Matrix K(N, 7), A(a_ij, 6, 6);
    std::vector<double> result, temp, x1, x2, x, T;
    size_t counter = 0;
    for (size_t i = 0; i < N; ++i) {
        result.push_back(u_0[i]);
    }
    T.push_back(a);
    while (t < b) {
        temp = get_vectors_part(result, counter * N, (counter + 1) * N);
        K.change_Col(0, 0, f(t, temp));
        x1 = temp + h * b1[0] * K.getCol(0);
        x2 = temp + h * b2[0] * K.getCol(0);
        for (size_t i = 1; i < 7; ++i) {
            x = temp;
            for (size_t j = 0; j < i; ++j) {
                x = x + h * A(i - 1, j) * K.getCol(j);
            }
            K.change_Col(i, 0, f(t + h * c[i - 1], x));
            x1 = x1 + h * b1[i] * K.getCol(i);
            x2 = x2 + h * b2[i] * K.getCol(i);
        }
        error = endless_norm(x1 - x2);
        h_opt = 0.9 * h * pow(epsilon / error, 0.2);
        if (h_opt < h / 2) {
            h = h_opt;
        } else {
            t += h;
            T.push_back(t);
            result.resize(size(result) + N);
            for (size_t i = 0; i < N; ++i) {
                result[(counter + 1) * N + i] = x1[i];
            }
            h = std::min(h_opt, h0);
            if (t + h > b) {
                h = b - t;
            }
            ++counter;
        }
    }
    std::pair<std::vector<double>, std::vector<double>> res;
    res.first = result;
    res.second = T;
    return res;
}