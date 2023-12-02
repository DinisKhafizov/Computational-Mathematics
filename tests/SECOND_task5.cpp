#include <gtest/gtest.h>

#include "../src/ODE/Dorman-Prince.hpp"
#include <fstream>

std::vector<double> f5(const double t, const std::vector<double> &x) {
    double mu = 0.012277471;
    double eta = 1 - mu;
    double A = (x[0] + mu) * (x[0] + mu) + x[2] * x[2], B = (x[0] - eta) * (x[0] - eta) + x[2] * x[2];
    A = sqrt(A * A * A);
    B = sqrt(B * B * B);
    std::vector<double> result(4);
    result[0] = x[1];
    result[1] = x[0] + 2 * x[3] - eta * (x[0] + mu) / A - mu * (x[0] - eta)/ B;
    result[2] = x[3];
    result[3] = x[2] - 2 * x[1] - eta * x[2] / A - mu * x[2] / B;
    return result;
}

void WritingToFile_sys_DP(const std::pair<std::vector<double>, std::vector<double>> &res, std::string filename) {
    std::ofstream fout;
    fout.open(filename);
    for (size_t j = 0; j < 4; ++j) {
        for (size_t i = 0, end = static_cast<size_t>(size(res.first)/4); i < end; ++i) {
            fout << res.first[i * 4 + j] << ";";
        }
        fout << "\n";
    }
    for (size_t i = 0, end = size(res.second); i < end; ++i) {
        fout << res.second[i] << ";";
    }
}

TEST(homework_2, num_5) {
    const double T = 17.0652165601579625588917206249;
    std::vector<double> u_0 = {0.994, 0, 0, -2.00158510637908252240537862224};
    std::pair<std::vector<double>, std::vector<double>> result = dorman_prince(0, 5 * T, u_0, f5, 1e-3, 1e-4);
    WritingToFile_sys_DP(result, "ODE_task5.txt");
}