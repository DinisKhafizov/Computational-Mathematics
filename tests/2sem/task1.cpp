#include <gtest/gtest.h>
#include <fstream>
#include <string>

#include "../src/2sem/funcs.hpp"

void WriteToFile(std::vector<double> &result, std::string filename, int M, int N);
const double y(const double x);

int main() {

    // Константы для записи файлов
    std::string filename_ugolok = "../../tasks/hw1/ugolok_CFL_", filename_laks = "../../tasks/hw1/laks_CFL_", end_file = ".txt";

    // Константы вычисления
    const double L = 20, T = 18, h = 0.5, c = 1;
    std::vector<double> CFL = {0.6, 1, 1.01}, TIMES = linspace(0, 18, 30);
    std::vector<double> tau(CFL.size()), result;
    int N = static_cast<int>(L/h);

    // Вычисление tau по CFL
    for (size_t i = 0, end = CFL.size(); i < end; ++i) {
        tau[i] = CFL[i] * h / c; 
    }

    // Вычисление уголком для всех CFL
    for (size_t i = 0, end = CFL.size(); i < end; ++i) {
        result = ugolok_left_task(y, h, tau[i], L, c, T, TIMES);
        WriteToFile(result, filename_ugolok + std::to_string(i) + end_file, TIMES.size(), N);
    }

    // Вычисление по схеме Лакса-Вендроффа для всех CFL
    for (size_t i = 0, end = CFL.size(); i < end; ++i) {
        result = laks_vendroff_task(y, h, tau[i], L, c, T, TIMES);
        WriteToFile(result, filename_laks + std::to_string(i) + end_file, TIMES.size(), N);
    }
    return 0;
}

const double y(const double x) {
    const double L = 20;
    const double res = sin(4 * M_PI * x / L);
    return res;
}

void WriteToFile(std::vector<double> &result, std::string filename, int M, int N) {
    std::ofstream file(filename);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            file << result[i * N + j] << ";";
        }
        file << "\n";
    }
    file.close();
}