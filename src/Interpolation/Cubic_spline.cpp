#include "Cubic_spline.hpp"

#include <fstream>
#include <iostream>
#include <string>

// интерполирующая функция, принимает на вход точки, значения функции в этих
// точках, значения вторых производных на концах отрезка (гран условия, по
// умолчанию == 0) возвращает вектор размера 4n - 4 - коэффициенты a, b, c, d,
// идущие подряд.
// s_i(x) = a_i + b_i(x - x_i) + c_i/2 (x-x_i)^2 + d_i/6 (x - x_i)^3
std::vector<double> Interpolation(std::vector<double> &x,
                                  std::vector<double> &y, double sec_der_begin,
                                  double sec_der_end) {
    // a, b, c, f - векторы для Трехдиагональной матрицы и прогонки (состоят из
    // h_i), A_koef, B_koef, C_koef, D_koef - векторы для хранения коэф.
    // интерполяции
    std::vector<double> a(size(x) - 1), b(size(x), 2), c(size(x) - 1), f(size(x)),
        C_koef(size(x)), B_koef(size(x) - 1), A_koef(size(x) - 1),
        D_koef(size(x) - 1);
    // вектор для записи результата
    std::vector<double> res(4 * size(x) - 4);
    // заносим в векторы трехдиагональной матрицы краевые значения.
    a.back() = 0;
    c[0] = 0;
    f[0] = sec_der_begin;
    f.back() = sec_der_end;
    // заносим в векторы трехдиагональной матрицы значения h_i
    for (int i = 0; i < static_cast<int>(size(c)) - 1; ++i) {
        c[i + 1] = x[i + 2] - x[i + 1];
        a[i] = x[i + 1] - x[i];
        f[i + 1] = 6 * ((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) -
                        (y[i + 1] - y[i]) / (x[i + 1] - x[i]));
    }
    // инициализируем трехдиагональную матрицу
    TridiagonalMatrix A(a, b, c);
    // прогонка
    C_koef = Sweep(A, f);
    // вычисляем коэффициенты интерполяции
    for (int i = 0; i < static_cast<int>(size(B_koef)); ++i) {
        A_koef[i] = y[i];
        D_koef[i] = (C_koef[i + 1] - C_koef[i]) / (x[i + 1] - x[i]);
        B_koef[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) +
                    (x[i + 1] - x[i]) * C_koef[i + 1] / 3 +
                    (x[i + 1] - x[i]) * C_koef[i] / 6;
    }
    // заносим результат и возрващаем его
    for (int i = 0; i < int(size(res) / 4); ++i) {
        res[i * 4] = A_koef[i];
        res[i * 4 + 1] = B_koef[i];
        res[i * 4 + 2] = C_koef[i + 1];
        res[i * 4 + 3] = D_koef[i];
    }
    return res;
}

void WritingToFile(std::vector<double> res) {
    std::ofstream fout;
    fout.open("Interpolation.txt");
    for (int i = 0; i < int(size(res) / 4) - 1; ++i) {
        fout << res[i * 4] << ";" << res[i * 4 + 1] << ";" << res[i * 4 + 2] 
             << ";" << res[i * 4 + 3]  << std::endl;
    }
    fout << res[(int(size(res) / 4) - 1) * 4] << ";"
         << res[(int(size(res) / 4) - 1) * 4 + 1] << ";"
         << res[(int(size(res) / 4) - 1) * 4 + 2] << ";"
         << res[(int(size(res) / 4) - 1) * 4 + 3];
}

void Interpolate_and_write(std::vector<double> &x, std::vector<double> &y,
                           double sec_der_begin, double sec_der_end) {
    std::vector<double> res = Interpolation(x, y, sec_der_begin, sec_der_end);
    WritingToFile(res);
}

// boundary conditions - f''(x_0) = f''(x_n) = 0
std::vector<double> cubic_spline_interopolation(const std::vector<double> &x, const std::vector<double> &y) {
    const size_t N = size(x) - 1;
    std::vector<double> up_diag(N - 2), diag(N - 1, 2), und_diag(N - 2), f(N - 1);
    std::vector<double> a(N), b(N), c(N - 1), d(N);
    for (size_t i = 0; i < N; ++i) {
        a[i] = y[i + 1];
    }
    for (size_t i = 0, end = N - 2; i < end; ++i) {
        up_diag[i] = (x[i + 2] - x[i + 1])/(x[i + 2] - x[i]);
        und_diag[i] = (x[i + 2] - x[i + 1])/(x[i + 3] - x[i + 1]);
        f[i] = ((y[i + 2] - y[i + 1])/(x[i + 2] - x[i + 1]) - (y[i + 1] - y[i])/(x[i + 1] - x[i]))/(x[i + 2] - x[i]);
        f[i] *= 6;
    }
    f[N - 2] = 6 * ((y[N] - y[N - 1])/(x[N] - x[N - 1]) - (y[N - 1] - y[N - 2])/(x[N - 1] - x[N - 2]))/(x[N] - x[N - 2]);
    TridiagonalMatrix K(und_diag, diag, up_diag);
    c = Sweep(K, f);
    c.resize(N);
    b[0] = c[0] * (x[1] - x[0]) / 3  + (y[1] - y[0])/(x[1] - x[0]);
    d[0] = c[0] / (x[1] - x[0]);
    for (size_t i = 1; i < N; ++i) {
        b[i] = (x[i + 1] - x[i]) * (c[i] / 3 + c[i - 1] / 6) + (y[i + 1] - y[i])/(x[i + 1] - x[i]);
        d[i] = (c[i] - c[i - 1])/(x[i + 1] - x[i]);
    }
    std::vector<double> res(4 * N);
    for (size_t i = 0; i < N; ++i) {
        res[i * 4] = a[i];
        res[i * 4 + 1] = b[i];
        res[i * 4 + 2] = c[i];
        res[i * 4 + 3] = d[i];
    }
    return res;
}
