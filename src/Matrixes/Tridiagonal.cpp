#include "Tridiagonal.hpp"

int TridiagonalMatrix::GetN() { return n; }

double TridiagonalMatrix::GetA(int i) { return A[i]; }
double TridiagonalMatrix::GetB(int i) { return B[i]; }
double TridiagonalMatrix::GetC(int i) { return C[i]; }

std::vector<double> Sweep(TridiagonalMatrix &Matrix,
                          const std::vector<double> &D) {
    double denom;  // переменная для хранения знаменателя
    const int N = Matrix.GetN();  // размер входной матрицы
    std::vector<double> x(N), p(N - 1),
        q(N - 1);  // х - результат прогонки, p q - векторы для вычисления х
    p[0] = -Matrix.GetC(0) / Matrix.GetB(0);
    q[0] = D[0] / Matrix.GetB(0);
    for (int i = 0; i < N - 2; ++i) {
        denom = Matrix.GetA(i) * p[i] + Matrix.GetB(i + 1);
        p[i + 1] = -Matrix.GetC(i + 1) / denom;
        q[i + 1] = (D[i + 1] - Matrix.GetA(i) * q[i]) / denom;
    }
    x[N - 1] = (D[N - 1] - Matrix.GetA(N - 2) * q[N - 2]) /
               (Matrix.GetA(N - 2) * p[N - 2] + Matrix.GetB(N - 1));
    for (int i = 1; i < N; ++i) {
        x[N - i - 1] = x[N - i] * p[N - i - 1] + q[N - i - 1];
    }
    return x;
}