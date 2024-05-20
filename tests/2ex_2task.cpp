#include <stdio.h>

#include <cmath>
#include <fstream>

#include "../src/Matrixes/StandardMatrix.hpp"
#include "../src/Matrixes/Tridiagonal.hpp"
#include "../src/UsefulOperations/VectorOperations.hpp"

// Граничные условия при x = 0, x = 1.
double phi_at_x(double y, double t, double lambda, double x) {
  return pow(-1, x) * sin(5 * M_PI * y) * exp(-50 * M_PI * M_PI * lambda * t);
}

// Начальные условия.
double phi_0(double x, double y) { return cos(M_PI * x) * sin(5 * M_PI * y); }

void writeToFile(std::string filename, Matrix &A, int flag) {
  std::ofstream file;
  if (flag) {
    file.open(filename);
    for (int i = 0; i < A.GetM(); ++i) {
      for (int j = 0, end = A.GetN() - 1; j < end; ++j) {
        file << A(i, j) << ";";
      }
      file << A(i, A.GetN() - 1) << "\n";
    }
  }
}

int main() {
  // Иницализация констант
  const double lambda = 1e-4;
  const double alpha = lambda, beta = 25 * lambda;
  const double X = 1, Y = 1;
  const double tau_frames = 0.1, tau = 1e-3, T = 1;
  std::vector<double> Nxs = {10, 20, 50, 100, 200, 500, 1000};
  std::vector<std::string> files =
  { "../../tasks/hw2/second_task_1.txt",
    "../../tasks/hw2/second_task_2.txt",
    "../../tasks/hw2/second_task_3.txt",
    "../../tasks/hw2/second_task_4.txt",
    "../../tasks/hw2/second_task_5.txt",
    "../../tasks/hw2/second_task_6.txt",
    "../../tasks/hw2/second_task_7.txt" };

  for (int _ = 0; _ < 7; ++_) {
    const double Nx = Nxs[_], Ny = Nxs[_], X = 1, Y = 1;
    const double hx = X / (Nx - 1), hy = Y / (Ny - 1);
    const double tau = 1e-3, T = 1;

    // Иницализация трехдиагональных матриц для расчета методом переменных
    // направлений (схема Писмена-Речфорда)
    std::vector<double> a_x(Nx - 1, -alpha / hx / hx),
        b_x(Nx, 2 / tau + 2 * alpha / hx / hx), c_x(Nx - 1, -alpha / hx / hx),
        d_x(Nx);
    a_x.back() = 0;
    b_x[0] = 1;
    b_x.back() = 1;
    c_x[0] = 0;
    std::vector<double> a_y(Ny - 1, -beta / hy / hy),
        b_y(Ny, 2 / tau + 2 * beta / hy / hy), c_y(Ny - 1, -beta / hy / hy),
        d_y(Ny);
    a_y.back() = 0;
    b_y[0] = 1;
    b_y.back() = 1;
    c_y[0] = 0;
    TridiagonalMatrix matrX(a_x, b_x, c_x), matrY(a_y, b_y, c_y);

    // Иницализация переменных
    double t = 0;
    int flag = 1;
    Matrix phi(Nx, Ny);
    for (int i = 0; i < Nx; ++i) {
      for (int j = 0; j < Ny; ++j) {
        phi(i, j) = phi_0(hx * i, hy * j);
      }
    }
    Matrix tempX(Nx, Ny);
    std::vector<double> temp_x0(Nx), temp_x1(Nx);
    std::vector<double> temp_y01(Ny, 0);

    for (; t < T; t += tau) {
      // Расчет матрицы на первом шаге
      for (int i = 0; i < Ny; ++i) {
        temp_x0[i] = phi_at_x(i * hy, t, lambda, 0);
        temp_x1[i] = phi_at_x(i * hy, t, lambda, 1);
      }
      tempX.setRow(0, temp_x0);
      tempX.setRow(Nx - 1, temp_x1);
      for (int i = 1, end = Ny - 1; i < end; ++i) {
        std::vector<double> temp_1 = phi.getRow(i - 1), temp_2 = phi.getRow(i),
                            temp_3 = phi.getRow(i + 1);
        std::vector<double> d =
            2 / tau * temp_2 + beta / hy / hy * (temp_3 - 2 * temp_2 + temp_1);
        d[0] = 0;
        d.back() = 0;
        std::vector<double> res_temp = Sweep(matrX, d);
        tempX.setRow(i, res_temp);
      }

      // Расчет матрицы на втором шаге
      phi.setCol(0, temp_y01);
      phi.setCol(Nx - 1, temp_y01);

      for (int i = 1, end = Ny - 1; i < end; ++i) {
        std::vector<double> temp_1 = tempX.getCol(i - 1),
                            temp_2 = tempX.getCol(i),
                            temp_3 = tempX.getCol(i + 1);
        std::vector<double> d =
            2 / tau * temp_2 + alpha / hx / hx * (temp_3 - 2 * temp_2 + temp_1);
        d[0] = phi_at_x(hy * i, t, lambda, 0);
        d.back() = phi_at_x(hy * i, t, lambda, 1);
        std::vector<double> res_temp = Sweep(matrY, d);
        phi.setCol(i, res_temp);
      }

      if (t >= 0.5) {
        writeToFile(files[_], phi, flag);
        flag = 0;
      }
    }
  }
}
