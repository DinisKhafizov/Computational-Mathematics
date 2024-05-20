#include <gtest/gtest.h>

#include <fstream>

#include "../src/Matrixes/Tridiagonal.hpp"
#include "../src/UsefulOperations/VectorOperations.hpp"

double getRho(double p, const double rho_0, const double c_f,
              const double p_0) {
  return rho_0 * (1 + c_f * (p - p_0));
}


std::ofstream file("../../tasks/hw2/task1.txt");

void writeToFile(std::vector<double> &x, double t) {
  file << t << ";";
  for (size_t i = 0, end = x.size() - 1; i < end; ++i) {
    file << x[i] << ";";
  }
  file << x.back() << "\n";
}

int main() {
  const size_t Nx = 100;
  const double dz = 10, L = 500;
  const double p_inj = 150e5, p_prod = 50e5;
  double tau = 3600, T = 3500e3, t, h = L / Nx;

  const double k = 10e-14, mu = 10e-3, phi = 0.2, c_f = 10e-9, rho_0 = 1000;
  const double p_01 = 120e5;
  int counter = 0;

  std::vector<double> rho(Nx), p(Nx, 100e5);
  std::vector<double> a(Nx, 1), b(Nx - 1), c(Nx - 1), d(Nx);

  int num_frames = 100;
  std::vector<double> record = linspace(0, T, num_frames);

  for (; t <= T; t += tau) {
    // rho
    for (size_t i = 0; i < Nx; ++i) {
      rho[i] = getRho(p[i], rho_0, c_f, p_01);
    }
    // c
    for (size_t i = 0, end = Nx - 2; i < end; ++i) {
      if (p[i] >= p[i + 1]) {
        c[i] = k * rho[i] / mu / h / h;
      } else {
        c[i] = k * rho[i + 1] / mu / h / h;
      }
    }
    c.back() = 0;
    // b
    b[0] = 0;
    for (size_t i = 1, end = Nx - 1; i < end; ++i) {
      if (p[i + 1] >= p[i]) {
        b[i] = k * rho[i + 1] / mu / h / h;
      } else {
        b[i] = k * rho[i] / mu / h / h;
      }
    }
    // a
    for (size_t i = 1, end = Nx - 1; i < end; ++i) {
      a[i] = -c[i - 1] - b[i] - phi * c_f * rho_0 / tau;
    }
    // d
    d[0] = p_inj;
    for (size_t i = 1, end = Nx - 1; i < end; ++i) {
      d[i] = -phi * c_f * rho_0 * p[i] / tau;
    }
    d.back() = p_prod;
    TridiagonalMatrix A(c, a, b);
    p = Sweep(A, d);
    if (counter < num_frames) {
      if (t > record[counter]) {
        writeToFile(p, t);
        ++counter;
      }
    }
  }
  file.close();
  return 0;
}
