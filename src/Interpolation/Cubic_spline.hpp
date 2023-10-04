#ifndef CUBIC_SPLINE
#define CUBIC_SPLNE
#include "../Matrixes/Tridiagonal.hpp"

std::vector<double> Sweep(TridiagonalMatrix &Matrix,
                          const std::vector<double> &D);

std::vector<double> Interpolation(std::vector<double> &x,
                                  std::vector<double> &y, double sec_der_begin,
                                  double sec_der_end);

void WritingToFile(std::vector<double> res);

void Interpolate_and_write(std::vector<double> &x, std::vector<double> &y,
                           double sec_der_begin, double sec_der_end);
std::vector<double> cubic_spline_interopolation(const std::vector<double> &x, const std::vector<double> &y);

#endif