#ifndef RUNGE_KUTTA_METHODS
#define RUNGE_KUTTA_METHODS

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "UsefulOperations/VectorOperations.hpp"

std::vector<double> runge_kutta_4(const std::vector<double> &t,
                                  double (*f)(const double, const double),
                                  const double u0, const double h);

std::vector<double> runge_kutta_4_sys(
    const std::vector<double> &t,
    std::vector<double> (*f)(const double, const std::vector<double> &),
    const std::vector<double> &u_0, const double h);

void WritingToFile_sys(const std::vector<double> &res, std::string filename,
                       int type);

#endif