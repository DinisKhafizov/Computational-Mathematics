#ifndef UGOLOK
#define UGOLOK

#include "../UsefulOperations/VectorOperations.hpp"
#include "../Matrixes/StandardMatrix.hpp"
#include <vector>
#include <cmath>

std::vector<double> laks_vendroff_task(const double (*y)(const double), double h, double tau, double L, double c, double time, std::vector<double> &TIMES);
std::vector<double> ugolok_left_task(const double (*y)(const double), double h, double tau, double L, double c, double time, std::vector<double> &TIMES);

std::vector<double> laks_vendroff(const double (*y)(const double), double h, double tau, double L, double c, double time);
std::vector<double> ugolok_left(const double (*y)(const double), double h, double tau, double L, double c, double time);
std::vector<double> hybrid(const double (*y)(const double), double h, double tau, double L, double c, double time);

#endif