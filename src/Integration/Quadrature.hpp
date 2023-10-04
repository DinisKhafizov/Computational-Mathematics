#ifndef QUADRATURE_FORMULAS
#define QUADRATURE_FORMULAS
#include <vector>
#include <iostream>
#include <cmath>
const double trapezoid(const double a, const double b, const double (*y)(const double));
const double simpson(const double a, const double b, const double (*y)(const double));
const double three_eights(const double a, const double b, const double (*y)(const double));
const double gauss_quadrature(const double a, const double b, const size_t deg, const double (*y)(const double));
#endif