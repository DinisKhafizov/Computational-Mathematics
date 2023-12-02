#ifndef NEWTON_EQATION
#define NEWTON_EQUATION
#include "UsefulOperations/VectorOperations.hpp"

const double newton_eq(const double (*y)(const double), const double (*y_der)(const double), const double x_0, const double tolerance);

#endif