#ifndef SIMPLE_ITERATION_METH
#define SIMPLE_ITERATION_METH

#include <iostream>
#include <vector>
#include "../Matrixes/StandardMatrix.hpp"
#include "../UsefulOperations/VectorOperations.hpp"

std::vector<double> SIM(const Matrix &A, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const double tau);

#endif