#ifndef GS
#define GS
#include <vector>
#include "../Matrixes/Pentadiagonal.hpp"
#include "../UsefulOperations/VectorOperations.hpp"

std::vector<double> GS_Penta(const PentadiagonalMatrix &A, const std::vector<double> &b, const std::vector<double> &x_0, double tolerance);

#endif