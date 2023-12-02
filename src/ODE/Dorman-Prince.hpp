#ifndef DORMAN_PRINCE
#define DORMAN_PRINCE

#include <vector>
#include "../UsefulOperations/VectorOperations.hpp"
#include "../Matrixes/StandardMatrix.hpp"

std::pair<std::vector<double>, std::vector<double>> dorman_prince(const double a, const double b, const std::vector<double> &u_0,
    std::vector<double> (*f)(const double, const std::vector<double> &),
    const double epsilon, const double h);

#endif