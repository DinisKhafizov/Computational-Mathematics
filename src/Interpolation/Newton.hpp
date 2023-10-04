#ifndef NEWTON
#define NEWTON
#include <fstream>

#include "../UsefulOperations/VectorOperations.hpp"

std::vector<double> newton_interpolation(std::vector<double> x,
                                         std::vector<double> y);

double get_polynoms_value(double x, const std::vector<double> &res,
                          const std::vector<double> &x_vals);

void get_and_write_Yvalues(std::string filename, const std::vector<double> &x,
                           const std::vector<double> &y,
                           const std::vector<double> &res,
                           const int number_of_points);

#endif