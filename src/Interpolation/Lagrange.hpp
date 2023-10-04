#ifndef LAGR
#define LAGR
#include <fstream>

#include "../UsefulOperations/VectorOperations.hpp"

class Lagrange_interpolation {
   private:
    size_t N;
    std::vector<double> x, product;

   public:
    Lagrange_interpolation(const std::vector<double> &x0, size_t N0);
    double get_Y_val(double x_val, const std::vector<double> &y_points);
    void calc_and_write_mesh(std::string filename,
                             const std::vector<double> &y_points,
                             const int number_of_points);
};

#endif