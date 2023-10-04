#include "Lagrange.hpp"

Lagrange_interpolation::Lagrange_interpolation(const std::vector<double> &x0,
                                               size_t N0) {
    x = x0;
    N = N0;
    product.resize(N0, 1);
    for (size_t i = 0; i < N0; ++i) {
        for (size_t j = 0; j < N0; ++j) {
            if (i == j) {
                continue;
            }
            product[i] *= (x[i] - x[j]);
        }
    }
}

double Lagrange_interpolation::get_Y_val(double x_val,
                                         const std::vector<double> &y_points) {
    double y = 0, prod;
    for (size_t i = 0; i < N; ++i) {
        prod = 1;
        for (size_t j = 0; j < N; ++j) {
            if (i == j) {
                continue;
            }
            prod *= (x_val - x[j]);
        }
        y += y_points[i] * prod / product[i];
    }
    return y;
}

void Lagrange_interpolation::calc_and_write_mesh(
    std::string filename, const std::vector<double> &y_points,
    const int number_of_points) {
    std::ofstream fout;
    fout.open(filename);
    std::vector<double> x_mesh = linspace(x[0], x.back(), number_of_points),
                        y_mesh(number_of_points);
    for (size_t i = 0, end = N - 1; i < end; ++i) {
        fout << x[i] << ";";
    }
    fout << x[N - 1] << std::endl;
    for (size_t i = 0, end = N - 1; i < end; ++i) {
        fout << y_points[i] << ";";
    }
    fout << y_points[N - 1] << std::endl;

    for (int i = 0, end = number_of_points - 1; i < end; ++i) {
        fout << x_mesh[i] << ";";
    }
    fout << x_mesh[number_of_points - 1] << std::endl;

    for (int i = 0, end = number_of_points - 1; i < end; ++i) {
        fout << get_Y_val(x_mesh[i], y_points) << ";";
    }
    fout << get_Y_val(x_mesh[number_of_points - 1], y_points);
}