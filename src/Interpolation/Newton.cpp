#include "Newton.hpp"

std::vector<double> newton_interpolation(std::vector<double> x,
                                         std::vector<double> y) {
    const size_t N = size(x);  // number of points
    std::vector<double> res(N);
    res[0] = y[0];
    for (size_t i = 0, end = N - 1; i < end; ++i) {
        y[0] = (y[1] - y[0]) / (x[1 + i] - x[0]);
        res[i + 1] = y[0];
        for (size_t j = 1, end = N - 1 - i; j < end; ++j) {
            y[j] = (y[j + 1] - y[j]) / (x[j + 1 + i] - x[j]);
        }
    }
    return res;
}

double get_polynoms_value(double x, const std::vector<double> &res,
                          const std::vector<double> &x_vals) {
    const size_t N = size(x_vals);
    double product = x - x_vals[0], y = res[0];
    for (size_t i = 1; i < N; ++i) {
        y += product * res[i];
        product *= (x - x_vals[i]);
    }
    return y;
}

void get_and_write_Yvalues(std::string filename, const std::vector<double> &x,
                           const std::vector<double> &y,
                           const std::vector<double> &res,
                           const int number_of_points) {
    std::ofstream fout;
    fout.open(filename);

    // Putting (x, y) scatter to file
    for (int i = 0, end = static_cast<int>(size(x)) - 1; i < end; ++i) {
        fout << x[i] << ";";
    }
    fout << x[size(x) - 1] << std::endl;

    for (int i = 0, end = static_cast<int>(size(x)) - 1; i < end; ++i) {
        fout << y[i] << ";";
    }
    fout << y[size(x) - 1] << std::endl;

    // Creating mesh, calculating y, pushing all to the file
    std::vector<double> x_mesh = linspace(x[0], x.back(), number_of_points);
    for (int i = 0, end = number_of_points - 1; i < end; ++i) {
        fout << x_mesh[i] << ";";
    }
    fout << x_mesh[number_of_points - 1] << std::endl;
    for (int i = 0, end = number_of_points - 1; i < end; ++i) {
        fout << get_polynoms_value(x_mesh[i], res, x) << ";";
    }
    fout << get_polynoms_value(x_mesh[number_of_points - 1], res, x);
}
