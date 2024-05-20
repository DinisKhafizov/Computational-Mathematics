#include "../src/2sem/funcs.hpp"
#include <fstream>

const double y(double x) {
    double result = 1;
    if (x < 0.5) {
        return 10;
    }
    return result;
}

void writeToFile(std::string filename, const std::vector<double> &x, int type) {
    std::ofstream file;
    if (type == 0) {
        file.open(filename);
    }
    else {
        file.open(filename, std::ofstream::app);
    }
    for (size_t i = 0, end = size(x) - 1; i < end; ++i) {
        file << x[i] << ";";
    }
    file << x.back() << "\n";
}

int main() {
    const double L = 1;
    const double c = 1;
    const double T = 1./3;
    const double h = 1e-2;
    const double tau = 5e-3;
    std::vector<double> ugolok = ugolok_left(y, h, tau, L, c, T);
    std::vector<double> laks = laks_vendroff(y, h, tau, L, c, T);
    std::vector<double> hybri = hybrid(y, h, tau, L, c, T);
    writeToFile("daywork.txt", ugolok, 0);
    writeToFile("daywork.txt", laks, 1);
    writeToFile("daywork.txt", hybri, 1);
    return 0;
}
