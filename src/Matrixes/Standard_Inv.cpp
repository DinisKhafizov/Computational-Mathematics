#include "StandardMatrix.hpp"
#include "../SimpleIteration/SI.hpp"

Matrix get_inv(Matrix &A) {
    const int N = A.GetN();
    int counter = 0;
    std::vector<double> b(N), result(N * N), x;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            b[j] = 0;
            if (i == j) {
                b[j] = 1;
            }
        }
        std::vector<double> nulls(N, 1);
        x = SIM(A, nulls, b, 1e-3, 1e-2);
        for (int j = 0; j < N; ++j) {
            nulls[counter + j * N] = x[j];
        }
        ++counter;
    }
    Matrix INV(result, N, N);
    return INV;
}
