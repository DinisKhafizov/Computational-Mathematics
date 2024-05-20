#ifndef TRIDIAG
#define TRIDIAG
#include <vector>

class TridiagonalMatrix {
   private:
    std::vector<double> A, B, C;  // A - under diag, В - diag, С - up, n -dim
    int n;

   public:
    TridiagonalMatrix(const std::vector<double> &A1,
                      const std::vector<double> &B1,
                      const std::vector<double> &C1)
        : A{A1}, B{B1}, C{C1}, n{static_cast<int>(size(B1))} {}

    int GetN();
    double GetA(int i);
    double GetB(int i);
    double GetC(int i);
};

std::vector<double> Sweep(TridiagonalMatrix &Matrix,
                          const std::vector<double> &D);

#endif