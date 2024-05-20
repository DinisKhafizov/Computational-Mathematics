#ifndef PENTA
#define PENTA

#include <vector>

class PentadiagonalMatrix {
   private:
    std::vector<double> A1, A2, B, C1, C2;  // A - under diag, В - diag, С - up, n -dim
    int n;

   public:
    PentadiagonalMatrix(const std::vector<double> &A1_, const std::vector<double> &A2_,
                      const std::vector<double> &B1_,
                      const std::vector<double> &C1_, const std::vector<double> &C2_)
        : A1{A1_}, A2{A2_}, B{B1_}, C1{C1_}, C2{C2_}, n{static_cast<int>(size(B1_))} {}

    int GetN();
    double getA1(int i);
    double getA2(int i);
    double getB(int i);
    double getC1(int i);
    double getC2(int i);
};

#endif