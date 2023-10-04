#include "Tridiagonal.hpp"

int TridiagonalMatrix::GetN() { return n; }

double TridiagonalMatrix::GetA(int i) { return A[i]; }
double TridiagonalMatrix::GetB(int i) { return B[i]; }
double TridiagonalMatrix::GetC(int i) { return C[i]; }
