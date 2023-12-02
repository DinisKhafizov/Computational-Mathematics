#ifndef VECTOR_OP
#define VECTOR_OP

#include <cmath>
#include <iostream>
#include <vector>

/*
Three sections:
1) Arithmetic operations
2) Norms counting
3) Specific operations with vector

*/

// Arithmetic: '+', '-', '*', '/', '='

std::vector<double> operator+(const std::vector<double> &a,
                              const std::vector<double> &b);
std::vector<double> operator+(const std::vector<double> &a, double b);
std::vector<double> operator+=(const std::vector<double> &a,
                               const std::vector<double> &b);

std::vector<double> operator-(const std::vector<double> &a,
                              const std::vector<double> &b);
std::vector<double> operator-(const std::vector<double> &a, double b);
std::vector<double> operator-=(const std::vector<double> &a,
                               const std::vector<double> &b);

double operator*(const std::vector<double> &a, const std::vector<double> &b);
std::vector<double> operator*(double x, const std::vector<double> &a);
std::vector<double> operator*(const std::vector<double> &a, double x);

std::vector<double> operator/(const std::vector<double> &a, double x);
std::vector<double> operator/(const std::vector<double> &a,
                              const std::vector<double> &b);
std::vector<double> operator/=(const std::vector<double> &x, const double a);

bool operator==(const std::vector<double> &a, std::vector<double> &b);

// norms: first, second, endless

double first_norm(const std::vector<double> &x);

double second_norm(const std::vector<double> &x);

double endless_norm(const std::vector<double> &x);

// specific

std::vector<double> elWise_Mult(const std::vector<double> &a,
                                const std::vector<double> &b);

std::vector<double> linspace(const double start, const double end,
                             const int points_number);
std::vector<double> get_vectors_part(const std::vector<double> &x, const int begin, const int end) ;
#endif