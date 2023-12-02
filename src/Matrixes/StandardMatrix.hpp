#ifndef SMATRIX
#define SMATRIX

#include <iostream>
#include <vector>


class Matrix {
private:
    std::vector<double> A;
    int M, N;
public:
    Matrix(const int m, const int n): A{std::vector<double>(m* n)}, M{m}, N{n} {}
    Matrix(const std::vector<double> &a, const int m, const int n): A{a}, N{n}, M{m} {}
    Matrix(): N{0}, M{0} {}

    double operator ()(int i, int j) const;
    double &operator ()(int i, int j);
    void operator*(const double x);
    Matrix operator*(const Matrix &B); //for tests only

    std::vector<double> operator*(const std::vector<double> &x) const;
    std::vector<double> partly_dot(const std::vector<double> &x);

    void change_Col(const int j, const int i_begin, const std::vector<double> &x);
    void add_Identity();
    std::vector<double> getRow(const int i)  ;

    std::vector<double> getCol(const int j, const int i_begin = 0);
    std::vector<double> getRow(const int i, const int j_end = 0)  const;

    std::vector<double> getCol(const int j) const;
    std::vector<double> get_vals() const;


    void transpose();


    Matrix transpose(const int what);

    void SetStringOnEnd(const std::vector<double> &newstr);

    void deleteStringOnEnd();

    int &GetN();
    int &GetM();
    int GetN() const;
    int GetM() const;

    std::vector<double> &get_vals();

    void show();
};

#endif