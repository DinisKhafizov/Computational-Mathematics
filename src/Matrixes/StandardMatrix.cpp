#include "StandardMatrix.hpp"



double Matrix::operator ()(int i, int j) const {
    return A[i * N + j];
}
double &Matrix::operator ()(int i, int j) {
    return A[i * N + j];
}

std::vector<double> Matrix::operator*(const std::vector<double> &x) const {
    std::vector<double> res(M);
    int k;
    for (size_t i = 0; i < M; ++i) {
        k = i * N;
        for (size_t j = 0; j < N; ++j) {
            res[i] += A[k + j] * x[j];
        }
    }
    return res;
}
Matrix Matrix::operator*(const Matrix &B) { //usable only for tests
    Matrix Res(0, B.GetN());
    std::vector<double> string_end(B.GetN()), row(N);
    double lol=0;
    for (size_t i = 0; i < M; ++i) {
        lol = 0;
        for (size_t k = 0; k < N; ++k) {
            row[k] = A[i*N + k];
        }
        for(size_t j = 0, end = B.GetN(); j < end; ++j) {
            lol = 0;
            for(size_t l = 0; l < size(row); ++l) {
                lol += row[l] * B.getCol(j)[l];
            }
            string_end[j] = lol;
        }
        Res.SetStringOnEnd(string_end);
    }
    return Res;
}


void Matrix::operator*(const double x) {
    for (size_t i = 0, end = size(A); i < end; ++i) {
        A[i] *= x;
    }
}


std::vector<double> Matrix::partly_dot(const std::vector<double> &x) {
    const int N_x = size(x);
    int k;
    std::vector<double> res(M);
    for (size_t i = 0; i < M; ++i) {
        k = i * N;
        for (size_t j = 0; j < N_x; ++j) {
            res[i] += A[k + j] * x[j]; 
        }
    }
    return res;

}

void Matrix::change_Col(const int j, const int i_begin, const std::vector<double> &x) {
    for (size_t i = i_begin; i < M; ++i) {
        A[j + i*N] = x[i - i_begin];
    }
}
void Matrix::add_Identity() {
    const int k = N + 1;
    for (size_t i = 0; i < N; ++i) {
        A[k * i] += 1;
    }
}

std::vector<double> Matrix::getRow(const int i) {
    std::vector<double> vec(N);
    const int k = i * N;
    for (int j = 0; j < N; ++j) {
        vec[j] = A[k + j];
    }
    return vec;
}

std::vector<double> Matrix::getCol(const int j, const int i_begin)  {
    std::vector<double> vec(M - i_begin);
    for (size_t i = i_begin; i < M; ++i) {
        vec[i - i_begin] = A[j + i*N];
    }
    return vec;
}
std::vector<double> Matrix::getRow(const int i, const int j_end) const {
    std::vector<double> vec(j_end);
    const int k = i * N;
    if (j_end == 0) {
        vec.resize(N);
        for (int j = 0; j < N; ++j) {
            vec[j] = A[k + j];
        }
        return vec;
    }
    for (int j = 0; j < j_end; ++j) {
        vec[j] = A[k + j];
    }
    return vec;
}

std::vector<double> Matrix::getCol(const int j) const {
    std::vector<double> vec(M);
    for (size_t i = 0; i < M; ++i) {
        vec[i] = A[j + i*N];
    }
    return vec;
}


void Matrix::transpose() {
    std::vector<double> res(M * N);
    int k;
    for (size_t i = 0; i < N; ++i) {
        k = i * M;
        for (size_t j = 0; j < M; ++j) {
            res[j + k] = A[i + j*N];
        }
    }
    A = res;
    k = N;
    N = M;
    M = k;
}


Matrix Matrix::transpose(const int what){
    std::vector<double> res(M * N);
    int k;
    for (size_t i = 0; i < N; ++i) {
        k = i * M;
        for (size_t j = 0; j < M; ++j) {
            res[j + k] = A[i + j*N];
        }
    }
    Matrix B(res, N, M);
    return B;
}

void Matrix::SetStringOnEnd(const std::vector<double> &newstr) {
    const int n = size(newstr), m = size(A);
    A.resize(m + n);
    for (int i = 0; i < n; ++i) {
        A[m + i] = newstr[i];
    }
    if (M == 0) {
        M = 1;
    }
    else {
        M += 1;
    }
    if (N == 0) {
        N = n;
    }
}

void Matrix::deleteStringOnEnd() {
    for (size_t i = 0; i < N; ++i) {
        A.pop_back();
    }
    M -= 1;
}

int &Matrix::GetN() {
    return(N);
}

int &Matrix::GetM() {
    return(M);
}
int Matrix::GetN() const {
    return(N);
}

int Matrix::GetM() const{
    return(M);
}

std::vector<double> &Matrix::get_vals() {
    return A;
}
std::vector<double> Matrix::get_vals() const {
    return A;
}

void Matrix::show() {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << A[i*N + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
