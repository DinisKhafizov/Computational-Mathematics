#include "VectorOperations.hpp"

std::vector<double> operator+(const std::vector<double> &a,
                              const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}
std::vector<double> operator+(const std::vector<double> &a, double b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] + b;
    }
    return res;
}

std::vector<double> operator+=(const std::vector<double> &a,
                               const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}

std::vector<double> operator-(const std::vector<double> &a,
                              const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] - b[i];
    }
    return res;
}
std::vector<double> operator-(const std::vector<double> &a, double b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] + b;
    }
    return res;
}

std::vector<double> operator-=(const std::vector<double> &a,
                               const std::vector<double> &b) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] - b[i];
    }
    return res;
}

double operator*(const std::vector<double> &a, const std::vector<double> &b) {
    double res = 0;
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res += a[i] * b[i];
    }
    return res;
}

std::vector<double> operator*(double x, const std::vector<double> &a) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] * x;
    }
    return res;
}

std::vector<double> operator*(const std::vector<double> &a, double x) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] * x;
    }
    return res;
}

std::vector<double> operator/(const std::vector<double> &a, double x) {
    std::vector<double> res(size(a));
    for (size_t i = 0, end = size(a); i < end; ++i) {
        res[i] = a[i] / x;
    }
    return res;
}
std::vector<double> operator/(const std::vector<double> &a,
                              const std::vector<double> &b) {
    const int n = size(a);
    std::vector<double> res(n);
    if (n != size(b)) {
        std::cout << "Dimensions of vectors must be equal! ('vec1/vec2' - "
                     "element-wise division)"
                  << std::endl;
        return res;
    }
    for (size_t i = 0; i < n; ++i) {
        res[i] = a[i] / b[i];
    }
    return res;
}

std::vector<double> operator/=(const std::vector<double> &x, const double a) {
    std::vector<double> res(size(x));
    for (size_t i = 0, end = size(x); i < end; ++i) {
        res[i] = x[i] / a;
    }
    return res;
}

bool operator==(const std::vector<double> &a, std::vector<double> &b) {
    if (size(a) != size(b)) {
        return false;
    }
    for (size_t i = 0, end = size(a); i < end; ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

// norms: first, second, endless

double first_norm(const std::vector<double> &x) {
    double Norm = 0;
    for (size_t i = 0, end = size(x); i < end; ++i) {
        Norm += std::abs(x[i]);
    }
    return Norm;
}

double second_norm(const std::vector<double> &x) { return sqrt(x * x); }

double endless_norm(const std::vector<double> &x) {
    double max = 0;
    for (size_t i = 0, end = size(x) - 1; i < end; ++i) {
        if (std::abs(x[i]) > max) {
            max = std::abs(x[i]);
        }
    }
    return max;
}

// specific

std::vector<double> elWise_Mult(
    const std::vector<double> &a,
    const std::vector<double>
        &b) {  // vect1 * vect2 == vect3, where vect3(i) == vect1(i) * vect2(i)
    const int n = size(a);
    std::vector<double> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = a[i] * b[i];
    }
    return res;
}

std::vector<double> linspace(const double start, const double end,
                             const int points_number) {
    const double diff = (end - start) / (points_number - 1);
    std::vector<double> res(points_number);
    for (int i = 0; i < points_number; ++i) {
        res[i] = start + diff * i;
    }
    return res;
}

std::vector<double> get_vectors_part(const std::vector<double> &x, const int begin, const int end) {
    std::vector<double> res(end - begin);
    for (size_t i = 0, end1 = size(res); i < end1; ++i) {
        res[i] = x[begin + i];
    }
    return res;
}
/*
Помещает вектор v1 в вектор v2 с индекса start
v1 = [1, 2]
v2 = [0, 0, 0, 0, 0]
putV1ToV2(v1, v2, 1) -> v2 = [0, 1, 2, 0, 0]
*/
void putV1ToV2(const std::vector<double> &v1, std::vector<double> &v2, int start) {
    for (size_t i = 0, end = size(v1); i < end; ++i) {
        v2[start + i] = v1[i];
    }
}