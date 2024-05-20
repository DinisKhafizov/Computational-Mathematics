#include "../src/Matrixes/StandardMatrix.hpp"
#include "../src/UsefulOperations/VectorOperations.hpp"
#include <cmath>
#include <fstream>

void WriteToFile(Matrix &A, double gamma);
double max_val_diag(Matrix &A);

int main() {

    // Константы 
    const double L = 10, gamma = 5. / 3, v_l = 0, v_r = 0, rho_l = 13, rho_r = 1.3, P_l = 10e5, P_r = 1e5;
    const double T = 0.02;
    const double h = 2e-1;
    double tau = 1e-7;

    // Переменные
    double t = tau, u, e, c, rho, max_lambda, CFL;
    Matrix W_1(3, static_cast<int>(2. * L/h)), W_2(3, static_cast<int>(2. * L/h));

    // Начальные условия
    for (int i = 0, end = W_1.GetN() / 2; i < end; ++i) {
        W_1(0, i) = rho_l;
        W_1(1, i) = 0;
        W_1(2, i) = P_l / (gamma - 1);
    }
    for (int i = W_1.GetN() / 2, end = W_1.GetN(); i < end; ++i) {
        W_1(0, i) = rho_r;
        W_1(1, i) = 0;
        W_1(2, i) = P_r / (gamma - 1);
    }
    // Основной цикл
    for (; t < T; t += tau) {
        for (int i = 1, end = W_1.GetN() - 1; i < end; ++i) {
            rho = W_1(0, i);        
            u = W_1(1, i) / rho;
            e = W_1(2, i) / rho;
            c = sqrt(gamma * (gamma - 1) * e);
            Matrix A_all({0, 1, 0, 
                          -u * u, 2 * u, gamma - 1, 
                          -gamma * u * e, gamma * e, u}, 3, 3),
                   Lambda_abs({fabs(u + c), 0, 0, 
                               0, fabs(u), 0, 
                               0, 0, fabs(u - c)}, 3, 3),
                   Omega_T({-u * c, c, gamma - 1, 
                            -c * c, 0, gamma - 1, 
                            u * c, -c, gamma - 1}, 3, 3), 
                   Omega_T_inv({1 / (2 * c * c), -1 / (c * c), 1/(2 * c * c), 
                                (c + u)/(2 * c * c), -u/(c * c), (-c + u)/(2 * c * c), 
                                1/(2 * gamma - 2), 0, 1/(2 * gamma - 2)}, 3, 3);
            max_lambda = max_val_diag(Lambda_abs);
            CFL = tau * max_lambda / h;
            if (CFL > 1e-2) {
                tau = 1e-2 * h / max_lambda;
            }
            std::vector<double> temp_i = W_1.getCol(i), temp_im1 = W_1.getCol(i - 1), temp_ip1 = W_1.getCol(i + 1);
            std::vector<double> c1 = A_all * (temp_ip1 - temp_im1);
            std::vector<double> c2 = Omega_T_inv * Lambda_abs * Omega_T * (temp_ip1 - 2 * temp_i + temp_im1);
            std::vector<double> new_col = temp_i - (tau / (2 * h)) * (c1 - c2);
            W_2.setCol(i, new_col);
        }
        W_2.setCol(0, W_2.getCol(1));
        W_2.setCol(W_2.GetN() - 1, W_2.getCol(W_2.GetN() - 2));
        W_1 = W_2;
        if (t >= 0.015) {
            break;
        }
    }
    WriteToFile(W_2, gamma);
    return 0;
}

void WriteToFile(Matrix &A, double gamma) {
    std::ofstream file("../../tasks/hw1/second_task.txt");
    for (int i = 0, end = A.GetN(); i < end; ++i) {
        file << A(0, i) << ";";
    }
    file << "\n";
    for (int i = 0, end = A.GetN(); i < end; ++i) {
        file << A(1, i) / A(0, i) << ";";
    }
    file << "\n";
    for (int i = 0, end = A.GetN(); i < end; ++i) {
        file << A(2, i) / A(0, i) << ";";
    }
    file << "\n";
    for (int i = 0, end = A.GetN(); i < end; ++i) {
        file << (gamma - 1) * A(2, i) << ";";
    }
    file << "\n";
}

double max_val_diag(Matrix &A) {
    double res = A(0, 0);
    if (res < A(1, 1)) {
        res = A(1, 1);
    }
    if (res < A(2, 2)) {
        res = A(2, 2);
    }
    return res;
}