#include "SI.hpp"

std::vector<double> SIM(const Matrix &A, const std::vector<double> &x_0, const std::vector<double> &b, 
const double tolerance, const double tau) {

	std::vector<double> x = x_0, r = A * x_0 - b;

	while(first_norm(r) > tolerance) {
        x = x - tau*r;
		r = A * x;
		r = r - b;
	}
	return x;
}