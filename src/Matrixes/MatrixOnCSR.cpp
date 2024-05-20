#include "MatrixOnCSR.hpp"

double CSR::operator()(int i, int j) const { 
	for(int k = rows[i]; k < rows[i + 1]; ++k) {
		if (cols[k] == j) {
			return values[k];
		}
	}
	return 0;
}

std::vector<double> CSR::operator*(const std::vector<double> &vec) const {
	std::vector<double> res(ROWS);	
	for (int i = 0; i < ROWS; ++i) {
		for (int j = rows[i]; j < rows[i + 1]; ++j) {
			res[i] += values[j] * vec[cols[j]];
		}
	}
	return res;
}
void CSR::operator *(double x) {
	for (int i = 0; i < size(values); ++i) {
		values[i] *= x;
	}
}

std::vector<double> CSR::getDiag() const { //only for square matrix!
	std::vector<double> res(COLS);
	for (size_t i = 0; i < COLS; ++i) {
		for (size_t k = rows[i]; k < rows[i + 1]; ++k) {
			if (cols[k] == i) { 
				res[i] = values[k];
			}
		}
	}
	return res;
}

CSR CSR::transpose() const {
        size_t Nonzero = size(values);
        std::vector<double> tmatrix_els(Nonzero);
        std::vector<int> tCols(Nonzero);
        std::vector<int> tRows(COLS + 1);
        for (size_t i = 0; i < Nonzero; ++i) tRows[cols[i] + 1]++;
        size_t S = 0;
        size_t tmp;
        for (size_t i = 1; i <= COLS; ++i) {
            tmp = tRows[i];
            tRows[i] = S;
            S += tmp;
        }
        size_t j1, j2, Col, RIndex, IIndex;
        double V;
        for (size_t i = 0; i < ROWS; ++i) {
            j1 = rows[i];
            j2 = rows[i+1];
            Col = i;
            for (size_t j = j1; j < j2; ++j) {
                V = values[j];
                RIndex = cols[j];
                IIndex = tRows[RIndex + 1];
                tmatrix_els[IIndex] = V;
                tCols[IIndex] = Col;
                tRows[RIndex + 1]++;
            }
        }
        return CSR(tmatrix_els, tCols, tRows, ROWS);
}

int CSR::GetN() const {
	return ROWS;
}
std::vector<double> &CSR::GetVal() {
	return values;
}
std::vector<int> &CSR::GetRow()  {
	return rows;
}
std::vector<int> &CSR::GetCol() {
	return cols;
}
std::vector<double> CSR::GetVal() const {
	return values;
}
std::vector<int> CSR::GetRow()  const{
	return rows;
}
std::vector<int> CSR::GetCol()  const{
	return cols;
}
void CSR::show() { //only for test
	CSR RES(values, cols, rows, COLS);
	for(size_t i = 0; i < ROWS; ++i) {
		for(size_t j = 0; j < COLS; ++j) {
			std::cout << RES(i, j) << " ";
		}
		std::cout << std::endl;
	}
}
	