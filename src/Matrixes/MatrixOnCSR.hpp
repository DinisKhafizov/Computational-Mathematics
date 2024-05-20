#ifndef CSRMATRIX
#define CSRMATRIX

#include<iostream>
#include<vector>

class CSR {
	int ROWS;
	int COLS;
	std::vector<double> values;
	std::vector<int> cols;
	std::vector<int> rows;
public:

	CSR(const int n, const int m): ROWS{m}, COLS{n}, values(m * n), cols(0), rows(1, 0) {} 
	CSR(const std::vector<double> &val, const std::vector<int> &col, const std::vector<int> &row, const int Ncols):values{val}, 
	cols{col}, rows{row}, COLS{Ncols}, ROWS{static_cast<int>(size(row)) - 1} {}

	double operator()(int i, int j) const;

	std::vector<double> operator*(const std::vector<double> &vec) const;
	void operator *(double x);

	std::vector<double> getDiag() const;
	CSR transpose() const;

	int GetN() const;
	std::vector<double> &GetVal();
	std::vector<int> &GetRow();
	std::vector<int> &GetCol();
	std::vector<double> GetVal() const;
	std::vector<int> GetRow()  const;
	std::vector<int> GetCol()  const;

	void show();
	
};
#endif //CSRMATRIX