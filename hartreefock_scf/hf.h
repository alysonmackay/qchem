#ifndef HF_H
#define HF_H

#include <string>
#include <vector>

#include "Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

struct Result {
	Matrix C;
	Vector epsi;
	int nocc;
	double Escf;
	std::vector<double> TEI_AO;
};

extern std::vector<int> ioff;
int compound(int a, int b); 

void print_matrix(const Matrix& mat);

Result run_scf(const std::string& sys_in, const std::string& basis_in, const std::string& root_in = "../input");

#endif
