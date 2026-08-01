#ifndef HF_H
#define HF_H

#include <string>
#include "Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

struct Result {
	Matrix C;
	Vector epsi;
};

Result run_scf(const std::string& sys_in, const std::string& basis_in, const std::string& root_in = "../input");

#endif
