#ifndef CCSD_H
#define CCSD_H

#include "mp2.h"

#include <string>
#include <vector>

using Vec1 = std::vector<double>; 
using Vec2 = std::vector<Vec1>; 
using Vec3 = std::vector<Vec2>;
using Tensor4 = std::vector<Vec3>;

struct Denoms {
	Matrix D1; 
	Tensor4 D2;
};

struct Amplitudes {
	Matrix t1;
	Tensor4 t2;
};

double run_ccsd(const std::string& sys_in, const std::string& basis_in, const std::string& root_in = "../input");

#endif
