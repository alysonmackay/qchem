#ifndef CCSD_H
#define CCSD_H

#include "orb_transform.h"

#include <string>

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
