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

struct ccsd {
	Matrix f;
	Amplitudes A;
	Tensor4 TEI_SO;
	int nso;
	int nelec;
	double Eccsd;
	double Ecc;
};

ccsd run_ccsd(const std::string& sys_in, const std::string& basis_in, const std::string& root_in = "../input");

#endif
