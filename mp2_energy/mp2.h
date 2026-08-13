#ifndef MP2_H
#define MP2_H 

#include "hf.h"

#include <string>
#include <vector>

std::vector<double> ao2mo(const std::vector<double>& TEI_AO, const Matrix& coeff, int nao);

double mp2_energy(const std::string& sys_in, const std::string& basis_in, const std::string& root_in = "../input");

#endif
