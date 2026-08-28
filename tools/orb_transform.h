#ifndef ORB_TRANSFORM_H
#define ORB_TRANSFORM_H 

#include "hf.h"
#include <vector>

using Vec1 = std::vector<double>; 
using Vec2 = std::vector<Vec1>; 
using Vec3 = std::vector<Vec2>;
using Tensor4 = std::vector<Vec3>;

int spatial(int k);
bool same_spin(int k, int l);

std::vector<double> ao2mo(const std::vector<double>& TEI_AO, const Matrix& coeff, int nao);
std::vector<double> ao2mo_brute(const std::vector<double>& TEI_AO, const Matrix& coeff, int nao);
Tensor4 mo2so(const std::vector<double>& TEI_MO, int nao);


#endif 
