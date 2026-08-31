#ifndef TRIPLES_H
#define TRIPLES_H

#include "ccsd_setup.h"

using Vec1 = std::vector<double>;
using Vec2 = std::vector<Vec1>;
using Vec3 = std::vector<Vec2>;
using Vec4 = std::vector<Vec3>;
using Vec5 = std::vector<Vec4>;
using Tensor6 = std::vector<Vec5>; // six-dimensional array for triples storage


Tensor6 build_D3(const Matrix& f, int nso, int nelec);
Tensor6 disconnected_triples(const Amplitudes& A, const Tensor4& TEI_SO, const Tensor6& D3, int nso, int nelec);
Tensor6 connected_triples(const Amplitudes& A, const Tensor4& TEI_SO, const Tensor6& D3, int nso, int nelec); 

double Etriples(const Tensor6& T3c, const Tensor6& T3d, const Tensor6& D3, int nso, int nelec);

#endif
