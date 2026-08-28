#ifndef CCSD_SETUP_H
#define CCSD_SETUP_H

#include "ccsd.h"

int spatial(int k);
bool same_spin(int k, int l);

Tensor4 mo2so(const std::vector<double>& TEI_MO, int nao);
Matrix build_f(const Matrix& h, const Tensor4& TEI_SO, int nso, int nelec);
Denoms build_denominators(const Matrix& f, int nso, int nelec);
Amplitudes initial_guess_amp(const Tensor4& TEI_SO, const Tensor4& D2, int nso, int nelec);
double Emp2_guess(const Tensor4& TEI_SO, const Amplitudes& A, int nso, int nelec);

#endif
