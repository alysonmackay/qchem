#ifndef CCSD_SETUP_H
#define CCSD_SETUP_H

#include "ccsd.h"

Matrix build_f(const Matrix& h, const Tensor4& TEI_SO, int nso, int nelec);
Denoms build_denominators(const Matrix& f, int nso, int nelec);
Amplitudes initial_guess_amp(const Tensor4& TEI_SO, const Tensor4& D2, int nso, int nelec);
double Emp2_guess(const Tensor4& TEI_SO, const Amplitudes& A, int nso, int nelec);

#endif
