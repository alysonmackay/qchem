#ifndef CC_INT_H
#define CC_INT_H

#include "ccsd.h"

// effective doubles
Tensor4 build_tau(const Amplitudes& A, int nso, int nelec);
Tensor4 build_tau_tilde(const Amplitudes& A, int nso, int nelec);

// CC intermediates
// two index (F)
Matrix build_Fae(const Matrix& fock, const Amplitudes& A, const Tensor4& TEI_SO, const Tensor4& tau_tilde, int nso, int nelec);
Matrix build_Fmi(const Matrix& fock, const Amplitudes& A, const Tensor4& TEI_SO, const Tensor4& tau_tilde, int nso, int nelec);
Matrix build_Fme(const Matrix& fock, const Amplitudes& A, const Tensor4& TEI_SO, int nso, int nelec);

// four index (W) 
Tensor4 build_Wmnij(const Amplitudes& A, const Tensor4& TEI_SO, const Tensor4& tau, int nso, int nelec);
Tensor4 build_Wabef(const Amplitudes& A, const Tensor4& TEI_SO, const Tensor4& tau, int nso, int nelec);
Tensor4 build_Wmbej(const Amplitudes& A, const Tensor4& TEI_SO, int nso, int nelec); 

// update cluster amplitudes
Amplitudes update_amp(const Matrix& fock, Amplitudes& init, const Denoms& D, const Tensor4& TEI_SO, const Matrix& Fae, const Matrix& Fmi, const Matrix& Fme, const Tensor4& tau, const Tensor4& Wmnij, const Tensor4& Wabef, const Tensor4& Wmbej, int nso, int nelec);

//calculate CCSD energy
double cc_energy(const Matrix& f, const Tensor4& TEI_SO, const Amplitudes& A, int nso, int nelec);

#endif 
