#include "ccsd.h"
#include "ccsd_setup.h"
#include "cc_int.h"

#include <iostream>
#include <vector>

using namespace std;


double run_ccsd(const string& sys_in, const string& basis_in, const string& root_in) {

	string sys = sys_in; 
	string basis = basis_in;
	string root = root_in;

	// run HF SCF to get converged MOs
	Result conv_orbs = run_scf(sys, basis);
	printf("\nConverged MO Coefficient Matrix:\n");
	print_matrix(conv_orbs.C);
	printf("\nConverged MO Energies:\n");
	cout << conv_orbs.epsi << endl;
	double Escf = conv_orbs.Escf;
	int nao = conv_orbs.C.rows();

	// AO to MO basis transformation (spatial) 
	vector<double> TEI_MO = ao2mo(conv_orbs.TEI_AO, conv_orbs.C, nao); // N^5 algorithm
	// spatial MO basis -> spin-orbital basis transform
	Tensor4 TEI_SO = mo2so(TEI_MO, nao);
	int nso = 2 * nao; // total number of spin orbitals (occ + vit) 
	int ndocc = conv_orbs.nocc;
	int nelec = 2 * ndocc; // number of occupied spin orbitals
	
	Matrix h = conv_orbs.C.transpose() * conv_orbs.H * conv_orbs.C; // Hamilontian in MO basis
	Matrix f = build_f(h, TEI_SO, nso, nelec); 	// spin-orbital Fock matrix

	Denoms D = build_denominators(f, nso, nelec);
	const Matrix& D1 = D.D1;
	const Tensor4& D2 = D.D2;
	Amplitudes A = initial_guess_amp(TEI_SO, D2, nso, nelec);
	double Emp2 = Emp2_guess(TEI_SO, A, nso, nelec);
	double Etot = Escf + Emp2;
	printf("\nEscf = %10.12f\n", Escf);
	printf("\nEmp2 = %10.12f\n", Emp2);
	printf("\nEtot = %10.12f\n", Etot);

	double Ecc = cc_energy(f, TEI_SO, A, nso, nelec);
	//printf("\nIter = 1  Ecc = %10.12f\n", Ecc);
	double E_old = Ecc;

	for(int iter=1; iter<=100; iter++) {
		// Calculate CC intermediates and effective doubles (DOI: 10.1063/1.460620)
		Tensor4 tau = build_tau(A, nso, nelec);
		Tensor4 tau_tilde = build_tau_tilde(A, nso, nelec);

		// two-index (F) intermediates
		Matrix Fae = build_Fae(f, A, TEI_SO, tau_tilde, nso, nelec);
		Matrix Fmi = build_Fmi(f, A, TEI_SO, tau_tilde, nso, nelec);
		Matrix Fme = build_Fme(f, A, TEI_SO, nso, nelec);

		// four-index (W) intermediates
		Tensor4 Wmnij = build_Wmnij(A, TEI_SO, tau, nso, nelec);
		Tensor4 Wabef = build_Wabef(A, TEI_SO, tau, nso, nelec);
		Tensor4 Wmbej = build_Wmbej(A, TEI_SO, nso, nelec);

		// update cluster amplitudes 
		Amplitudes A_new = update_amp(f, A, D, TEI_SO, Fae, Fmi, Fme, tau, Wmnij, Wabef, Wmbej, nso, nelec);
		A = A_new;
	
		Ecc = cc_energy(f, TEI_SO, A, nso, nelec);
		Etot = Escf + Ecc;
		double dE = Ecc - E_old; 
		//double rms = (D - D_old).norm(); 

		printf("%02d %20.12f %20.12f %20.12f\n", iter, Ecc, Etot, dE); 
		if(fabs(dE) < 1e-12) {
			printf("\nThe CCSD energy has converged!\n");
			break;
		}
		E_old = Ecc;
	}
	return Ecc;
}
