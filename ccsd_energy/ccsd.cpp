#include "ccsd.h"

#include <iostream>
#include <vector>

using namespace std;

int spatial(int k) {return k / 2; }
bool same_spin(int k, int l) {return k%2 == l%2; } // spin parity

Tensor4 mo2so(const vector<double>& TEI_MO, int nao) {
	int nso = 2 * nao; // spin orbitals 
	Tensor4 TEI_SO(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));

	for(int p=0; p < nso; p++)
		for(int q=0; q < nso; q++)
			for(int r=0; r < nso; r++) 
				for(int s=0; s < nso; s++) {
					int pr = compound(spatial(p),spatial(r));
					int qs = compound(spatial(q),spatial(s));
					int pqrs = compound(pr,qs);
					double direct = TEI_MO[pqrs] * same_spin(p,r) * same_spin(q,s); // Dirac <pq|rs>
					int ps = compound(spatial(p),spatial(s));
					int qr = compound(spatial(q),spatial(r));
					int psqr = compound(ps,qr);
					double exchange = TEI_MO[psqr] * same_spin(p,s) * same_spin(q,r); // <pq|sr>
					TEI_SO[p][q][r][s] = direct - exchange; // antisymmetrized <pq||rs> 
				}
	return TEI_SO;
}

Matrix build_f(const Matrix& h, const Tensor4& TEI_SO, int nso, int nelec) {
	Matrix f = Matrix::Zero(nso,nso); 
	for(int p=0; p < nso; p++)
		for(int q=0; q < nso; q++) {
			f(p,q) = h(spatial(p),spatial(q)) * same_spin(p,q);
			for(int m=0; m < nelec; m++)
				f(p,q) += TEI_SO[p][m][q][m];
		}
	return f;
}

Denoms build_denominators(const Matrix& f, int nso, int nelec) {
	Denoms D;
	D.D1 = Matrix::Zero(nso,nso);
	D.D2 = Tensor4(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));

	for(int i=0; i < nelec; i++) 
		for(int a=nelec; a < nso; a++) {
			D.D1(i,a) = f(i,i) - f(a,a);
		}

	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++)
			for(int a=nelec; a < nso; a++)
				for(int b=nelec; b < nso; b++) {
					D.D2[i][j][a][b] = f(i,i) + f(j,j) - f(a,a) - f(b,b);
				}
	return D;
}

Amplitudes initial_guess_amp(const Tensor4& TEI_SO, const Tensor4& D2, int nso, int nelec) {
	// build initial-guess cluster amplitudes 
	Amplitudes A;
	A.t1 = Matrix::Zero(nso,nso); // first order guess for singles -> 0 
	A.t2 = Tensor4(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));

	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++) 
			for(int a=nelec; a < nso; a++)
				for(int b=nelec; b < nso; b++) {
					A.t2[i][j][a][b] = TEI_SO[i][j][a][b] / D2[i][j][a][b];
				}
	return A;
}

double Emp2_guess(const Tensor4& TEI_SO, const Amplitudes& A, int nso, int nelec) {
	double Emp2 = 0.0; 
	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++) 
			for(int a=nelec; a < nso; a++)
				for(int b=nelec; b < nso; b++) {
					Emp2 += (0.25) * TEI_SO[i][j][a][b] * A.t2[i][j][a][b]; 
				}
	//printf("\nCluster amplitude-guess MP2 Energy: %.12f\n", Emp2);
	return Emp2;
}

double cc_energy(const Matrix& f, const Tensor4& TEI_SO, const Amplitudes& A, int nso, int nelec) {
	double Ecc = 0.0; 
	for(int i=0; i < nelec; i++)
		for(int a=nelec; a < nso; a++) {
			Ecc += f(i,a) * A.t1(i,a);
			for(int j=0; j < nelec; j++) 
				for(int b=nelec; b < nso; b++) {
					Ecc += (0.25) * TEI_SO[i][j][a][b] * A.t2[i][j][a][b];
					Ecc += (0.50) * TEI_SO[i][j][a][b] * A.t1(i,a) * A.t1(j,b);
				}
			}
	//printf("\n Coupled Cluster Energy: %.12f\n", Ecc);
	return Ecc;
}

Tensor4 build_tau(const Amplitudes& A, int nso, int nelec) {
	Tensor4 tau(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));
	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++) 
			for(int a=nelec; a < nso; a++)
				for(int b=nelec; b < nso; b++) {
					tau[i][j][a][b] = A.t2[i][j][a][b] + A.t1(i,a) * A.t1(j,b) -  A.t1(i,b) * A.t1(j,a);
				}
	return tau;
}

Tensor4 build_tau_tilde(const Amplitudes& A, int nso, int nelec) {
	Tensor4 tau_tilde(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));
	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++) 
			for(int a=nelec; a < nso; a++)
				for(int b=nelec; b < nso; b++) {
					tau_tilde[i][j][a][b] = A.t2[i][j][a][b] + 0.50 * (A.t1(i,a) * A.t1(j,b) -  A.t1(i,b) * A.t1(j,a));
				}
	return tau_tilde;
}

// CC Intermediates
Matrix build_Fae(const Matrix& fock, const Amplitudes& A, const Tensor4& TEI_SO, const Tensor4& tau_tilde, int nso, int nelec) {
	Matrix Fae = Matrix::Zero(nso,nso);
	for(int a=nelec; a < nso; a++)
		for(int e=nelec; e < nso; e++) {
			Fae(a,e) = (1 - (a == e)) * fock(a,e);
			for(int m=0; m < nelec; m++) {
				Fae(a,e) += -0.5 * fock(m,e) * A.t1(m,a);
				for(int f=nelec; f < nso; f++) {
					Fae(a,e) += A.t1(m,f) * TEI_SO[m][a][f][e];
					for(int n=0; n < nelec; n++)
						Fae(a,e) += -0.5 * tau_tilde[m][n][a][f] * TEI_SO[m][n][e][f];
				}
			}
		}
	return Fae;
}

Matrix build_Fmi(const Matrix& fock, const Amplitudes& A, const Tensor4& TEI_SO, const Tensor4& tau_tilde, int nso, int nelec) {
	Matrix Fmi = Matrix::Zero(nso,nso);
	for(int m=0; m < nelec; m++)
		for(int i=0; i < nelec; i++) {
			Fmi(m,i) = (1 - (m == i)) * fock(m,i);
			for(int e=nelec; e < nso; e++) {
				Fmi(m,i) += 0.5 * A.t1(i,e) * fock(m,e);
				for(int n=0; n < nelec; n++) {
					Fmi(m,i) += A.t1(n,e) * TEI_SO[m][n][i][e];
					for(int f=nelec; f < nso; f++)
						Fmi(m,i) += 0.5 * tau_tilde[i][n][e][f] * TEI_SO[m][n][e][f];
				}
			}
		}
	return Fmi;
}

Matrix build_Fme(const Matrix& fock, const Amplitudes& A, const Tensor4& TEI_SO, int nso, int nelec) {
	Matrix Fme = Matrix::Zero(nso,nso);
	for(int m=0; m < nelec; m++)
		for(int e=nelec; e < nso; e++) {
			Fme(m,e) = fock(m,e);
			for(int n=0; n < nelec; n++)
				for(int f=nelec; f < nso; f++)
					Fme(m,e) += A.t1(n,f) * TEI_SO[m][n][e][f]; 
		}
	return Fme;
}

//- Occupied: i, j, k, l, m, n
//- Virtual: a, b, c, d, e, f
Tensor4 build_Wmnij(const Amplitudes& A, const Tensor4& TEI_SO, const Tensor4& tau, int nso, int nelec) {
	Tensor4 Wmnij(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));
	for(int m=0; m < nelec; m++)
		for(int n=0; n < nelec; n++)
			for(int i=0; i < nelec; i++) 
				for(int j=0; j < nelec; j++) {
					Wmnij[m][n][i][j] = TEI_SO[m][n][i][j]; 
					for(int e=nelec; e < nso; e++) {
						Wmnij[m][n][i][j] += A.t1(j,e) * TEI_SO[m][n][i][e] - A.t1(i,e) * TEI_SO[m][n][j][e];
						for(int f=nelec; f < nso; f++) 
							Wmnij[m][n][i][j] += 0.25 * tau[i][j][e][f] * TEI_SO[m][n][e][f];
					}
				}
	return Wmnij;
}

Tensor4 build_Wabef(const Amplitudes& A, const Tensor4& TEI_SO, const Tensor4& tau, int nso, int nelec) {
	Tensor4 Wabef(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));
	for(int a=nelec; a < nso; a++)
		for(int b=nelec; b < nso; b++)
			for(int e=nelec; e < nso; e++) 
				for(int f=nelec; f < nso; f++) {
					Wabef[a][b][e][f] = TEI_SO[a][b][e][f]; 
					for(int m=0; m < nelec; m++) {
						Wabef[a][b][e][f] += -A.t1(m,b) * TEI_SO[a][m][e][f] + A.t1(m,a) * TEI_SO[b][m][e][f];
						for(int n=0; n < nelec; n++) 
							Wabef[a][b][e][f] += 0.25 * tau[m][n][a][b] * TEI_SO[m][n][e][f];
					}
				}
	return Wabef;
}

Tensor4 build_Wmbej(const Amplitudes& A, const Tensor4& TEI_SO, int nso, int nelec) { 
	Tensor4 Wmbej(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));
	for(int m=0; m < nelec; m++)
		for(int b=nelec; b < nso; b++)
			for(int e=nelec; e < nso; e++)
				for(int j=0; j < nelec; j++) {
					Wmbej[m][b][e][j] = TEI_SO[m][b][e][j];
					for(int f=nelec; f < nso; f++)
						Wmbej[m][b][e][j] += A.t1(j,f) * TEI_SO[m][b][e][f];
					for(int n=0; n < nelec; n++) {
						Wmbej[m][b][e][j] += -A.t1(n,b) * TEI_SO[m][n][e][j];
						for(int f=nelec; f < nso; f++) 
							Wmbej[m][b][e][j] += -(0.5 * A.t2[j][n][f][b] + A.t1(j,f) * A.t1(n,b)) * TEI_SO[m][n][e][f];
					}
				}
	return Wmbej;
}

// update cluster amplitudes
Amplitudes update_amp(const Matrix& fock, Amplitudes& init, const Denoms& D, const Tensor4& TEI_SO, const Matrix& Fae, const Matrix& Fmi, const Matrix& Fme, const Tensor4& tau, const Tensor4& Wmnij, const Tensor4& Wabef, const Tensor4& Wmbej, int nso, int nelec) { 
	Amplitudes A;
	A.t1 = Matrix::Zero(nso,nso);
	A.t2 = Tensor4(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));

	// T1 amplitudes
	for(int i=0; i < nelec; i++)
		for(int a=nelec; a < nso; a++) {
			A.t1(i,a) = fock(i,a);
			for(int e=nelec; e < nso; e++)
				A.t1(i,a) += init.t1(i,e) * Fae(a,e);
			for(int m=0; m < nelec; m++) {
				A.t1(i,a) += -init.t1(m,a) * Fmi(m,i);
				for(int e=nelec; e < nso; e++) {
					A.t1(i,a) += init.t2[i][m][a][e] * Fme(m,e);
					for(int n=0; n < nelec; n++)
						A.t1(i,a) += -0.5 * init.t2[m][n][a][e] * TEI_SO[n][m][e][i];
					for(int f=nelec; f < nso; f++)
						A.t1(i,a) += -0.5 * init.t2[i][m][e][f] * TEI_SO[m][a][e][f];
				}
			}
			for(int n=0; n < nelec; n++)
				for(int f=nelec; f < nso; f++)
					A.t1(i,a) += -init.t1(n,f) * TEI_SO[n][a][i][f];
			A.t1(i,a) = A.t1(i,a) / D.D1(i,a); // normalize
		}

	// T2 amplitudes
	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++) 
			for(int a=nelec; a < nso; a++)
				for(int b=nelec; b < nso; b++) {
					A.t2[i][j][a][b] += TEI_SO[i][j][a][b];
					for(int e=nelec; e < nso; e++) {
						double bracket = Fae(b,e);
						for(int m=0; m < nelec; m++)
							bracket -= 0.5 * init.t1(m,b) * Fme(m,e);
						A.t2[i][j][a][b] += init.t2[i][j][a][e] * bracket;
					}
					for(int e=nelec; e < nso; e++) {
						double bracket = Fae(a,e);
						for(int m=0; m < nelec; m++)
							bracket -= 0.5 * init.t1(m,a) * Fme(m,e);
						A.t2[i][j][a][b] += -init.t2[i][j][b][e] * bracket;
					}
					for(int m=0; m < nelec; m++) {
						double bracket = Fmi(m,j);
						for(int e=nelec; e < nso; e++)
							bracket += 0.5 * init.t1(j,e) * Fme(m,e);
						A.t2[i][j][a][b] -= init.t2[i][m][a][b] * bracket;
					}
					for(int m=0; m < nelec; m++) {
						double bracket = Fmi(m,i);
						for(int e=nelec; e < nso; e++)
							bracket += 0.5 * init.t1(i,e) * Fme(m,e);
						A.t2[i][j][a][b] -= -init.t2[j][m][a][b] * bracket;
					}
					for(int m=0; m < nelec; m++)
						for(int n=0; n < nelec; n++)
							A.t2[i][j][a][b] += 0.5 * tau[m][n][a][b] * Wmnij[m][n][i][j];
					for(int e=nelec; e < nso; e++)
						for(int f=nelec; f < nso; f++)
							A.t2[i][j][a][b] += 0.5 * tau[i][j][e][f] * Wabef[a][b][e][f];
					for(int m=0; m < nelec; m++)
						for(int e=nelec; e < nso; e++) 
							A.t2[i][j][a][b] += init.t2[i][m][a][e] * Wmbej[m][b][e][j] - init.t1(i,e) * init.t1(m,a) * TEI_SO[m][b][e][j];
					for(int m=0; m < nelec; m++)
						for(int e=nelec; e < nso; e++) 
							A.t2[i][j][a][b] -= init.t2[j][m][a][e] * Wmbej[m][b][e][i] - init.t1(j,e) * init.t1(m,a) * TEI_SO[m][b][e][i];
					for(int m=0; m < nelec; m++)
						for(int e=nelec; e < nso; e++) 
							A.t2[i][j][a][b] -= init.t2[i][m][b][e] * Wmbej[m][a][e][j] - init.t1(i,e) * init.t1(m,b) * TEI_SO[m][a][e][j];
					for(int m=0; m < nelec; m++)
						for(int e=nelec; e < nso; e++) 
							A.t2[i][j][a][b] += init.t2[j][m][b][e] * Wmbej[m][a][e][i] - init.t1(j,e) * init.t1(m,b) * TEI_SO[m][a][e][i];
					for(int e=nelec; e < nso; e++)
						A.t2[i][j][a][b] += init.t1(i,e) * TEI_SO[a][b][e][j];
					for(int e=nelec; e < nso; e++)
						A.t2[i][j][a][b] += -init.t1(j,e) * TEI_SO[a][b][e][i];
					for(int m=0; m < nelec; m++)
						A.t2[i][j][a][b] -= init.t1(m,a) * TEI_SO[m][b][i][j];
					for(int m=0; m < nelec; m++)
						A.t2[i][j][a][b] -= -init.t1(m,b) * TEI_SO[m][a][i][j];
					A.t2[i][j][a][b]  = A.t2[i][j][a][b]  / D.D2[i][j][a][b]; // normalize
				}
	return A;
}


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
