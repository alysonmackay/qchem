#include "cc_int.h"
// Equations adapted from Stanton 1991. DOI: 10.1063/1.460620

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
