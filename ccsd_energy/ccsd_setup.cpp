#include "ccsd_setup.h"

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
