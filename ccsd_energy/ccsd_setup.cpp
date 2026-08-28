#include "ccsd_setup.h"

using namespace std;

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
	return Emp2;
}
