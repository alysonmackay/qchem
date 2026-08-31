// perturbative triples correction to CCSD --> CCSD(T) 
#include "triples.h"

using namespace std;

Tensor6 build_D3(const Matrix& f, int nso, int nelec) {
	Tensor6 D3 = Tensor6(nso, Vec5(nso, Vec4(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))))));

	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++)
			for(int k=0; k < nelec; k++)
				for(int a=nelec; a < nso; a++)
					for(int b=nelec; b < nso; b++) 
						for(int c=nelec; c < nso; c++) {
							D3[i][j][k][a][b][c] = f(i,i) + f(j,j) + f(k,k) - f(a,a) - f(b,b) - f(c,c);
						}
	return D3;
}

Tensor6 disconnected_triples(const Amplitudes& A, const Tensor4& TEI_SO, const Tensor6& D3, int nso, int nelec) {
	Tensor6 T3d = Tensor6(nso, Vec5(nso, Vec4(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))))));

	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++)
			for(int k=0; k < nelec; k++)
				for(int a=nelec; a < nso; a++)
					for(int b=nelec; b < nso; b++) 
						for(int c=nelec; c < nso; c++) {
							T3d[i][j][k][a][b][c] += A.t1(i,a) * TEI_SO[j][k][b][c];
							T3d[i][j][k][a][b][c] -= A.t1(j,a) * TEI_SO[i][k][b][c];
							T3d[i][j][k][a][b][c] -= A.t1(k,a) * TEI_SO[j][i][b][c];

							T3d[i][j][k][a][b][c] -= A.t1(i,b) * TEI_SO[j][k][a][c];
							T3d[i][j][k][a][b][c] += A.t1(j,b) * TEI_SO[i][k][a][c];
							T3d[i][j][k][a][b][c] += A.t1(k,b) * TEI_SO[j][i][a][c];

							T3d[i][j][k][a][b][c] -= A.t1(i,c) * TEI_SO[j][k][b][a];
							T3d[i][j][k][a][b][c] += A.t1(j,c) * TEI_SO[i][k][b][a];
							T3d[i][j][k][a][b][c] += A.t1(k,c) * TEI_SO[j][i][b][a];
							T3d[i][j][k][a][b][c] = T3d[i][j][k][a][b][c] / D3[i][j][k][a][b][c];
						}
	return T3d;
}

Tensor6 connected_triples(const Amplitudes& A, const Tensor4& TEI_SO, const Tensor6& D3, int nso, int nelec) {
	Tensor6 T3c = Tensor6(nso, Vec5(nso, Vec4(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))))));
	Tensor6 e_sum = Tensor6(nso, Vec5(nso, Vec4(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))))));
	Tensor6 m_sum = Tensor6(nso, Vec5(nso, Vec4(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))))));

	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++)
			for(int k=0; k < nelec; k++)
				for(int a=nelec; a < nso; a++)
					for(int b=nelec; b < nso; b++) 
						for(int c=nelec; c < nso; c++) {
							for(int e=nelec; e < nso; e++) {
								e_sum[i][j][k][a][b][c] += A.t2[j][k][a][e] * TEI_SO[e][i][b][c];
								e_sum[i][j][k][a][b][c] -= A.t2[i][k][a][e] * TEI_SO[e][j][b][c];
								e_sum[i][j][k][a][b][c] -= A.t2[j][i][a][e] * TEI_SO[e][k][b][c];

								e_sum[i][j][k][a][b][c] -= A.t2[j][k][b][e] * TEI_SO[e][i][a][c];
								e_sum[i][j][k][a][b][c] += A.t2[i][k][b][e] * TEI_SO[e][j][a][c];
								e_sum[i][j][k][a][b][c] += A.t2[j][i][b][e] * TEI_SO[e][k][a][c];

								e_sum[i][j][k][a][b][c] -= A.t2[j][k][c][e] * TEI_SO[e][i][b][a];
								e_sum[i][j][k][a][b][c] += A.t2[i][k][c][e] * TEI_SO[e][j][b][a];
								e_sum[i][j][k][a][b][c] += A.t2[j][i][c][e] * TEI_SO[e][k][b][a];
							}
							for(int m=0; m < nelec; m++) {
								m_sum[i][j][k][a][b][c] += A.t2[i][m][b][c] * TEI_SO[m][a][j][k];
								m_sum[i][j][k][a][b][c] -= A.t2[j][m][b][c] * TEI_SO[m][a][i][k];
								m_sum[i][j][k][a][b][c] -= A.t2[k][m][b][c] * TEI_SO[m][a][j][i];

								m_sum[i][j][k][a][b][c] -= A.t2[i][m][a][c] * TEI_SO[m][b][j][k];
								m_sum[i][j][k][a][b][c] += A.t2[j][m][a][c] * TEI_SO[m][b][i][k];
								m_sum[i][j][k][a][b][c] += A.t2[k][m][a][c] * TEI_SO[m][b][j][i];

								m_sum[i][j][k][a][b][c] -= A.t2[i][m][b][a] * TEI_SO[m][c][j][k];
								m_sum[i][j][k][a][b][c] += A.t2[j][m][b][a] * TEI_SO[m][c][i][k];
								m_sum[i][j][k][a][b][c] += A.t2[k][m][b][a] * TEI_SO[m][c][j][i];
							}
							T3c[i][j][k][a][b][c] = e_sum[i][j][k][a][b][c] - m_sum[i][j][k][a][b][c];
							T3c[i][j][k][a][b][c] = T3c[i][j][k][a][b][c] / D3[i][j][k][a][b][c];
						}
	return T3c;
}

double Etriples(const Tensor6& T3c, const Tensor6& T3d, const Tensor6& D3, int nso, int nelec) {
	double ET=0.0;

	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++)
			for(int k=0; k < nelec; k++)
				for(int a=nelec; a < nso; a++)
					for(int b=nelec; b < nso; b++) 
						for(int c=nelec; c < nso; c++) {
							double outer = T3c[i][j][k][a][b][c] * D3[i][j][k][a][b][c];
							double connected = outer * T3c[i][j][k][a][b][c];
							double disconnected = outer * T3d[i][j][k][a][b][c];
							ET += connected + disconnected;
							}
	ET = ET / 36;
	return ET;
}

