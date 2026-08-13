#include "mp2.h" // for spatial ao2mo transform; temporary 

#include <iostream>
#include <chrono>
#include <vector>

using namespace std;
using Vec1 = vector<double>; 
using Vec2 = vector<Vec1>; 
using Vec3 = vector<Vec2>;
using Tensor4 = vector<Vec3>;

int spatial(int k) {return k / 2; }
bool same_spin(int k, int l) {return k%2 == l%2; } 


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


int main(int argc, char* argv[]) 
{
	// run HF SCF to get converged MOs
	Result conv_orbs = run_scf(argv[1], argv[2]);
	printf("\nConverged MO Coefficient Matrix:\n");
	print_matrix(conv_orbs.C);
	printf("\nConverged MO Energies:\n");
	cout << conv_orbs.epsi << endl;
	double E_scf = conv_orbs.Escf;
	int nao = conv_orbs.C.rows();

	// AO to MO basis transformation (spatial) 
	auto t0 = chrono::steady_clock::now();
	vector<double> TEI_MO = ao2mo(conv_orbs.TEI_AO, conv_orbs.C, nao); // N^5 algorithm
	auto t1 = chrono::steady_clock::now();
	chrono::duration<double> dt = t1 - t0;
	printf("\nTime for AO2MO Transformation (N^5): %.6f s\n", dt.count());

	// spatial MO basis -> spin-orbital basis transform
	Tensor4 TEI_SO = mo2so(TEI_MO, nao);
	
	int nso = 2 * nao;
	int ndocc = conv_orbs.nocc;
	int nelec = 2 * ndocc; // number of occupied spin orbitals
	
	Matrix h = conv_orbs.C.transpose() * conv_orbs.H * conv_orbs.C; // Hamilontian in MO basis
	Matrix f = Matrix::Zero(nso,nso); 

	for(int p=0; p < nso; p++)
		for(int q=0; q < nso; q++) {
			f(p,q) = h(spatial(p),spatial(q)) * same_spin(p,q);
			for(int m=0; m < nelec; m++)
				f(p,q) += TEI_SO[p][m][q][m];
		}
	
	print_matrix(f); 

	return 0;
}
