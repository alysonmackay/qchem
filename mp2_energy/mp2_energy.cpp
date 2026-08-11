#include "hf.h"

#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;

static vector<int> ioff;
static int compound(int a, int b) {
	return a>b ? ioff[a] + b : ioff[b] + a;
}
// transform TEI_AO to MO basis 
vector<double> ao2mo_transform(const vector<double>& TEI_AO, const Matrix& coeff, int nao) {
	int size = TEI_AO.size();
	vector<double> TEI_MO(size);
	int i, j, k, l, ijkl; // MO indices 
	int p, q, r, s, pq, rs, pqrs; // AO indices
	for(i=0, ijkl=0; i < nao; i++) {
		for(j=0; j <= i; j++) {
			for(k=0; k <= i; k++) {
				for(l=0; l <= (i==k ? j : k); l++, ijkl++) { 

					for(p=0; p < nao; p++) {
						for(q=0; q < nao; q++) { 
							pq = compound(p,q);
							for(r=0; r < nao; r++) {
								for(s=0; s < nao; s++) {
									rs = compound(r,s); 
									pqrs = compound(pq, rs); 
									
									TEI_MO[ijkl] += coeff(p,i) * coeff(q,j) * coeff(r,k) * coeff(s,l) * TEI_AO[pqrs];
								}
							}
						}
					}
				}
			}
		}
        }
	return TEI_MO;
}



int main(int argc, char* argv[]) {

	// run HF SCF to get converged MOs
	Result conv_orbs = run_scf(argv[1], argv[2]);
	printf("\nConverged MO Coefficient Matrix:\n");
	print_matrix(conv_orbs.C);
	printf("\nConverged MO Energies:\n");
	cout << conv_orbs.epsi << endl;
	double E_scf = conv_orbs.Escf; 

	// read two-electron integrals in AO basis
	int nao = conv_orbs.C.rows();
	ifstream two_elec("../input/" + string(argv[1]) + "/" + argv[2] + "/eri.dat"); // kinda messy 
	int M = nao*(nao+1)/2; // number of distinict pairs
	ioff.resize(M);
	ioff[0] = 0; 
	for(int k=1; k < M; k++)
		ioff[k] = ioff[k-1] +k;
	int size = compound(M-1, M-1) + 1; 
	vector<double> TEI_AO(size);
	int mu, nu, la, si; // AO indices 
	double x;
	while(two_elec >> mu >> nu >> la >> si >> x) { 
		TEI_AO[ compound( compound(mu-1,nu-1), compound(la-1,si-1) ) ] = x; 
	}

	auto t0 = chrono::steady_clock::now(); 
	vector<double> TEI_MO = ao2mo_transform(TEI_AO, conv_orbs.C, nao);
	auto t1 = chrono::steady_clock::now();
	chrono::duration<double> dt = t1 - t0;
	printf("\nTime for AO2MO Transformation: %.12f s\n", dt.count());
	
	double Emp2 = 0.0;
	int ndocc = conv_orbs.nocc; // assume only closed shell for now
	int i, j, a, b, ia, ja, jb, ib, iajb, ibja; // i and j are doubly occupied; a and b unoccupied
	for(i=0; i < ndocc; i++) {
		for(a=ndocc; a < nao; a++) {
			ia = compound(i,a);
			for(j=0; j < ndocc; j++) { 
				ja = compound(j,a); 
				for(b=ndocc; b < nao; b++) {
					jb = compound(j,b); 
					ib = compound(i,b); 
					iajb = compound(ia,jb);
					ibja = compound(ib,ja);
					Emp2 += TEI_MO[iajb] * (2 * TEI_MO[iajb] - TEI_MO[ibja]) / (conv_orbs.epsi[i] + conv_orbs.epsi[j] - conv_orbs.epsi[a] - conv_orbs.epsi[b]);
					
				}
			}
		}
	}
	printf("\nSCF Energy: %.12f\n", E_scf); 
	printf("\nMP2 Energy: %.12f\n", Emp2);
	double E_tot = E_scf + Emp2; 
	printf("\nTotal Energy: %.12f\n", E_tot);
	
	return 0;
}
