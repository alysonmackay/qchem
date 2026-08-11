#include "hf.h"

#include <iostream>
#include <fstream>

using namespace std;

static vector<int> ioff;
static int compound(int a, int b) {
	return a>b ? ioff[a] + b : ioff[b] + a;
}

int main(int argc, char* argv[]) {

	// run HF SCF to get converged MOs
	Result conv_orbs = run_scf(argv[1], argv[2]);
	printf("\nConverged MO Coefficient Matrix:\n");
	print_matrix(conv_orbs.C);
	printf("\nConverged MO Energies:\n");
	cout << conv_orbs.epsi << endl;

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
	int mu, nu, la, si = 0; // AO indices 
	double x;
	while(two_elec >> mu >> nu >> la >> si >> x) { 
		TEI_AO[ compound( compound(mu-1,nu-1), compound(la-1,si-1) ) ] = x; 
	}

	// transform TEI_AO to MO basis 
	vector<double> TEI_MO(size);
	// first attempt: Noddy algorithm (N^8)
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
									
									TEI_MO[ijkl] += conv_orbs.C(p,i) * conv_orbs.C(q,j) * conv_orbs.C(r,k) * conv_orbs.C(s,l) * TEI_AO[pqrs];
								}
							}
						}
					}
				}
			}
		}
        }

	double Emp2 = 0.0;
	int ndocc = conv_orbs.nocc; // assume only closed 
        int a, b, ia, ja, jb, ib, iajb, ibja; // i and j are doubly occupied; a and b unoccupied
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
	printf("\nMP2 Energy:");
	cout << Emp2 << endl;
	
	return 0;
}
