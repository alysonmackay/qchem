#include "mp2.h"

#include <iostream>
#include <chrono>

using namespace std;

double mp2_energy(const string& sys_in, const string& basis_in, const string& root_in) {

	// run HF SCF to get converged MOs
	Result conv_orbs = run_scf_verbose(sys_in, basis_in);
	double E_scf = conv_orbs.Escf;
	int nao = conv_orbs.C.rows();

	// AO to MO basis transformation 
	auto t0 = chrono::steady_clock::now(); 
	vector<double> TEI_MO = ao2mo(conv_orbs.TEI_AO, conv_orbs.C, nao); // N^5 algorithm
	auto t1 = chrono::steady_clock::now();
	chrono::duration<double> dt = t1 - t0;
	printf("\nTime for AO2MO Transformation (N^5): %.6f s\n", dt.count());
	
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

	return Emp2;
}
