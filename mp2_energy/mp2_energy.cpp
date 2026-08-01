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
	int i, j, k, l = 0;
	double x;
	while(two_elec >> i >> j >> k >> l >> x) { 
		TEI_AO[ compound( compound(i-1,j-1), compound(k-1,l-1) ) ] = x; 
	}

	// transform TEI_AO to MO basis 
	vector<double> TEI_MO(size);
	// first attempt: Noddy algorithm (N^8)
	
	return 0;
}
