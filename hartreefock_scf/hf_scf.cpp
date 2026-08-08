#include "hf.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

string sys, basis, root;
vector<int> ioff; 
int nao;

string path(const string& name) {
	string dir = root + "/" + sys + "/" + basis + "/";
	string filename = dir + name + ".dat"; 
	return filename;
}

void fill(Matrix& mat, const string& path) {
	ifstream in(path);
	int i, j = 0;
	double x;
	while(in >> i >> j >> x) {
		mat(i-1, j-1) = x; 
		mat(j-1, i-1) = x;
	}
}

void print_matrix(const Matrix& mat) {
	for(int start = 0; start < mat.rows(); start += 10) {
		int end = min<int>(start + 10, mat.rows());
		printf("\n");
		for(int c=start; c < end; c++) printf("%12d", c+1);
		printf("\n\n");
		for(int r=0; r < mat.rows(); r++) {
     			printf("%5d", r+1);
     			for(int c=start; c < end; c++) printf("%12.7f", mat(r,c));
     			printf("\n");
		}
	}
	printf("\n");
}

int compound(int a, int b) {
	return a>b ? ioff[a] + b : ioff[b] + a;
}

Matrix build_Fock(const Matrix& H, const Matrix& D, const vector<double>& eri) {
	Matrix F(H.rows(), H.cols());
	for(int i=0; i < nao; i++)
		for(int j=0; j < nao; j++) { 
			F(i,j) = H(i,j);
			int ij = compound(i,j);
			for(int k=0; k < nao; k++) {
				int ik = compound(i,k);
				for(int l=0; l < nao; l++) {
     					int kl = compound(k,l); 
     					int ijkl = compound(ij,kl); 
     					int jl = compound(j,l); 
     					int ikjl = compound(ik,jl); 
     					F(i,j) += D(k,l) * (2.0 * eri[ijkl] - eri[ikjl]);
     				}
			}
		}
	return F;
}

Matrix build_density(const Matrix& C, int nocc) {
	Matrix D = Matrix::Zero(nao,nao);
	for(int i=0; i < nao; i++)
		for(int j=0; j < nao; j++)
			for(int m=0; m < nocc; m++)
				D(i, j) += C(i,m) * C(j,m);
	return D;
}

Result run_scf(const string& sys_in, const string& basis_in, const string& root_in) {

	sys = sys_in; 
	basis = basis_in;
	root = root_in;


	printf("\nMolecule: %s\nBasis: %s\n", sys.c_str(), basis.c_str());

	// read nuclear repulsion energy
	ifstream enuc_data(path("enuc"));
	double enuc;
	enuc_data >> enuc;
	printf("\nNuclear repulsion energy = %20.12f\n", enuc);

	// one-electron integrals
	ifstream overlap(path("s")); 	// read overlap in AO basis
	int i, j = 0;
	double x; 
	while(overlap >> i >> j >> x) {
		if(i > nao) nao = i; // nao = number of basis functions 
	}
	overlap.close();

	Matrix S(nao,nao); 		// overlap
	fill(S, path("s"));
	printf("\nOverlap Integrals:\n");
	print_matrix(S);
	Matrix T(nao,nao); 		// kinetic energy
	fill(T, path("t"));
	printf("\nKinetic-Energy Integrals:\n");
	print_matrix(T);
	Matrix V(nao,nao);		// nuclear-attraction
	fill(V, path("v"));
	printf("\nNuclear Attraction Integrals:\n");
	print_matrix(V);

	// core Hamiltonian
	Matrix H = T + V; 
	printf("\nCore Hamiltonian:\n");
	print_matrix(H);
	
	// two-electron integrals
	ifstream two_elec(path("eri"));
	int M = nao*(nao+1)/2; // number of distinict pairs
	ioff.resize(M);
	ioff[0] = 0; 
	for(int k=1; k < M; k++)
		ioff[k] = ioff[k-1] +k;
	int size = compound(M-1, M-1) + 1; 
	vector<double> eri(size);
	int k, l = 0;
	while(two_elec >> i >> j >> k >> l >> x) { 
		eri[ compound( compound(i-1,j-1), compound(k-1,l-1) ) ] = x; 
	}

	// orthogonalization matrix
	Vector Lambda(nao);
	Eigen::SelfAdjointEigenSolver<Matrix> solver(S);
	Matrix evals = solver.eigenvalues();
	Matrix evecs = solver.eigenvectors();
	for(int h=0; h < nao; h++) {
		Lambda(h) = pow(evals(h),(-1.0/2.0));
	}
	Matrix S_half = evecs * Lambda.asDiagonal() * evecs.transpose();
	printf("\nS^-1/2 Matrix:\n");
	print_matrix(S_half);

	// build inital guess density matrix
	Matrix F = S_half.transpose() * H * S_half;
	printf("\nInitial Fock Matrix:\n");
	print_matrix(F);
	// diagonalize Fock matrix 
	Eigen::SelfAdjointEigenSolver<Matrix> diag(F);
	Matrix epsi = diag.eigenvalues(); // epsilon0; orbital energies 
	Matrix Cp = diag.eigenvectors(); // C' coefficent matrix
	Matrix C = S_half * Cp; // backtransform evecs to original AO basis
	printf("\nInitial C Matrix:\n");
	print_matrix(C);
	
	ifstream geom(path("geom"));
	int natom;
	double zvals, coords;
	int total = 0;
	geom >> natom;
	for(int i=0; i < natom; i++) {
		geom >> zvals >> coords >> coords >> coords;
		total += zvals;
	}
	int nocc = total / 2; 
	printf("\nNumber of occupied orbitals = %d\n", nocc);

	Matrix D = build_density(C, nocc);
	printf("\nInitial Density Matrix:\n");
	print_matrix(D);
	printf("%-10s%10s%20s%20s%20s\n", "Iter", "E(elec)", "E(tot)", "Delta(E)", "RMS(D)"); 
	double E_elec0 = (D.array() * (2 * H.array())).sum();
	double E_tot0 = E_elec0 + enuc; 
	printf("00 %20.12f %20.12f\n", E_elec0, E_tot0);

	// begin SCF iterations until energy convergence 
	double E_old = E_elec0;
	for(int iter=1; iter <= 100; iter++) { 
		F = build_Fock(H, D, eri);

		// new orbs 
		Matrix Fp = S_half.transpose() * F * S_half; 
		Eigen::SelfAdjointEigenSolver<Matrix> solv(Fp);
		epsi = solv.eigenvalues();
		C = S_half * solv.eigenvectors();
		Matrix D_old = D; 
		D = build_density(C, nocc); 

		double E_elec = (D.array() * (H.array() + F.array())).sum();
		double E_tot = E_elec + enuc;
		double dE = E_elec - E_old; 
		double rms = (D - D_old).norm(); 

		printf("%02d %20.12f %20.12f %20.12f %20.12f\n", iter, E_elec, E_tot, dE, rms); 
		if(fabs(dE) < 1e-12 && rms < 1e-11) {
			printf("\nThe SCF energy has converged!\n");
			break;
		} 
		E_old = E_elec;
	}

	// MO basis Fock matrix 
	Matrix F_MO = C.transpose() * F * C;
	printf("\nFock Matrix in MO Basis:\n");
	print_matrix(F_MO);

	Result MO; 
	MO.C = C; 
	MO.epsi = epsi;
	MO.nocc = nocc;

	return MO;

}
