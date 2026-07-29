#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ratio>

#include "Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

void fill(Matrix& mat, const string& path) {
	ifstream in(path);
	int i, j = 0;
	double x;
	while(in >> i >> j >> x) {
		mat(i-1, j-1) = x; 
		mat(j-1, i-1) = x;
	}
}

vector<int> ioff;

int compound(int a, int b) {
	return a>b ? ioff[a] + b : ioff[b] + a;
}

int main(int argc, char* argv[]) 
{
	string sys = argv[1];
	string basis = argv[2];

	string dir = "input/" + sys + "/" + basis + "/"; 
	string spath = dir + "s.dat";
	string tpath = dir + "t.dat";
	string vpath = dir + "v.dat"; 
	string npath = dir + "enuc.dat";
	string two_path = dir + "eri.dat";
	string geom_path = dir + "geom.dat";

	// read nuclear repulsion energy
	ifstream enuc_data(npath);
	double enuc;
	enuc_data >> enuc;
	printf("Nuclear repulsion energy: %10.12f\n", enuc);

	// one-electron integrals
	ifstream overlap(spath); 	// read overlap in AO basis
	int i, j, nao = 0;
	double x; 
	while(overlap >> i >> j >> x) {
		if(i > nao) nao = i; // nao = number of basis functions 
	}
	overlap.close();

	Matrix S(nao,nao); 		// overlap
	fill(S, spath);
	Matrix T(nao,nao); 		// kinetic energy
	fill(T, tpath);
	Matrix V(nao,nao);		// nuclear-attraction
	fill(V, vpath);

	// core Hamiltonian
	Matrix H = T + V; 
	cout << H << endl;
	
	// two-electron integrals
	ifstream two_elec(two_path);
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
	//cout << S_half << endl; // TODO: fix print formating 
	
	// build inital guess density matrix
	Matrix F0 = S_half.transpose() * H * S_half;
	cout << F0 << endl; 
	// diagonalize Fock matrix 
	Eigen::SelfAdjointEigenSolver<Matrix> diag(F0);
	Matrix epsi0 = diag.eigenvalues(); // epsilon0; orbital energies 
	Matrix C0p = diag.eigenvectors(); // C0' coefficent matrix
	Matrix C0 = S_half * C0p; // transform evecs to original AO basis
	
	ifstream geom(geom_path);
	int natom;
	double zvals, coords;
	int total = 0;
	geom >> natom;
	for(int i=0; i < natom; i++) {
		geom >> zvals >> coords >> coords >> coords;
		total += zvals;
	}
	int nocc = total / 2; 
	cout << nocc << endl;

	Matrix D = Matrix::Zero(nao,nao);
	for(int mu=0; mu < nao; ++mu)
		for(int nu=0; nu < nao; ++nu)
			for(int m=0; m < nocc; ++m)
				D(mu, nu) += C0(mu,m) * C0(nu,m);
	cout << D << endl; // TODO: fix print formatting
	
	// compute initial SCF energy
	Matrix F = H;
	double elec0 = (D.array() * (H.array() + F.array())).sum();
	double E_tot = elec0 + enuc;
	cout << elec0 << endl; // TODO: fix print formatting
	
	// compute new Fock matrix using last density
	for(int i=0; i < nao; i++)
		for(int j=0; j < nao; j++) { 
			F(i,j) = H(i,j); 
			for(int k=0; k < nao; k++)
				for(int l=0; l < nao; l++) {
     					int ij = compound(i,j);
     					int kl = compound(k,l); 
     					int ijkl = compound(ij,kl); 
     					int ik = compound(i,k); 
     					int jl = compound(j,l); 
     					int ikjl = compound(ik,jl); 
     					F(i,j) += D(k,l) * (2.0 * eri[ijkl] - eri[ikjl]);
     				}
		}
	cout << F << endl;
				
	return 0;
}
