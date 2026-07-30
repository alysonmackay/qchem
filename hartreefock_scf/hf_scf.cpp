#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

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

int nao; 

Matrix build_Fock(const Matrix& H, const Matrix& D, const vector<double>& eri) {
	Matrix F(H.rows(), H.cols());
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
	int i, j = 0;
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
	Matrix F = S_half.transpose() * H * S_half;
	cout << F << endl; 
	// diagonalize Fock matrix 
	Eigen::SelfAdjointEigenSolver<Matrix> diag(F);
	Matrix epsi = diag.eigenvalues(); // epsilon0; orbital energies 
	Matrix Cp = diag.eigenvectors(); // C' coefficent matrix
	Matrix C = S_half * Cp; // backtransform evecs to original AO basis
	
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

	printf("Iter		E(elec) 	E(tot)	    		Delta(E)		RMS(D)\n"); 
	Matrix D = build_density(C, nocc);
	double E_elec0 = (D.array() * (H.array() + H.array())).sum();
	double E_tot0 = E_elec0 + enuc; 
	printf("00 %20.12f %20.12f\n", E_elec0, E_tot0);
	
	double E_old = 0.0;
	for(int iter=1; iter <= 100; iter++) { 
		F = build_Fock(H, D, eri);

		// new orbs 
		Matrix Fp = S_half.transpose() * F * S_half; 
		Eigen::SelfAdjointEigenSolver<Matrix> solv(Fp);
		Matrix C = S_half * solv.eigenvectors();
		Matrix D_old = D; 
		D = build_density(C, nocc); 

		double E_elec = (D.array() * (H.array() + F.array())).sum();
		double E_tot = E_elec + enuc;
		double dE = E_elec - E_old; 
		double rms = (D - D_old).norm(); 

		printf("%02d %20.12f %20.12f %20.12f %20.12f\n", iter, E_elec, E_tot, dE, rms); 
		if(fabs(dE) < 1e-12 && rms < 1e-11) break; 
		E_old = E_elec;
	}

	return 0;
}
