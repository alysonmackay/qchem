#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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

int compound(int a, int b) {
	return a>b ? a*(a+1)/2 + b : b*(b+1)/2 + a; // TODO: refactor to use pre-computed lookup arrays
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

	// read nuclear repulsion energy
	ifstream enuc_data(npath);
	double enuc;
	enuc_data >> enuc;
	printf("Nuclear repulsion energy: %10.12f\n", enuc);

	// one-electron integrals
	ifstream overlap(spath); 	// read overlap in AO basis
	int i, j, n = 0;
	double x; 
	while(overlap >> i >> j >> x) {
		if(i > n) n = i; // n = number of basis functions 
	}
	overlap.close();

	Matrix S(n,n); 		// overlap
	fill(S, spath);
	Matrix T(n,n); 		// kinetic energy
	fill(T, tpath);
	Matrix V(n,n);		// nuclear-attraction
	fill(V, vpath);

	// core Hamiltonian
	Matrix H = T + V; 
	cout << H << endl;
	
	// two-electron integrals
	ifstream two_elec(two_path);
	int M = n*(n+1)/2; // number of distinict pairs
	int size = compound(M-1, M-1) + 1; 

	vector<double> eri(size);
	int k, l = 0;
	while(two_elec >> i >> j >> k >> l >> x) { 
		eri[ compound( compound(i-1,j-1), compound(k-1,l-1) ) ] = x; 
	}

	return 0;
}
