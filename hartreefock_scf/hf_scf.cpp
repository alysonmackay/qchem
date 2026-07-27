#include <iostream>
#include <fstream>
#include <string>

#include "Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

int main(int argc, char* argv[]) 
{
	string sys = argv[1];
	string basis = argv[2];

	string dir = "input/" + sys + "/" + basis + "/"; 
	string spath = dir + "s.dat";
	string tpath = dir + "t.dat";
	string vpath = dir + "v.dat"; 
	string npath = dir + "enuc.dat";

	// read nuclear repulsion energy
	ifstream enuc_data(npath);
	double enuc;
	enuc_data >> enuc;
	printf("Nuclear repulsion energy: %10.12f\n", enuc);

	// read overlap in AO basis
	ifstream overlap(spath);
	int r, c, n = 0;
	double x; 
	while(overlap >> r >> c >> x) {
		if(r > n) n = r; // n = number of basis functions 
	}
	Matrix S(n,n);
	overlap.clear();
	overlap.seekg(0);
	while(overlap >> r >> c >> x) {
		S(r-1, c-1) = x;
		S(c-1, r-1) = x;
	}

	// read kinetic energy
	Matrix T(n,n);
	ifstream kinetic(tpath);
	while(kinetic >> r >> c >> x) {
		T(r-1, c-1) = x;
		T(c-1, r-1) = x;
	}

	// read nuclear-attraction integrals
	Matrix V(n,n);
	ifstream nuclear_att(vpath);
	while(nuclear_att >> r >> c >> x) {
		V(r-1, c-1) = x;
		V(c-1, r-1) = x; // TODO: create helper function for filling matrices 
	}

	// core Hamiltonian
	Matrix H = T + V; 
	cout << H << endl;

	return 0;
}
