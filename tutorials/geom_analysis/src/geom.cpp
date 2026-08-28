#include "molecule.h"
#include "masses.h"

#include <iostream> 
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath> 

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

int main(int argc, char* argv[])
{	// argument handling
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " <filename> [charge]\n"; 
		return 1; 
	} 
	const char* filename = argv[1];
	int charge = 0; // default if not supplied 
	
	if (argc >= 3) {
		try { 
			size_t pos; 
			charge = stoi(argv[2], &pos); 
			if (pos != string(argv[2]).size()) {
				cerr << "Warning. Trailing characters ignored\n"; 
			}
		} catch (const invalid_argument&) {
			cerr << "Error: charge argument '" << argv[2] << "' is not a valid integer\n"; 
			return 1;
		} catch (const out_of_range&) {
			cerr << "Error: charge argument '" << argv[2] << "' is out of range for an int\n"; 
			return 1;
		}
	}

	// main program 
	Molecule mol(filename, charge);

	cout << "Number of atoms: " << mol.natom << endl; 
	cout << "Input Cartesian coordinates:\n"; 
	mol.print_geom(); 

	cout << "Interatomic distances (bohr):\n"; 
	for(int i=0; i < mol.natom; i++)
     		for(int j=0; j < i; j++)
			printf("%d %d %8.5f\n", i, j, mol.bond(i,j));

	cout << "\nBond angles:\n";
	for(int i=0; i < mol.natom; i++) {
		for(int j=0; j < i; j++) {
     			for(int k=0; k < j; k++) { 
     				if(mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0)
     					printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i,j,k)*(180.0/acos(-1.0)));
     			}
     		}
	}

	cout << "\nOut-of-plane angles:\n"; 
	for(int i=0; i < mol.natom; i++) {
		for(int k=0; k < mol.natom; k++) {
     			for(int j=0; j < mol.natom; j++) {
     				for(int l=0; l < j; l++) {
     					if(i!=j && i!=k && i!=l && j!=k && k!=l && mol.bond(i,k) < 4.0 && mol.bond(k,j) < 4.0 && mol.bond(k,l) < 4.0)
     						printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.oop(i,j,k,l)*(180.0/acos(-1.0)));
     				}
     			}
     		}
	}

	cout << "\nTorsion angles:\n";
	for(int i=0; i < mol.natom; i++) {
		for(int j=0; j < i; j++) {
			for(int k=0; k < j; k++) {
     				for(int l=0; l < k; l++) {
     					if(mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0 && mol.bond(k,l) < 4.0)
     						printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.torsion(i,j,k,l)*(180.0/acos(-1.0)));
     				}
			}
		}
	}
	// finding center of mass 
	double M = 0.0;
	for(int i=0; i < mol.natom; i++) M += atomic_mass[(int) mol.zvals[i]];

	double x_cm=0.0;
	double y_cm=0.0;
	double z_cm=0.0;
	double mi; 
	for(int i=0; i < mol.natom; i++) {
		mi = atomic_mass[(int) mol.zvals[i]];
		x_cm += mi * mol.geom[i][0];
		y_cm += mi * mol.geom[i][1];
		z_cm += mi * mol.geom[i][2];
	}
	x_cm /= M;
	y_cm /= M;
	z_cm /= M;
	printf("\nMolecular center of mass: %12.8f %12.8f %12.8f\n", x_cm, y_cm, z_cm);

	mol.translate(-x_cm, -y_cm, -z_cm);

	// moment of inertia tensor
	Matrix I(3,3); // invoke 3x3 matrix
	
	// fill tensor (diagonals first)
	for(int i=0; i < mol.natom; i++) {
		mi = atomic_mass[(int) mol.zvals[i]];
		I(0,0) += mi * (mol.geom[i][1]*mol.geom[i][1] + mol.geom[i][2]*mol.geom[i][2]);
		I(1,1) += mi * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][2]*mol.geom[i][2]);
		I(2,2) += mi * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][1]*mol.geom[i][1]);

		I(0,1) -= mi * mol.geom[i][0]*mol.geom[i][1];
		I(0,2) -= mi * mol.geom[i][0]*mol.geom[i][2];
		I(1,2) -= mi * mol.geom[i][1]*mol.geom[i][2];
	}
	// symmetry 
	I(1,0) = I(0,1);
	I(2,0) = I(0,2);
	I(2,1) = I(1,2);

	cout << fixed << setprecision(7);
	cout << "\nMomemnt of inertia tensor (amu bohr^2):\n";
	cout << I << endl;
	
	// principal moment of inertia (diagonalization of I)
	Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
	Matrix evecs = solver.eigenvectors();
	Matrix evals = solver.eigenvalues();

	cout << "\nPrincipal moments of inertia (amu * bohr^2):\n";
	cout << evals << endl;

	double conv = 0.529177249 * 0.529177249;
	cout << "\nPrincipal moments of inertia (amu * AA^2):\n";
	cout << evals * conv << endl;

	conv = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
	cout << "\nPrincipal moments of inertia (g * cm^2):\n";
	cout << evals * conv << endl;  // TODO: fix precision problem here 
	
	// assigning rotor types 
	if(mol.natom == 2) cout << "\nMolecule is diatomic\n";
	else if(evals(0) < 1e-4) cout << "\nMolecule is linear\n";
	else if(fabs(evals(0) - evals(1)) < 1e-4 && (fabs(evals(1) - evals(2)) < 1e-4))
		cout << "\nMolecule is a spherical top\n";
	else if(fabs(evals(0) - evals(1)) < 1e-4 && (fabs(evals(1) - evals(2)) > 1e-4))
		cout << "\nMolecule is an oblate symmetric top\n";
	else if(fabs(evals(0) - evals(1)) > 1e-4 && (fabs(evals(1) - evals(2)) < 1e-4))
		cout << "\nMolecule is a prolate symmetric top \n";
	else cout << "\nMolecule is an asymmetric top\n";

	// compute the rotational constants 
	double _pi = acos(-1.0);
	conv = 6.6260755E-34/(8.0 * _pi * _pi);
	conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
	conv *= 1e-6;
	cout << "\nRotational constants (MHz):\n";
	cout << "\tA = " << conv/evals(0) << "\t B = " << conv/evals(1) << "\t C = " << conv/evals(2) << endl;

	conv = 6.6260755E-34/(8.0 * _pi * _pi);
	conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
	conv /= 2.99792458E10;
	cout << "\nRotational constants (cm-1):\n";
	cout << "\tA = " << conv/evals(0) << "\t B = " << conv/evals(1) << "\t C = " << conv/evals(2) << endl;

	return 0;
}
