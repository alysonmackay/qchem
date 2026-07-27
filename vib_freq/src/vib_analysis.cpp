#include "molecule.h"
#include "masses.h"

#include <iostream> 
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])
{ // argument handling
	if (argc < 4) {
		cerr << "Usage: " << argv[0] << " <filename> [charge] <hessian>\n"; 
		return 1; 
	} 
	const char* filename = argv[1];
	int charge = 0; // default if not supplied 
	
	if (argc >= 4) {
		try { 
			size_t pos; 
			charge = stoi(argv[2], &pos); // convert charge arg to int
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

	// reading hessian matrix
	ifstream hessian(argv[3]);
	int hess_natom; 
	hessian >> hess_natom;
	// confirm that hessian matches geom 
	if(mol.natom != hess_natom) {
		printf("Error. Number of atoms in hessian does not match geometry coordinates");
		return 1;
	}

	// converting from 3x(3N^2) to 3N x 3N 
	double **H = new double* [mol.natom*3];
	for(int i=0; i < mol.natom*3; i++)
		H[i] = new double[mol.natom*3]; 

	for(int i=0; i < mol.natom*3; i++) {
     		for(int j=0; j < mol.natom; j++) {
			hessian >> H[i][3*j] >> H[i][3*j+1] >> H[i][3*j+2];
     		}
	}

	// replace matrix with mass-weighted version 
	for(int i=0; i < mol.natom*3; i++) {
		for(int j=0; j < mol.natom*3; j++) {
			double mass_i = atomic_mass[mol.zvals[i/3]];
			double mass_j = atomic_mass[mol.zvals[j/3]];
			H[i][j] = H[i][j] / sqrt(mass_i * mass_j);
     		}
	}
	// print matrix
	for(int i=0; i < mol.natom*3; i++) {
		for(int j=0; j < mol.natom*3; j++)
     			printf("%13.7f", H[i][j]);
		printf("\n"); 
	}

	return 0;
}
