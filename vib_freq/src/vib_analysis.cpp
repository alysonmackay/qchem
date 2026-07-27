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

	return 0;
}
