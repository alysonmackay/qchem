#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) 
{
	// read nuclear repulsion energy
	ifstream enuc_data(argv[1]);
	double enuc;
	enuc_data >> enuc;

	printf("Nuclear repulsion energy: %10.12f\n", enuc);

	return 0;
}
