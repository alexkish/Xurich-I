#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector>

#include "genStruct.hh"

using namespace std;

/*struct dataHolder 
{
	double *testDbl1;
	vector<double> * testDbl2;
	int *testInt;
};
*/
int main()
{
	genStruct * binFile = new genStruct("test_Struct.dat");

	double *A = new double;
	int *B = new int;
	bool *C = new bool;
	vector<double> *vectD = new vector<double>;

	binFile->AddField(A,"A","1,1");
	binFile->AddField(B,"B","1,1");
	binFile->AddField(C,"C","1,1");
	binFile->AddField(vectD,"vectD","1,1");

	binFile->MakeHeader();
	
	// first iteration
	*A = 3.1415;
	*B = 2;
	*C = true;
	vectD->clear();
	vectD->push_back(1.1);
	vectD->push_back(22.445);
	binFile->Write();

	// second iteration
	*A = 33.6;
	*B = 455;
	*C = false;
	vectD->clear();
	vectD->push_back(323.234234);
	vectD->push_back(57.2323);
	vectD->push_back(102.334);
	vectD->push_back(27.1);
	vectD->push_back(-44.2211);
	binFile->Write();


	delete binFile;
	delete A;
	delete B;
	delete C;
	delete vectD;
	return 0;
}
