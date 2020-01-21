#ifndef genStruct_h
#define genStruct_h 1

#include<iostream>
#include<stdio.h>
#include<vector>
#include<string.h>
#include<fstream>

using namespace std;
//typedef std::vector<bool>::reference BoolRef;
class genField;

class genStruct
{

	public:
		genStruct(string fileName);
		~genStruct();
	
	public:
		void Clear();
		void Write();
		void MakeHeader();
		void MakeFooter();
		void AddField(double* x, string fieldName,string matrixDimens);
		void AddField(vector<double> * x,string fieldName,string matrixDimens);
		void AddField(int* x, string fieldName,string matrixDimens);
		void AddField(short * x, string fieldName,string matrixDimens);
		void AddField(long long * x, string fieldName,string matrixDimens);
		void AddField(unsigned int * x, string fieldName,string matrixDimens);
		void AddField(unsigned short * x, string fieldName,string matrixDimens);
		void AddField(unsigned long long * x, string fieldName,string matrixDimens);
		void AddField(bool * x, string fieldName,string matrixDimens);
//		void AddField(vector<bool> * x, string fieldName,string matrixDimens);
//		void AddField(BoolRef x, string fieldName,string matrixDimens);
		void AddField(char * x, string fieldName,string matrixDimens);
	
	private:
		vector<genField> fieldVect;
		ofstream * theFile;
		unsigned long int numEvents;
		bool headerMade;
		bool footerMade;
	

};

#endif
