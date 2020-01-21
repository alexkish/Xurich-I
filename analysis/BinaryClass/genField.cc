#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>

#include "genField.hh"
//typedef std::vector<bool>::reference BoolRef;
// Overloading the hell out of the constructor for any forseable data types

genField::genField(double* x,string fieldName,string matrixDimens, ofstream* writingFile) 
{
	inputData.dbData = x;
	Name = fieldName;
	headerPart = fieldName + string(";double;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 1;
	isVect = false;
}

genField::genField(vector<double> * x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.vect_dbData = x;
	Name = fieldName;
	headerPart = fieldName + string(";vector<double>;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 2;
	maxSize = 0;
	isVect = true;
}

genField::genField(int* x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.intData = x;
	Name = fieldName;
	headerPart = fieldName + string(";int32;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 3;
	isVect = false;
}

genField::genField(short * x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.shortData = x;
	Name = fieldName;
	headerPart = fieldName + string(";int16;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 4;
	isVect = false;
}
genField::genField(long long * x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.longlongData = x;
	Name = fieldName;
	headerPart = fieldName + string(";int64;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 5;
	isVect = false;
}
genField::genField(unsigned int * x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.uintData = x;
	Name = fieldName;
	headerPart = fieldName + string(";uint32;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 6;
	isVect = false;
}
genField::genField(unsigned short * x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.ushortintData = x;
	Name = fieldName;
	headerPart = fieldName + string(";uint32;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 7;
	isVect = false;
}
genField::genField(unsigned long long * x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.ulonglongData = x;
	Name = fieldName;
	headerPart = fieldName + string(";uint64;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 8;
	isVect = false;
}
genField::genField(bool * x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.boolData = x;
	Name = fieldName;
	headerPart = fieldName + string(";logical;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 9;
	isVect = false;
}
/*
genField::genField(BoolRef x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.vect_boolData = x;
	Name = fieldName;
	headerPart = fieldName + string(";vector<logical>;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 10;
	maxSize = 0;
	isVect = true;
}

genField::genField(vector<bool> * x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.vect_boolData = x;
	Name = fieldName;
	headerPart = fieldName + string(";vector<logical>;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 10;
	maxSize = 0;
	isVect = true;
}*/
genField::genField(char * x, string fieldName,string matrixDimens, ofstream* writingFile)
{
	inputData.charData = x;
	Name = fieldName;
	headerPart = fieldName + string(";char;") + matrixDimens+string(";");
	outF = writingFile;
	dataType = 11;
	isVect = false;
}

genField::~genField()
{

}

void genField::Write()
{
	switch (dataType) 
	{
		case 1 :
			outF->write( (char *) inputData.dbData, sizeof *(inputData.dbData));
			break;
		case 2 :
			// VECTOR
			vectLength = (inputData.vect_dbData)->size();
			
		        if(vectLength==0)
			  {
			    (inputData.vect_dbData)->push_back(0.);
			    vectLength = (inputData.vect_dbData)->size();
			  }
			
			if (vectLength>maxSize)
			{
				maxSize = vectLength;
			}
			outF->write((char *) &vectLength, sizeof vectLength);
			outF->write((char *) &((*(inputData.vect_dbData))[0]),vectLength*sizeof(double));
			inputData.vect_dbData->clear();
			break;
		case 3 :
			outF->write( (char *) inputData.intData, sizeof *(inputData.intData));
			break;
		case 4 :
			outF->write( (char *) inputData.shortData, sizeof *(inputData.shortData));
			break;
		case 5 :
			outF->write( (char *) inputData.longlongData, sizeof *(inputData.longlongData));
			break;
		case 6 :
			outF->write( (char *) inputData.uintData, sizeof *(inputData.uintData));
			break;
		case 7 :
			outF->write( (char *) inputData.ushortintData, sizeof *(inputData.ushortintData));
			break;
		case 8 :
			outF->write( (char *) inputData.ulonglongData, sizeof *(inputData.ulonglongData));
			break;
		case 9 :
			outF->write( (char *) inputData.boolData, sizeof *(inputData.boolData));
			break;
		case 10 :  
/*			// VECTOR
			vectLength = (inputData.vect_boolData)->size();
			if (vectLength>maxSize)
			{
				maxSize = vectLength;
			}
			outF->write((char *) &vectLength, sizeof vectLength);
			outF->write((char *) &((*(inputData.vect_boolData))[0]),vectLength*sizeof(bool));
			inputData.vect_boolData->clear();  */
			break;			
		case 11 :
			outF->write( (char *) inputData.charData, sizeof *(inputData.charData));
			break; 
		default :
			break;
	}	

}

void genField::Clear()
{
	switch (dataType) 
	{
		case 1 :
			*(inputData.dbData) = 0;
			break;
		case 2 :
			inputData.vect_dbData->clear();
			break;
		case 3 :
			*(inputData.intData) = 0;
			break;
		case 4 :
			*(inputData.shortData) = 0;
			break;
		case 5 :
			*(inputData.longlongData) = 0;
			break;
		case 6 :
			*(inputData.uintData) = 0;
			break;
		case 7 :
			*(inputData.ushortintData) = 0;
			break;
		case 8 :
			*(inputData.ulonglongData) = 0;
			break;
		case 9 :
			*(inputData.boolData) = 0;
			break;
		case 10 :  
			//inputData.vect_boolData->clear();
			break;			
		case 11 :
			*(inputData.charData) = 0;
			break; 
		default :
			break;
	}	
	
}

string genField::GetHeaderPart()
{
	return headerPart;
}


// that should be it... but need to overload the constructor
