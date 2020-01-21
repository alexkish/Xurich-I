#include "genStruct.hh"
#include "genField.hh"
#include <sstream>

//typedef std::vector<bool>::reference BoolRef;

genStruct::genStruct(string fileName)
{
	theFile = new ofstream(fileName.c_str(), ios::out|ios::binary);
	Clear();
	numEvents = 0;
	headerMade = false;
	footerMade = false;
}

genStruct::~genStruct()
{
	if (!footerMade)
		MakeFooter();
	theFile->close();
	delete theFile;
}

void genStruct::Clear()
{
	vector<genField>::iterator itFieldVect;
	for(itFieldVect = fieldVect.begin(); itFieldVect!=fieldVect.end();itFieldVect++)
	{
		itFieldVect->Clear();
	}
}

void genStruct::Write()
{
	vector<genField>::iterator itFieldVect;
	for(itFieldVect = fieldVect.begin(); itFieldVect!=fieldVect.end();itFieldVect++)
	{
		itFieldVect->Write();
	}
	numEvents++;
}

void genStruct::MakeHeader()
{
	if (!headerMade)
	{
		string fileHeader = "";
		vector<genField>::iterator itFieldVect;
		for(itFieldVect = fieldVect.begin(); itFieldVect!=fieldVect.end();itFieldVect++)
		{
			//headerPart = (*itFieldVect).GetHeaderPart();
			//headerPart = itFieldVect->GetHeaderPart();
			fileHeader = fileHeader + itFieldVect->GetHeaderPart();
		}
		unsigned short fileHeaderSize = fileHeader.size();
		theFile->write( (char *) &fileHeaderSize, sizeof fileHeaderSize);
		theFile->write(fileHeader.c_str(), fileHeaderSize);
		headerMade = true;
	}
}


void genStruct::MakeFooter()
{
	footerMade = true;
	std::stringstream hStreamDist;	
	std::string fileFooter = "";
	vector<genField>::iterator itFieldVect;
	for(itFieldVect = fieldVect.begin(); itFieldVect!=fieldVect.end();itFieldVect++)
	{
		if (itFieldVect->isVect)
		{
			hStreamDist << itFieldVect->maxSize;
			fileFooter = fileFooter + itFieldVect->Name + string(";") + string(hStreamDist.str()) + 
						 string(";");
			hStreamDist.clear();
			hStreamDist.str("");
		}
	
	}
	
	unsigned short fileFooterSize = fileFooter.size();
	theFile->write(fileFooter.c_str(), fileFooterSize);
	theFile->write((char *) &numEvents, sizeof numEvents);
	theFile->write((char *) &fileFooterSize, sizeof fileFooterSize);
}

// ---------------------------------------------------------- //
//  overloaded addfields below, one for every data type
// ---------------------------------------------------------- //

void genStruct::AddField(double* x, string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}

void genStruct::AddField(vector<double> * x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}

void genStruct::AddField(int* x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}

void genStruct::AddField(short * x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}
void genStruct::AddField(long long * x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}
void genStruct::AddField(unsigned int * x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}
void genStruct::AddField(unsigned short * x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}
void genStruct::AddField(unsigned long long * x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}
void genStruct::AddField(bool * x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}
//void genStruct::AddField(vector<bool> * x,string fieldName,string matrixDimens)
/*
void genStruct::AddField(BoolRef x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
} */
void genStruct::AddField(char * x,string fieldName,string matrixDimens)
{
	if (headerMade)
		cout << "TOO LATE TO ADD FIELDS, YOU ALREADY MADE THE HEADER!!" << endl;
	else
	{
		fieldVect.push_back(genField(x,fieldName,matrixDimens,theFile));
	}
}

