/*
	dat2tree
	
	Usage: dat2tree -i <inputfile> [-o <outputfile>] [-bhv]
		where <inputfile> is of the *.dat format (from c++ genStruct or matlab WriteBinary.m), 
		and <outputfile> is a root binary.
		
		Options:
			-i : Input file name (of .dat format).
			-o : (optional, with argument) Output file name.  Default: Take the input file name, 
				 remove the '.dat' extension, and append '.root'.
			-b : (optional, no argument) Specifies that the input file is written in a big-endian 
				 format.  Default: little-endian.
			-h : (optional, no argument) Print help info and quit.
			-v : (optional, no argument) Print the version of the program and quit.
		
		Examples:
			dat2tree -i x01_20110522T1993.dat -o x01_20110522T1993.root
			dat2tree -i x01_20110522T1993_00001.dat -o x01_20110522T1993_00001.root -b
			dat2tree -v
	
	v01: Initial try, no support for vectors yet.
	
	2012.02.06  -- A. Manalaysay
*/

#include <iostream>
#include <stdint.h>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include <vector>
#include <sstream>
#include <unistd.h>

using namespace std;

// Function declarations to swap bytes -- need this for reading big-endian files. These functions 
// are included in up-to-date c++ libraries, but don't seem to exist on idarkx.
int16_t  _bswap(int16_t n);  // "bswap" means "byte swap".
uint16_t _bswap(uint16_t n);
int32_t  _bswap(int32_t n);
uint32_t _bswap(uint32_t n);
float    _bswap(float n);
double   _bswap(double n);

// Help declaration
void dat2treeHelp();

int main(int argc,char** argv)
{
	string file_in;
	string file_out;
	bool flag_setInFile = false;
	bool flag_setOutFile = false;
	bool flag_bigEndian = false;
	char c=0;

// Read and parse input flags and command-line options
	while((c=getopt(argc,argv,"i:o:bhv"))!=-1)
	{
		switch(c)
		{
			case 'i':
				file_in = optarg;
				flag_setInFile = true;
				break;
			case 'o':
				file_out = optarg;
				flag_setOutFile = true;
				break;
			case 'b':
				flag_bigEndian = true;
				cout << "Assuming input file is big-endian..." << endl;
				break;
			case 'h':
				dat2treeHelp();
				return 0;
			case 'v':
				cout << endl << "dat2tree version 1.0" << endl
					 << "	genStruct *.dat to root converter" << endl << endl;
				return 0;
			default:
				cout << "Unrecognized option: " << c << endl;
				return 0;
		};
	};
	
	if(!flag_setInFile)
	{
		cout << "Usage: dat2tree -i <inputfile> [-o <outputfile>] [-bhv]" << endl;
		return 1;
	};
	if(!flag_setOutFile)
	{
		int32_t pointPosition = file_in.find_last_of('.');
		string file_out_try;
		if(pointPosition!=string::npos)
		{
			file_out = file_in.substr(0,pointPosition) + ".root";
		}
		else
		{
			file_out = file_in + ".root";
		}
	};
		
// Open the file
	fstream BinaryFile;
	BinaryFile.open(file_in.c_str(),ios::binary|ios::in);
	if(!BinaryFile.is_open())
	{
		cout << "Could not open the file: " << file_in << endl << endl;
		return 1;
	};
	
// Read the header length and header
	uint16_t HeaderSize;
	BinaryFile.seekg(0,ios::beg);
	BinaryFile.read((char*)&HeaderSize,sizeof(uint16_t));
	if(flag_bigEndian)
	{
		HeaderSize = _bswap(HeaderSize);
	};
	if(HeaderSize>2000)
	{
		cout << "WARNING: size of header is read as >2000 bytes.  It is possible that you" << endl
			 << "are using the wrong endian; try toggling the -b option." << endl;
	};
	char HeaderChar[HeaderSize];
	string Header;
	BinaryFile.read(HeaderChar,HeaderSize*sizeof(char));
	Header = HeaderChar;
	Header.resize(HeaderSize);
	if(Header[0]!=';')
	{
		Header = ";" + Header;
	};	
	
// Search header for semicolon deliminators
	int32_t delim_lastSemi;
	delim_lastSemi = Header.find_last_of(";");
	if(delim_lastSemi>HeaderSize)
	{
		cout << "Bad file header." << endl << endl;
		return 1;
	};
	
	vector<int32_t> HeaderDelims;
	int32_t HeaderDelimPos = -1;
	while(HeaderDelimPos<HeaderSize)
	{
		HeaderDelimPos = Header.find(';',HeaderDelimPos+1);
		HeaderDelims.push_back(HeaderDelimPos);
	};
	
// Parse the header, taking information between semicolons
	vector<string> FieldNames;
	vector<string> FieldTypes;
	vector<string> FieldDim1s;
	vector<int32_t> FieldDim1;
	vector<string> FieldDim2s;
	vector<int32_t> FieldDim2;
	int32_t commaPosition;
	string FieldDimsString;
	string FieldDimIndvS;
	int32_t FieldDimIndv;
	stringstream FieldDimsStream;
	
	for(int i=0;i<(HeaderDelims.size()-1);i+=3)
	{
		// Read info for field name
		FieldNames.push_back(Header.substr(HeaderDelims[i]+1,HeaderDelims[i+1]-HeaderDelims[i]-1));
		// Read info for field type
		FieldTypes.push_back(Header.substr(HeaderDelims[i+1]+1,HeaderDelims[i+2]-HeaderDelims[i+1]-1));
		// Read info for field dimensions and parse
		FieldDimsString = Header.substr(HeaderDelims[i+2]+1,HeaderDelims[i+3]-HeaderDelims[i+2]-1);
		commaPosition = FieldDimsString.find(",");
		FieldDimIndvS = FieldDimsString.substr(0,commaPosition);
		FieldDimsStream << FieldDimIndvS;
		FieldDimsStream >> FieldDimIndv;
		FieldDim1.push_back(FieldDimIndv);
		FieldDim1s.push_back(FieldDimIndvS);
		FieldDimsStream.clear();
		FieldDimIndvS = FieldDimsString.substr(commaPosition+1,FieldDimsString.size());
		FieldDimsStream << FieldDimIndvS;
		FieldDimsStream >> FieldDimIndv;
		FieldDim2.push_back(FieldDimIndv);
		FieldDim2s.push_back(FieldDimIndvS);
		FieldDimsStream.clear();		
	};
	// I now have the following vectors containing information about the fields ("fields"=="branches")
	// FieldNames (string)
	// FieldTypes (string)
	// FieldDim1  (int32_t)
	// FieldDim1s (string)
	// FieldDim2  (int32_t)
	// FieldDim2s (string)

// Create the array of pointers that will hold the data pointers
	int32_t NumFields = FieldNames.size();
	char** FieldDataPointersArray = (char**)malloc(NumFields);
	string FieldLeafList;

// Construct the root file and the tree
	TFile hfile(file_out.c_str(),"RECREATE","File with ROOT tree");
	TTree d("d","Data_tree"); // "d" for "data", what I call it in matlab anyways.

	vector<int32_t> FieldSizeBytes;
	int32_t FieldsTotalSizeBytes = 0;
	
// Initialize the ROOT tree
	for(int32_t i=0;i<NumFields;i++)
	{
	//Construct the leaflist, with appropriate array dimensions (but without the type specifier)
		if((FieldDim1[i]<1)||(FieldDim2[i]<1))
		{
			cout << "Field '" << FieldNames[i] << "' has improper dimensions." << endl << endl;
			if(i>0)
			{
				for(int32_t i_e=0;i_e<i;i_e++)
				{
					free(FieldDataPointersArray[i_e]);
				};
			};
			free(FieldDataPointersArray);
			return 1;
		}
		else if((FieldDim1[i]==1)&&(FieldDim2[i]==1))
		{
			FieldLeafList = FieldNames[i];
		}
		else if((FieldDim1[i]==1)&&(FieldDim2[i]>1))
		{
			FieldLeafList = FieldNames[i] + "[" + FieldDim2s[i] + "]";
		}
		else if((FieldDim1[i]>1)&&(FieldDim2[i]==1))
		{
			FieldLeafList = FieldNames[i] + "[" + FieldDim1s[i] + "]";
		}
		else
		{
			FieldLeafList = FieldNames[i] + "[" + FieldDim1s[i] + "][" + FieldDim2s[i] + "]";
		};
		
	// Apply the appopriate data type and allocate the memory in the pointer array.
		if(FieldTypes[i]=="char")
		{
			FieldSizeBytes.push_back(sizeof(char)*FieldDim1[i]*FieldDim2[i]);
			FieldsTotalSizeBytes += FieldSizeBytes[i];
			FieldDataPointersArray[i] = (char*)malloc(FieldSizeBytes[i]);
			FieldLeafList = FieldLeafList + "/B";
			d.Branch(FieldNames[i].c_str(),(char*)FieldDataPointersArray[i],FieldLeafList.c_str());
		}
		else if(FieldTypes[i]=="logical")
		{
			FieldSizeBytes.push_back(sizeof(bool)*FieldDim1[i]*FieldDim2[i]);
			FieldsTotalSizeBytes += FieldSizeBytes[i];
			FieldDataPointersArray[i] = (char*)malloc(FieldSizeBytes[i]);
			FieldLeafList = FieldLeafList + "/O";
			d.Branch(FieldNames[i].c_str(),(bool*)FieldDataPointersArray[i],FieldLeafList.c_str());
		}
		else if(FieldTypes[i]=="int16")
		{
			FieldSizeBytes.push_back(sizeof(int16_t)*FieldDim1[i]*FieldDim2[i]);
			FieldsTotalSizeBytes += FieldSizeBytes[i];
			FieldDataPointersArray[i] = (char*)malloc(FieldSizeBytes[i]);
			FieldLeafList = FieldLeafList + "/S";
			d.Branch(FieldNames[i].c_str(),(int16_t*)FieldDataPointersArray[i],FieldLeafList.c_str());
		}
		else if(FieldTypes[i]=="uint16")
		{
			FieldSizeBytes.push_back(sizeof(uint16_t)*FieldDim1[i]*FieldDim2[i]);
			FieldsTotalSizeBytes += FieldSizeBytes[i];
			FieldDataPointersArray[i] = (char*)malloc(FieldSizeBytes[i]);
			FieldLeafList = FieldLeafList + "/s";
			d.Branch(FieldNames[i].c_str(),(uint16_t*)FieldDataPointersArray[i],FieldLeafList.c_str());
		}
		else if(FieldTypes[i]=="int32")
		{
			FieldSizeBytes.push_back(sizeof(int32_t)*FieldDim1[i]*FieldDim2[i]);
			FieldsTotalSizeBytes += FieldSizeBytes[i];
			FieldDataPointersArray[i] = (char*)malloc(FieldSizeBytes[i]);
			FieldLeafList = FieldLeafList + "/I";
			d.Branch(FieldNames[i].c_str(),(int32_t*)FieldDataPointersArray[i],FieldLeafList.c_str());
		}
		else if(FieldTypes[i]=="uint32")
		{
			FieldSizeBytes.push_back(sizeof(uint32_t)*FieldDim1[i]*FieldDim2[i]);
			FieldsTotalSizeBytes += FieldSizeBytes[i];
			FieldDataPointersArray[i] = (char*)malloc(FieldSizeBytes[i]);
			FieldLeafList = FieldLeafList + "/i";
			d.Branch(FieldNames[i].c_str(),(uint32_t*)FieldDataPointersArray[i],FieldLeafList.c_str());
		}
		else if(FieldTypes[i]=="float")
		{
			FieldSizeBytes.push_back(sizeof(float)*FieldDim1[i]*FieldDim2[i]);
			FieldsTotalSizeBytes += FieldSizeBytes[i];
			FieldDataPointersArray[i] = (char*)malloc(FieldSizeBytes[i]);
			FieldLeafList = FieldLeafList + "/F";
			d.Branch(FieldNames[i].c_str(),(float*)FieldDataPointersArray[i],FieldLeafList.c_str());
		}
		else if(FieldTypes[i]=="double")
		{
			FieldSizeBytes.push_back(sizeof(double)*FieldDim1[i]*FieldDim2[i]);
			FieldsTotalSizeBytes += FieldSizeBytes[i];
			FieldDataPointersArray[i] = (char*)malloc(FieldSizeBytes[i]);
			FieldLeafList = FieldLeafList + "/D";
			d.Branch(FieldNames[i].c_str(),(double*)FieldDataPointersArray[i],FieldLeafList.c_str());
		}
		else
		{
			cout << "I don't recognize data type '" << FieldTypes[i] << "', please use supported types only." 
				 << endl << endl;
			return 1;
		};
	};
		
// Calcuate the number of events in the .dat file
	uint32_t FileSize;
	BinaryFile.seekg(0,ios::end);
	FileSize = BinaryFile.tellg();
	BinaryFile.seekg(2+HeaderSize,ios::beg);
	uint32_t NumEvents = (FileSize-HeaderSize-2)/FieldsTotalSizeBytes;

// Read through the .dat file and fill the tree	
	for(uint32_t i_0=0;i_0<NumEvents;i_0++) // loop over events
	{
		for(uint32_t i_1=0;i_1<NumFields;i_1++) // loop over fields
		{
			BinaryFile.read(FieldDataPointersArray[i_1],FieldSizeBytes[i_1]);

			if(flag_bigEndian)
			{
				for(uint32_t i_2=0;i_2<(FieldDim1[i_1]*FieldDim2[i_1]);i_2++)
				{
					if(FieldTypes[i_1]=="int16")
					{
						((int16_t*)FieldDataPointersArray[i_1])[i_2] = _bswap(((int16_t*)FieldDataPointersArray[i_1])[i_2]);
					}
					else if(FieldTypes[i_1]=="uint16")
					{
						((uint16_t*)FieldDataPointersArray[i_1])[i_2] = _bswap(((uint16_t*)FieldDataPointersArray[i_1])[i_2]);
					}
					else if(FieldTypes[i_1]=="int32")
					{
						((int32_t*)FieldDataPointersArray[i_1])[i_2] = _bswap(((int32_t*)FieldDataPointersArray[i_1])[i_2]);
					}
					else if(FieldTypes[i_1]=="uint32")
					{
						((uint32_t*)FieldDataPointersArray[i_1])[i_2] = _bswap(((uint32_t*)FieldDataPointersArray[i_1])[i_2]);
					}
					else if(FieldTypes[i_1]=="float")
					{
						((float*)FieldDataPointersArray[i_1])[i_2] = _bswap(((float*)FieldDataPointersArray[i_1])[i_2]);
					}
					else if(FieldTypes[i_1]=="double")
					{
						((double*)FieldDataPointersArray[i_1])[i_2] = _bswap(((double*)FieldDataPointersArray[i_1])[i_2]);
					}
				}
			}
		};
		d.Fill();
	};
	
// Write to the root file, close it, close the .dat file, free them from heap memory
//
	hfile.Write();
	hfile.Close();
	BinaryFile.close();
	for(int32_t i=0;i<NumFields;i++)
	{
		free(FieldDataPointersArray[i]);
	};
	free(FieldDataPointersArray);

// Return zero for happy ending
	return 0;	
};
// ---------------------END OF MAIN FUNCTION--------------------- //

// ....
// ....

// -------------------=======================-------------------- //
// ---------------------AUXILIARY FUNCTIONS---------------------- //
// -------------------=======================-------------------- //
int16_t  _bswap(int16_t n)
{
	int16_t nn = n;
	char* a = (char*)&nn;
	char b[sizeof(int16_t)];
	for(int32_t i=0;i<sizeof(int16_t);i++)
	{
		b[i] = a[sizeof(int16_t)-1-i];
	};
	return *((int16_t*)b);	
};
uint16_t _bswap(uint16_t n)
{
	uint16_t nn = n;
	char* a = (char*)&nn;
	char b[sizeof(uint16_t)];
	for(int32_t i=0;i<sizeof(uint16_t);i++)
	{
		b[i] = a[sizeof(uint16_t)-1-i];
	};
	return *((uint16_t*)b);	
};
int32_t  _bswap(int32_t n)
{
	int32_t nn = n;
	char* a = (char*)&nn;
	char b[sizeof(int32_t)];
	for(int32_t i=0;i<sizeof(int32_t);i++)
	{
		b[i] = a[sizeof(int32_t)-1-i];
	};
	return *((int32_t*)b);	
};
uint32_t _bswap(uint32_t n)
{
	uint32_t nn = n;
	char* a = (char*)&nn;
	char b[sizeof(uint32_t)];
	for(int32_t i=0;i<sizeof(uint32_t);i++)
	{
		b[i] = a[sizeof(uint32_t)-1-i];
	};
	return *((int32_t*)b);	
}
float    _bswap(float n)
{
	float nn = n;
	char* a = (char*)&nn;
	char b[sizeof(float)];
	for(int32_t i=0;i<sizeof(float);i++)
	{
		b[i] = a[sizeof(float)-1-i];
	};
	return *((float*)b);	
};
double   _bswap(double n)
{
	double nn = n;
	char* a = (char*)&nn;
	char b[sizeof(double)];
	for(int32_t i=0;i<sizeof(double);i++)
	{
		b[i] = a[sizeof(double)-1-i];
	};
	return *((double*)b);	
};

void dat2treeHelp()
{
	cout << endl
		 << "dat2tree" << endl << endl
		 << "Usage: dat2tree -i <inputfile> [-o <outputfile>] [-bhv]" << endl
		 << "	where <inputfile> is of the *.dat format (from c++ genStruct or matlab WriteBinary.m), " << endl
		 << "	and <outputfile> is a root binary." << endl << endl
		 << "	Options:" << endl
		 << "		-i : Input file name (of .dat format)." << endl
		 << "		-o : (optional, with argument) Output file name.  Default: Take the input file name," << endl
		 << "			 remove the '.dat' extension, and append '.root'." << endl
		 << "		-b : (optional, no argument) Specifies that the input file is written in a big-endian" << endl
		 << "			 format.  Default: little-endian." << endl
		 << "		-h : (optional, no argument) Print help info and quit." << endl
		 << "		-v : (optional, no argument) Print the version of the program and quit." << endl << endl
		 << "	Examples:" << endl
		 << "		dat2tree -i x01_20110522T1993.dat -o x01_20110522T1993.root" << endl
		 << "		dat2tree -i x01_20110522T1993_00001.dat -o x01_20110522T1993_00001.root -b" << endl
		 << "		dat2tree -v" << endl << endl
		 << "v01: Initial try, no support for vectors yet." << endl << endl
		 << "2012.02.06 -- A. Manalaysay" << endl << endl;
};