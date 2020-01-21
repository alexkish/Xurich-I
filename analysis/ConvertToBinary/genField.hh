#ifndef genField_h
#define genField_h 1

using namespace std;
//typedef std::vector<bool>::reference BoolRef;
class genField
{
	public: // constructors, destructor
		genField(double* x, string fieldName,string matrixDimens, ofstream* writingFile);
		genField(vector<double> * x,string fieldName,string matrixDimens, ofstream* writingFile);
		genField(int* x, string fieldName,string matrixDimens, ofstream* writingFile);
		genField(short * x, string fieldName,string matrixDimens, ofstream* writingFile);
		genField(long long * x, string fieldName,string matrixDimens, ofstream* writingFile);
		genField(unsigned int * x, string fieldName,string matrixDimens, ofstream* writingFile);
		genField(unsigned short * x, string fieldName,string matrixDimens, ofstream* writingFile);
		genField(unsigned long long * x, string fieldName,string matrixDimens, ofstream* writingFile);
		genField(bool * x, string fieldName,string matrixDimens, ofstream* writingFile);
//		genField(vector<bool> * x, string fieldName,string matrixDimens, ofstream* writingFile);
//		genField(BoolRef x, string fieldName,string matrixDimens, ofstream* writingFile);
		genField(char * x, string fieldName,string matrixDimens, ofstream* writingFile);

		~genField();
		
	public: // public member functions
		void Write();
		string GetHeaderPart();
		void Clear();

	public: // public variables
		bool isVect;
		int maxSize;
		string Name;
		
	private: //private variables
		union fieldData {		
			double * dbData;
			vector<double> * vect_dbData;
			int * intData;
			short * shortData;
			long long * longlongData;
			unsigned int * uintData;
			unsigned short int * ushortintData;
			unsigned long long int * ulonglongData;
			bool * boolData;
			//vector<bool> * vect_boolData;
			//BoolRef vect_boolData;
			char * charData;
		} inputData;
		string headerPart;
		ofstream* outF;
		int dataType;
		unsigned short vectLength;
};

#endif
