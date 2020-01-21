//====================================================================
// ROOT macro to read binary files from Xuerich procesed data
// and to produce a ROOT tree out of it
//====================================================================
// Teresa Marrodan Undagoitia, 27/02/09
// Modified to be compiled directly (not with CINT from ROOT) 18/06/09
// Working with a first data file 01/09/09
// Modified because Tree branches do not accept matrizes 08/10/09
// Modified to correct an error found by Aaron on the S2 branches 23/11/2010
// Comments added to make the program easier for others 26/11/2010
// added part to convert data with TAC and NaI area 02/03/2011

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include <vector>

using namespace std;

int main(int argc,char** argv){

  if (argc == 3){    
    cout << "Data file : " << argv[1] << "; Output file : " << argv[2]<< endl;
  }
  else {
    cout << "You have to give a file name and output name!" << endl;
    exit(0);
  }
    

// Compile: g++ `root-config --glibs` `root-config --cflags` 
// V2_Comp_ROOTreader.C -o V2_Comp_ROOTreader
// Run: ./V2_Comp_ROOTreader Name_of_file.dat


// === === === Reading the header === === === 

  // The first two bites contain the length of the header in a 'unsigned 
  // short' variable. The header has variable name, variable type and 
  // dimensions of the data matrix, all separated by colons or semicolons.
  // Example: 'channel1;int16;1,100;channel2;int16;1,100;'
 
  char HeaderSize[100]={};
  char Header[1000]={};
  char H_Characters[1000]={};
  unsigned short H_Size=0;
  vector<string> VariableName;
  vector<string> VariableType;
  vector<int> Dimm1;
  vector<int> Dimm2;

  VariableName.clear();
  VariableType.clear();
  Dimm1.clear();
  Dimm2.clear();

  //cout << "Size of double = " << sizeof(double)<< endl;
  
  vector<string> Types;
  Types.push_back("char");
  Types.push_back("logical");
  Types.push_back("int16");
  Types.push_back("int32");
  Types.push_back("uint32");
  Types.push_back("double"); 
  int TypesNum[6]={1, 1, 2, 4, 4, 8};
  

// Open the file --------- ---------
  fstream BinaryFile;
  BinaryFile.open(argv[1],ios::binary|ios::in); //File is an argument
  //BinaryFile.open("x01_20090817T1711.dat",ios::binary|ios::in);
  if (BinaryFile.is_open()){cout << "File open!"<< endl;}
 

// Read the header size and header --------- ---------
  BinaryFile.seekg(0, ios::beg);
  BinaryFile.read((char*)&H_Size,2);// to read directly a short variable
  cout << "HeaderSize= "<< H_Size << endl;
  
  BinaryFile.seekg(2, ios::beg); //jump two bytes beyond the beginning
  BinaryFile.read(Header,H_Size); //read bytes corresponding to the header
  cout << "Header = "<< Header << endl;
  //ofstream myFile ("Output.txt");// Write the header to a file
  //myFile.write(Header,H_Size);


  // Search for ";" or "," in the header --------- ---------
  BinaryFile.seekg(2, ios::beg);
  
  int PosLim=2; 
  int counter=0;
  int counter1=0;
  int counter2=0;
  int counter3=0;
  int Dimm_1=0;
  int Dimm_2=0;

  for(int i=3;i<=2+H_Size;i++){
    BinaryFile.read(H_Characters,1);    
    char bla[100]={};
   
    if (H_Characters[0] == 59 || H_Characters[0] == 44){
      //cout << i << "-" << H_Characters << endl;
      BinaryFile.seekg(PosLim, ios::beg);
      BinaryFile.read(bla,(i-PosLim-1));
      PosLim=i;
      
      if((counter%4)==0){
	VariableName.push_back(bla);
      }      counter+=1;

      if((counter1%4)==1){
	VariableType.push_back(bla);
      }      counter1+=1;

      if((counter2%4)==2){
	Dimm_1=atoi(bla);
      	Dimm1.push_back(Dimm_1);
      }      counter2+=1;

      if((counter3%4)==3){
	Dimm_2=atoi(bla);
      	Dimm2.push_back(Dimm_2);
      }      counter3+=1;

    }
    BinaryFile.seekg(i, ios::beg);
  }
  
  cout << endl << "Content of the vector VariableName: " << endl
       << "----------" << endl; 
  for (int i=0; i< VariableName.size(); i++) {
   // cout << "Te quiero mucho, mi Maus! " << i << endl;
    cout << "VariableName "<< ": " << VariableName[i] << endl;
    cout << "VariableType "<< ": " << VariableType[i] << endl;
    cout << "Dimm1 "<< ": " << Dimm1[i] << endl;
    cout << "Dimm2 " << ": " << Dimm2[i] << endl;
  }
  cout << "----------" << endl;
  //cout << "Size of VariableName = " << VariableName.size() <<endl; 


  // Size of the file --------- ----------
  int begin,end;
  BinaryFile.seekg(0, ios::beg);
  begin = BinaryFile.tellg();
  BinaryFile.seekg (0, ios::end);
  end = BinaryFile.tellg();
  cout << "The file size is: " << (end-begin) << " bytes long " << endl;


// Compare variable types and assign byte-length --------- ----------
  int Comp[VariableName.size()][Types.size()]; 
  int Length[VariableName.size()];  
  for(int i=0;i<VariableName.size();i++){
    for(int j=0;j<(Types.size());j++){
      Comp[i][j]= VariableType[i].compare(Types[j]);
      if(Comp[i][j]==0){Length[i]=TypesNum[j];
      }
    }
  }
  cout << "-----------"<<endl;
 

  // Calculate the number of events (total_length/event_length)  ---------
  int N_events = 0;
  int E_length = 0; // Total length of one event
  for (int i=0;i<VariableName.size();i++){
    E_length += Length[i]*Dimm1[i]*Dimm2[i];
  }
  N_events = (end-begin-H_Size-2)/(E_length);
  cout<< "N_events: "<< N_events << endl;



  // === === === === === File and tree definition === === === === 

  TFile *hfile = new TFile(argv[2],"RECREATE","File with ROOT tree");
  TTree *tree = new TTree("tree","XuerichData_Tree"); // Tree definition

 

 // --- Create the tree if there are S2s in the data  --------- ----------

  if (VariableName[0]=="evtNum" && VariableName[1]=="pA_S2" && 
      VariableName[2]=="pH_S2" && VariableName[3]=="S2_mask"&& 
      VariableName[4]=="S2"&& VariableName[5]=="t0lf_S2"&& 
      VariableName[6]=="t50lf_S2" && VariableName[7]=="t50rf_S2" && 
      VariableName[8]=="t0rf_S2" && VariableName[9]=="numS2" && 
      VariableName[10]=="pA_S1" && VariableName[11]=="pH_S1" && 
      VariableName[12]=="t0l_S1" && VariableName[13]=="t10l_S1" && 
      VariableName[14]=="t50l_S1" && VariableName[15]=="t1_S1" && 
      VariableName[16]=="t50r_S1" && VariableName[17]=="t10r_S1" && 
      VariableName[18]=="t0r_S1" && VariableName[19]=="S1_mask" && 
      VariableName[20]=="numS1" && VariableName[21]=="S1" && 
      VariableName[22]=="dT" && VariableName[23]=="BS_ave" && 
      VariableName[24]=="BS_std" && VariableName[25]=="Sat"){

    cout  << "----------" << endl
	  << "XUERICH DATA FORMAT: with Field (S2)" <<endl
	  << "----------" << endl;

 // Declaration of the variables of the tree  --------- ----------

      unsigned int evtNum[1]={};

      double pA_S2_f[2]={};
      double pA_S2_s[2]={};
      double pA_S2_t[2]={};
      double pH_S2_f[2]={};
      double pH_S2_s[2]={};
      double pH_S2_t[2]={};

      bool S2_mask[3]={};
      double S2[3]={};
      double t0lf_S2[3]={};
      double t50lf_S2[3]={};
      double t50rf_S2[3]={};
      double t0rf_S2[3]={};
      double numS2[1]={};

      double pA_S1_f[2]={};
      double pA_S1_s[2]={};
      double pH_S1_f[2]={};
      double pH_S1_s[2]={};

      double t0l_S1[2]={};
      double t10l_S1[2]={};
      double t50l_S1[2]={};
      double t1_S1[2]={};
      double t50r_S1[2]={};
      double t10r_S1[2]={};
      double t0r_S1[2]={};
      bool S1_mask[2]={};
      double numS1[1]={};
      double S1[2]={};
      double dT[3]={};
      double BS_ave[2]={};
      double BS_std[2]={};
      bool Sat[2]={};

      // Declaration of the branches of the tree "name", address, variable/type ----- ----- 

      TBranch *MyBranch_1 = tree->Branch("evtNum", &evtNum, 
					 "evtNum[1]/i");
      TBranch *MyBranch_2_f = tree->Branch("pA_S2_f", &pA_S2_f, 
					 "pA_S2_f[2]/D");
      TBranch *MyBranch_2_s = tree->Branch("pA_S2_s", &pA_S2_s, 
					 "pA_S2_s[2]/D");
      TBranch *MyBranch_2_t = tree->Branch("pA_S2_t", &pA_S2_t, 
					 "pA_S2_t[2]/D");
      TBranch *MyBranch_3_f = tree->Branch("pH_S2_f", &pH_S2_f, 
					 "pH_S2_f[2]/D");
      TBranch *MyBranch_3_s = tree->Branch("pH_S2_s", &pH_S2_s, 
					 "pH_S2_s[2]/D");
      TBranch *MyBranch_3_t = tree->Branch("pH_S2_t", &pH_S2_t, 
					 "pH_S2_t[2]/D");
      TBranch *MyBranch_4 = tree->Branch("S2_mask", &S2_mask, 
					 "S2_mask[3]/O");
      TBranch *MyBranch_5 = tree->Branch("S2", &S2, 
					 "S2[3]/D"); 
      TBranch *MyBranch_6 = tree->Branch("t0lf_S2", &t0lf_S2, 
					 "t0lf_S2[3]/D"); 
      TBranch *MyBranch_7 = tree->Branch("t50lf_S2", &t50lf_S2, 
					 "t50lf_S2[3]/D"); 
      TBranch *MyBranch_8 = tree->Branch("t50rf_S2", &t50rf_S2, 
					 "t50rf_S2[3]/D"); 
      TBranch *MyBranch_9 = tree->Branch("t0rf_S2", &t0rf_S2, 
					 "t0rf_S2[3]/D"); 
      TBranch *MyBranch_10 = tree->Branch("numS2", &numS2, 
					 "numS2[1]/D");
      TBranch *MyBranch_11_f = tree->Branch("pA_S1_f", &pA_S1_f, 
					"pA_S1_f[2]/D"); 
      TBranch *MyBranch_11_s = tree->Branch("pA_S1_s", &pA_S1_s, 
					"pA_S1_s[2]/D"); 
      TBranch *MyBranch_12_f = tree->Branch("pH_S1_f", &pH_S1_f, 
					"pH_S1_f[2]/D"); 
      TBranch *MyBranch_12_s = tree->Branch("pH_S1_s", &pH_S1_s, 
					"pH_S1_s[2]/D"); 
      TBranch *MyBranch_13 = tree->Branch("t0l_S1", &t0l_S1, 
					"t0l_S1[2]/D"); 
      TBranch *MyBranch_14 = tree->Branch("t10l_S1", &t10l_S1, 
					"t10l_S1[2]/D"); 
      TBranch *MyBranch_15 = tree->Branch("t50l_S1", &t50l_S1, 
					  "t50l_S1[2]/D"); 
      TBranch *MyBranch_16 = tree->Branch("t1_S1", &t1_S1, 
					  "t1_S1[2]/D"); 
      TBranch *MyBranch_17 = tree->Branch("t50r_S1", &t50r_S1, 
					  "t50r_S1[2]/D"); 
      TBranch *MyBranch_18 = tree->Branch("t10r_S1", &t10r_S1, 
					  "t10r_S1[2]/D");
      TBranch *MyBranch_19 = tree->Branch("t0r_S1", &t0r_S1, 
					  "t0r_S1[2]/D");
      TBranch *MyBranch_20 = tree->Branch("S1_mask", &S1_mask, 
					  "S1_mask[3]/O");
      TBranch *MyBranch_21 = tree->Branch("numS1", &numS1, 
					  "numS1[1]/D");
      TBranch *MyBranch_22 = tree->Branch("S1", &S1, 
					  "S1[2]/D");
      TBranch *MyBranch_23 = tree->Branch("dT", &dT, 
					  "dT[3]/D");
      TBranch *MyBranch_24 = tree->Branch("BS_ave", &BS_ave, 
					  "BS_ave[2]/D");
      TBranch *MyBranch_25 = tree->Branch("BS_std", &BS_std, 
					  "BS_std[2]/D");
      TBranch *MyBranch_26 = tree->Branch("Sat", &Sat, 
					  "Sat[2]/O");


      // --- Loop over the data to fill the variable with the numbers ------ ------

      for(int i=0;i<N_events;i++){ // --- loop over number of event ---
	if((i%1000)==0){cout << "Event = " << i <<endl;}
	int Variable_Length[1000000]={};

	for (int j=0; j<(VariableName.size());j++){// --- loop over the variables ---
	  for(int k=0; k<Dimm1[j];k++){ // --- loop over the first dimension
	    for(int l=0; l<Dimm2[j];l++){ // --- loop over the second dimension

	      // --- initiate the variable to avoid errors --- 
	      unsigned int evtNum_v=0;
	      double pA_S2_v=-1;
	      double pH_S2_v=-1;
	      bool S2_mask_v=-1;
	      double S2_v=-1;
	      double t0lf_S2_v=-1;
	      double t50lf_S2_v=-1;
	      double t50rf_S2_v=-1;
	      double t0rf_S2_v=-1;
	      double numS2_v=-1;
	      double pA_S1_v=-1;
	      double pH_S1_v=-1;
	      double t0l_S1_v=-1;
	      double t10l_S1_v=-1;
	      double t50l_S1_v=-1;
	      double t1_S1_v=-1;
	      double t50r_S1_v=-1;
	      double t10r_S1_v=-1;
	      double t0r_S1_v=-1;
	      bool S1_mask_v=-1;
	      double numS1_v=-1;
	      double S1_v=-1;
	      double dT_v=-1;
	      double BS_ave_v=-1;
	      double BS_std_v=-1;
	      bool Sat_v=-1;


	      Variable_Length[j]=(Length[j])*(Dimm1[j])*(Dimm2[j]);
	      Variable_Length[j]= Variable_Length[j]+Variable_Length[j-1];


	      BinaryFile.seekg(2 + H_Size 
			       + (k*Dimm2[j]*(Length[j]))
			       + (l*(Length[j]))
			       + (Variable_Length[j-1]) 
			       + (i*E_length), ios::beg);

	 
	      if(VariableName[j]=="evtNum"){
		BinaryFile.read((char*)(&evtNum_v),Length[j]);
		//cout << "----------" << endl
		//     << "Event Nummer = " << evtNum_v << endl;
		evtNum[l]=evtNum_v;
	      }

	      if(VariableName[j]=="pA_S2"){
	      	BinaryFile.read((char*)(&pA_S2_v),Length[j]);
		if(l==0){pA_S2_f[k]=pA_S2_v; // First S2
		}
		if(l==1){pA_S2_s[k]=pA_S2_v; // Second S2
		}
		if(l==2){pA_S2_t[k]=pA_S2_v; // Third S2
		}
	      }

	      if(VariableName[j]=="pH_S2"){
	      	BinaryFile.read((char*)(&pH_S2_v),Length[j]);
		if(l==0){pH_S2_f[k]=pH_S2_v; // First S2
		}
		if(l==1){pH_S2_s[k]=pH_S2_v; // Second S2
		}
		if(l==2){pH_S2_t[k]=pH_S2_v; // Third S2
		}
	      }

	      if(VariableName[j]=="S2_mask"){
	      	BinaryFile.read((char*)(&S2_mask_v),Length[j]);
		S2_mask[l]=S2_mask_v;
	      }

	      if(VariableName[j]=="S2"){
	      	BinaryFile.read((char*)(&S2_v),Length[j]);
		S2[l]=S2_v;
	      }

	      if(VariableName[j]=="t0lf_S2"){
	      	BinaryFile.read((char*)(&t0lf_S2_v),Length[j]);
		t0lf_S2[l]=t0lf_S2_v;
	      }

	      if(VariableName[j]=="t50lf_S2"){
	      	BinaryFile.read((char*)(&t50lf_S2_v),Length[j]);
		t50lf_S2[l]=t50lf_S2_v;
	      }
 
	      if(VariableName[j]=="t50rf_S2"){
	      	BinaryFile.read((char*)(&t50rf_S2_v),Length[j]);
		t50rf_S2[l]=t50rf_S2_v;
	      }

	      if(VariableName[j]=="t0rf_S2"){
	      	BinaryFile.read((char*)(&t0rf_S2_v),Length[j]);
		t0rf_S2[l]=t0rf_S2_v;
	      }

	      if(VariableName[j]=="numS2"){
		BinaryFile.read((char*)(&numS2_v),Length[j]);
	      	numS2[l]=numS2_v;
		//if(numS2_v!=0){cout<< "numS2_v = " << numS2_v <<endl;}
	      }

	      if(VariableName[j]=="numS1"){
	      	BinaryFile.read((char*)(&numS1_v),Length[j]);
		numS1[l]=numS1_v;
	      } 

	      if(VariableName[j]=="pA_S1"){
	      	BinaryFile.read((char*)(&pA_S1_v),Length[j]);	       
		if(l==0){pA_S1_f[k]=pA_S1_v; //First S1
	   //cout <<"k = "<<k<<" "<<l<<" S1 area_f : "<<pA_S1_f[l]<< endl;
		}
		if(l==1){pA_S1_s[k]=pA_S1_v; //Second S1
	   //cout <<"k = "<<k<<" "<<l<<" S1 Pulse area_s : "<<pA_S1_s[l]<< endl;
		}		
	      }

	      if(VariableName[j]=="pH_S1"){
	      	BinaryFile.read((char*)(&pH_S1_v),Length[j]);
		if(l==0){pH_S1_f[k]=pH_S1_v; //First S1
		}
		if(l==1){pH_S1_s[k]=pH_S1_v; //Second S1
		}
	      }

	      if(VariableName[j]=="t0l_S1"){
		BinaryFile.read((char*)(&t0l_S1_v),Length[j]);
	      	t0l_S1[l]=t0l_S1_v;
	      }

	      if(VariableName[j]=="t10l_S1"){
		BinaryFile.read((char*)(&t10l_S1_v),Length[j]);
	      	t10l_S1[l]=t10l_S1_v;
	      }

	      if(VariableName[j]=="t50l_S1"){
		BinaryFile.read((char*)(&t50l_S1_v),Length[j]);
	      	t50l_S1[l]=t50l_S1_v;
	      }

	      if(VariableName[j]=="t1_S1"){
		BinaryFile.read((char*)(&t1_S1_v),Length[j]);
	      	t1_S1[l]=t1_S1_v;
		//if(numS1[0]==2){
		  //cout <<" l = "<<l<<" "<<" t1_S1 : "<<t1_S1[l]<< endl;}
	      }

	      if(VariableName[j]=="t0r_S1"){
		BinaryFile.read((char*)(&t0r_S1_v),Length[j]);
	      	t0r_S1[l]=t0r_S1_v;
	      }

	      if(VariableName[j]=="t10r_S1"){
		BinaryFile.read((char*)(&t10r_S1_v),Length[j]);
	      	t10r_S1[l]=t10r_S1_v;
	      }

	      if(VariableName[j]=="t50r_S1"){
		BinaryFile.read((char*)(&t50r_S1_v),Length[j]);
	      	t50r_S1[l]=t50r_S1_v;
	      }

	      if(VariableName[j]=="S1_mask"){
	      	BinaryFile.read((char*)(&S1_mask_v),Length[j]);
		S1_mask[l]=S1_mask_v;
	      } 
 
	      if(VariableName[j]=="S1"){
	      	BinaryFile.read((char*)(&S1_v),Length[j]);
		S1[l]=S1_v;
	      //cout <<" k = "<<k<<"; l = "<<l<<" ---- S1 : "<<S1[l]<< endl;
	      } 

	      if(VariableName[j]=="dT"){
	      	BinaryFile.read((char*)(&dT_v),Length[j]);
		dT[l]=dT_v;
	      } 

	      if(VariableName[j]=="BS_ave"){
	      	BinaryFile.read((char*)(&BS_ave_v),Length[j]);
		BS_ave[l]=BS_ave_v;
	      } 

	      if(VariableName[j]=="BS_std"){
	      	BinaryFile.read((char*)(&BS_std_v),Length[j]);
		BS_std[l]=BS_std_v;
	      } 

	      if(VariableName[j]=="Sat"){
	      	BinaryFile.read((char*)(&Sat_v),Length[j]);
		Sat[l]=Sat_v;
	      } 
 
	    }
	  }
	}// end loop over the variables


	// --- Filling the tree with the variables declared above ------ ------
	tree->Fill(); 
	
	//if (i==2000) break;

      }// end loop over the events
  }


  else{


  // === === === === Here one phase including TAC and NaI === === === === 

 if (VariableName[0]=="evtNum" && 
     VariableName[1]=="pA_S1" && VariableName[2]=="pH_S1" && 
     VariableName[3]=="t0l_S1" && VariableName[4]=="t10l_S1" && 
     VariableName[5]=="t50l_S1" && VariableName[6]=="t1_S1" && 
     VariableName[7]=="t50r_S1" && VariableName[8]=="t10r_S1" && 
     VariableName[9]=="t0r_S1" && VariableName[10]=="S1_mask" && 
     VariableName[11]=="numS1" && VariableName[12]=="S1" && 
     VariableName[13]=="BS_ave" && VariableName[14]=="BS_std" && 
     VariableName[15]=="Sat" && VariableName[16]=="TAC_ave" && 
     VariableName[17]=="NaI_Area"){
   
   cout  << "----------" << endl
	 << "XUERICH DATA FORMAT: Only S1 including TAC and NaI" <<endl
	 << "----------" << endl;
   


      unsigned int evtNum[1]={};
      double pA_S1_f[2]={};
      double pA_S1_s[2]={};
      double pH_S1_f[2]={};
      double pH_S1_s[2]={};
      double t0l_S1[2]={};
      double t10l_S1[2]={};
      double t50l_S1[2]={};
      double t1_S1[2]={};
      double t50r_S1[2]={};
      double t10r_S1[2]={};
      double t0r_S1[2]={};
      bool S1_mask[2]={};
      double numS1[1]={};
      double S1[2]={};
      double BS_ave[2]={};
      double BS_std[2]={};
      bool Sat[2]={};
      double TAC_ave[1]={};
      double NaI_Area[1]={};



      TBranch *MyBranch_1 = tree->Branch("evtNum", &evtNum, 
					 "evtNum[1]/i");
      TBranch *MyBranch_11_f = tree->Branch("pA_S1_f", &pA_S1_f, 
					"pA_S1_f[2]/D"); 
      TBranch *MyBranch_11_s = tree->Branch("pA_S1_s", &pA_S1_s, 
					"pA_S1_s[2]/D"); 
      TBranch *MyBranch_12_f = tree->Branch("pH_S1_f", &pH_S1_f, 
					"pH_S1_f[2]/D"); 
      TBranch *MyBranch_12_s = tree->Branch("pH_S1_s", &pH_S1_s, 
					"pH_S1_s[2]/D"); 
      TBranch *MyBranch_13 = tree->Branch("t0l_S1", &t0l_S1, 
					"t0l_S1[2]/D"); 
      TBranch *MyBranch_14 = tree->Branch("t10l_S1", &t10l_S1, 
					"t10l_S1[2]/D"); 
      TBranch *MyBranch_15 = tree->Branch("t50l_S1", &t50l_S1, 
					  "t50l_S1[2]/D"); 
      TBranch *MyBranch_16 = tree->Branch("t1_S1", &t1_S1, 
					  "t1_S1[2]/D"); 
      TBranch *MyBranch_17 = tree->Branch("t50r_S1", &t50r_S1, 
					  "t50r_S1[2]/D"); 
      TBranch *MyBranch_18 = tree->Branch("t10r_S1", &t10r_S1, 
					  "t10r_S1[2]/D");
      TBranch *MyBranch_19 = tree->Branch("t0r_S1", &t0r_S1, 
					  "t0r_S1[2]/D");
      TBranch *MyBranch_20 = tree->Branch("S1_mask", &S1_mask, 
					  "S1_mask[3]/O");
      TBranch *MyBranch_21 = tree->Branch("numS1", &numS1, 
					  "numS1[1]/D");
      TBranch *MyBranch_22 = tree->Branch("S1", &S1, 
					  "S1[2]/D");
      TBranch *MyBranch_24 = tree->Branch("BS_ave", &BS_ave, 
					  "BS_ave[2]/D");
      TBranch *MyBranch_25 = tree->Branch("BS_std", &BS_std, 
					  "BS_std[2]/D");
      TBranch *MyBranch_26 = tree->Branch("Sat", &Sat, 
					  "Sat[2]/O");
      TBranch *MyBranch_27 = tree->Branch("TAC_ave", &TAC_ave, 
					  "TAC_ave[1]/D");
      TBranch *MyBranch_28 = tree->Branch("NaI_Area", &NaI_Area, 
					  "NaI_Area[1]/D");

      for(int i=0;i<N_events;i++){ // loop over number of event
	if((i%1000)==0){cout << "Event = " << i <<endl;}
	int Variable_Length[1000000]={};

	for (int j=0; j<(VariableName.size());j++){// loop over the variables
	  for(int k=0; k<Dimm1[j];k++){ // loop over the first dimension
	    for(int l=0; l<Dimm2[j];l++){ // loop over the second dimension

	      unsigned int evtNum_v=0;
	      double pA_S1_v=-1;
	      double pH_S1_v=-1;
	      double t0l_S1_v=-1;
	      double t10l_S1_v=-1;
	      double t50l_S1_v=-1;
	      double t1_S1_v=-1;
	      double t50r_S1_v=-1;
	      double t10r_S1_v=-1;
	      double t0r_S1_v=-1;
	      bool S1_mask_v=-1;
	      double numS1_v=-1;
	      double S1_v=-1;
	      double BS_ave_v=-1;
	      double BS_std_v=-1;
	      bool Sat_v=-1;
	      double TAC_ave_v=-1;
	      double NaI_Area_v=-1;

	      Variable_Length[j]=(Length[j])*(Dimm1[j])*(Dimm2[j]);
	      Variable_Length[j]= Variable_Length[j]+Variable_Length[j-1];


	      BinaryFile.seekg(2 + H_Size 
			       + (k*Dimm2[j]*(Length[j]))
			       + (l*(Length[j]))
			       + (Variable_Length[j-1]) 
			       + (i*E_length), ios::beg);

	    
	      if(VariableName[j]=="evtNum"){
		BinaryFile.read((char*)(&evtNum_v),Length[j]);
		//cout << "----------" << endl
		//     << "Event Nummer = " << evtNum_v << endl;
		evtNum[l]=evtNum_v;
	      }

	      if(VariableName[j]=="pA_S1"){
	      	BinaryFile.read((char*)(&pA_S1_v),Length[j]);	       
		if(l==0){pA_S1_f[k]=pA_S1_v; //First S1
	   //cout <<"k = "<<k<<" "<<l<<" S1 area_f : "<<pA_S1_f[l]<< endl;
		}
		if(l==1){pA_S1_s[k]=pA_S1_v; //Second S1
	   //cout <<"k = "<<k<<" "<<l<<" S1 Pulse area_s : "<<pA_S1_s[l]<< endl
		}		
	      }

	      if(VariableName[j]=="pH_S1"){
	      	BinaryFile.read((char*)(&pH_S1_v),Length[j]);
		if(l==0){pH_S1_f[k]=pH_S1_v; //First S1
		}
		if(l==1){pH_S1_s[k]=pH_S1_v; //Second S1
		}
	      }

	      if(VariableName[j]=="t0l_S1"){
		BinaryFile.read((char*)(&t0l_S1_v),Length[j]);
	      	t0l_S1[l]=t0l_S1_v;
	      }

	      if(VariableName[j]=="t10l_S1"){
		BinaryFile.read((char*)(&t10l_S1_v),Length[j]);
	      	t10l_S1[l]=t10l_S1_v;
	      }

	      if(VariableName[j]=="t50l_S1"){
		BinaryFile.read((char*)(&t50l_S1_v),Length[j]);
	      	t50l_S1[l]=t50l_S1_v;
	      }

	      if(VariableName[j]=="t1_S1"){
		BinaryFile.read((char*)(&t1_S1_v),Length[j]);
	      	t1_S1[l]=t1_S1_v;
	      }

	      if(VariableName[j]=="t0r_S1"){
		BinaryFile.read((char*)(&t0r_S1_v),Length[j]);
	      	t0r_S1[l]=t0r_S1_v;
	      }

	      if(VariableName[j]=="t10r_S1"){
		BinaryFile.read((char*)(&t10r_S1_v),Length[j]);
	      	t10r_S1[l]=t10r_S1_v;
	      }

	      if(VariableName[j]=="t50r_S1"){
		BinaryFile.read((char*)(&t50r_S1_v),Length[j]);
	      	t50r_S1[l]=t50r_S1_v;
	      }

	      if(VariableName[j]=="S1_mask"){
	      	BinaryFile.read((char*)(&S1_mask_v),Length[j]);
		S1_mask[l]=S1_mask_v;
	      } 

	      if(VariableName[j]=="numS1"){
	      	BinaryFile.read((char*)(&numS1_v),Length[j]);
		numS1[l]=numS1_v;
	      } 
 
	      if(VariableName[j]=="S1"){
	      	BinaryFile.read((char*)(&S1_v),Length[j]);
		S1[l]=S1_v;
	      //cout <<" k = "<<k<<"; l = "<<l<<" ---- S1 : "<<S1[l]<< endl;
	      } 


	      if(VariableName[j]=="BS_ave"){
	      	BinaryFile.read((char*)(&BS_ave_v),Length[j]);
		BS_ave[l]=BS_ave_v;
	      } 

	      if(VariableName[j]=="BS_std"){
	      	BinaryFile.read((char*)(&BS_std_v),Length[j]);
		BS_std[l]=BS_std_v;
	      } 

	      if(VariableName[j]=="Sat"){
	      	BinaryFile.read((char*)(&Sat_v),Length[j]);
		Sat[l]=Sat_v;
	      } 

	      if(VariableName[j]=="TAC_ave"){
	      	BinaryFile.read((char*)(&TAC_ave_v),Length[j]);
		TAC_ave[l]=TAC_ave_v;
	      } 

	      if(VariableName[j]=="NaI_Area"){
	      	BinaryFile.read((char*)(&NaI_Area_v),Length[j]);
		NaI_Area[l]=NaI_Area_v;
	      } 
	     
	    }
	  }
	}// end loop over the variables

	tree->Fill(); 
       

      }// end loop over the events

 }//end if of S1 + TAC + NaI	

 else{

  // === === === === Here only S1 data (no S2, no TAC, NaI) === === === === 

 if (VariableName[0]=="evtNum" && 
     VariableName[1]=="pA_S1" && VariableName[2]=="pH_S1" && 
     VariableName[3]=="t0l_S1" && VariableName[4]=="t10l_S1" && 
     VariableName[5]=="t50l_S1" && VariableName[6]=="t1_S1" && 
     VariableName[7]=="t50r_S1" && VariableName[8]=="t10r_S1" && 
     VariableName[9]=="t0r_S1" && VariableName[10]=="S1_mask" && 
     VariableName[11]=="numS1" && VariableName[12]=="S1" && 
     VariableName[13]=="BS_ave" && VariableName[14]=="BS_std" && 
     VariableName[15]=="Sat"){
   
   cout  << "----------" << endl
	 << "XUERICH DATA FORMAT: 0-Field (without S2)" <<endl
	 << "----------" << endl;
   


      unsigned int evtNum[1]={};
      double pA_S1_f[2]={};
      double pA_S1_s[2]={};
      double pH_S1_f[2]={};
      double pH_S1_s[2]={};
      double t0l_S1[2]={};
      double t10l_S1[2]={};
      double t50l_S1[2]={};
      double t1_S1[2]={};
      double t50r_S1[2]={};
      double t10r_S1[2]={};
      double t0r_S1[2]={};
      bool S1_mask[2]={};
      double numS1[1]={};
      double S1[2]={};
      double BS_ave[2]={};
      double BS_std[2]={};
      bool Sat[2]={};



      TBranch *MyBranch_1 = tree->Branch("evtNum", &evtNum, 
					 "evtNum[1]/i");
      TBranch *MyBranch_11_f = tree->Branch("pA_S1_f", &pA_S1_f, 
					"pA_S1_f[2]/D"); 
      TBranch *MyBranch_11_s = tree->Branch("pA_S1_s", &pA_S1_s, 
					"pA_S1_s[2]/D"); 
      TBranch *MyBranch_12_f = tree->Branch("pH_S1_f", &pH_S1_f, 
					"pH_S1_f[2]/D"); 
      TBranch *MyBranch_12_s = tree->Branch("pH_S1_s", &pH_S1_s, 
					"pH_S1_s[2]/D"); 
      TBranch *MyBranch_13 = tree->Branch("t0l_S1", &t0l_S1, 
					"t0l_S1[2]/D"); 
      TBranch *MyBranch_14 = tree->Branch("t10l_S1", &t10l_S1, 
					"t10l_S1[2]/D"); 
      TBranch *MyBranch_15 = tree->Branch("t50l_S1", &t50l_S1, 
					  "t50l_S1[2]/D"); 
      TBranch *MyBranch_16 = tree->Branch("t1_S1", &t1_S1, 
					  "t1_S1[2]/D"); 
      TBranch *MyBranch_17 = tree->Branch("t50r_S1", &t50r_S1, 
					  "t50r_S1[2]/D"); 
      TBranch *MyBranch_18 = tree->Branch("t10r_S1", &t10r_S1, 
					  "t10r_S1[2]/D");
      TBranch *MyBranch_19 = tree->Branch("t0r_S1", &t0r_S1, 
					  "t0r_S1[2]/D");
      TBranch *MyBranch_20 = tree->Branch("S1_mask", &S1_mask, 
					  "S1_mask[3]/O");
      TBranch *MyBranch_21 = tree->Branch("numS1", &numS1, 
					  "numS1[1]/D");
      TBranch *MyBranch_22 = tree->Branch("S1", &S1, 
					  "S1[2]/D");
      TBranch *MyBranch_24 = tree->Branch("BS_ave", &BS_ave, 
					  "BS_ave[2]/D");
      TBranch *MyBranch_25 = tree->Branch("BS_std", &BS_std, 
					  "BS_std[2]/D");
      TBranch *MyBranch_26 = tree->Branch("Sat", &Sat, 
					  "Sat[2]/O");


      for(int i=0;i<N_events;i++){ // loop over number of event
	if((i%1000)==0){cout << "Event = " << i <<endl;}
	int Variable_Length[1000000]={};

	for (int j=0; j<(VariableName.size());j++){// loop over the variables
	  for(int k=0; k<Dimm1[j];k++){ // loop over the first dimension
	    for(int l=0; l<Dimm2[j];l++){ // loop over the second dimension

	      unsigned int evtNum_v=0;
	      double pA_S1_v=-1;
	      double pH_S1_v=-1;
	      double t0l_S1_v=-1;
	      double t10l_S1_v=-1;
	      double t50l_S1_v=-1;
	      double t1_S1_v=-1;
	      double t50r_S1_v=-1;
	      double t10r_S1_v=-1;
	      double t0r_S1_v=-1;
	      bool S1_mask_v=-1;
	      double numS1_v=-1;
	      double S1_v=-1;
	      double BS_ave_v=-1;
	      double BS_std_v=-1;
	      bool Sat_v=-1;


	      Variable_Length[j]=(Length[j])*(Dimm1[j])*(Dimm2[j]);
	      Variable_Length[j]= Variable_Length[j]+Variable_Length[j-1];


	      BinaryFile.seekg(2 + H_Size 
			       + (k*Dimm2[j]*(Length[j]))
			       + (l*(Length[j]))
			       + (Variable_Length[j-1]) 
			       + (i*E_length), ios::beg);

	    
	      if(VariableName[j]=="evtNum"){
		BinaryFile.read((char*)(&evtNum_v),Length[j]);
		//cout << "----------" << endl
		//     << "Event Nummer = " << evtNum_v << endl;
		evtNum[l]=evtNum_v;
	      }

	      if(VariableName[j]=="pA_S1"){
	      	BinaryFile.read((char*)(&pA_S1_v),Length[j]);	       
		if(l==0){pA_S1_f[k]=pA_S1_v; //First S1
	   //cout <<"k = "<<k<<" "<<l<<" S1 area_f : "<<pA_S1_f[l]<< endl;
		}
		if(l==1){pA_S1_s[k]=pA_S1_v; //Second S1
	   //cout <<"k = "<<k<<" "<<l<<" S1 Pulse area_s : "<<pA_S1_s[l]<< endl
		}		
	      }

	      if(VariableName[j]=="pH_S1"){
	      	BinaryFile.read((char*)(&pH_S1_v),Length[j]);
		if(l==0){pH_S1_f[k]=pH_S1_v; //First S1
		}
		if(l==1){pH_S1_s[k]=pH_S1_v; //Second S1
		}
	      }

	      if(VariableName[j]=="t0l_S1"){
		BinaryFile.read((char*)(&t0l_S1_v),Length[j]);
	      	t0l_S1[l]=t0l_S1_v;
	      }

	      if(VariableName[j]=="t10l_S1"){
		BinaryFile.read((char*)(&t10l_S1_v),Length[j]);
	      	t10l_S1[l]=t10l_S1_v;
	      }

	      if(VariableName[j]=="t50l_S1"){
		BinaryFile.read((char*)(&t50l_S1_v),Length[j]);
	      	t50l_S1[l]=t50l_S1_v;
	      }

	      if(VariableName[j]=="t1_S1"){
		BinaryFile.read((char*)(&t1_S1_v),Length[j]);
	      	t1_S1[l]=t1_S1_v;
	      }

	      if(VariableName[j]=="t0r_S1"){
		BinaryFile.read((char*)(&t0r_S1_v),Length[j]);
	      	t0r_S1[l]=t0r_S1_v;
	      }

	      if(VariableName[j]=="t10r_S1"){
		BinaryFile.read((char*)(&t10r_S1_v),Length[j]);
	      	t10r_S1[l]=t10r_S1_v;
	      }

	      if(VariableName[j]=="t50r_S1"){
		BinaryFile.read((char*)(&t50r_S1_v),Length[j]);
	      	t50r_S1[l]=t50r_S1_v;
	      }

	      if(VariableName[j]=="S1_mask"){
	      	BinaryFile.read((char*)(&S1_mask_v),Length[j]);
		S1_mask[l]=S1_mask_v;
	      } 

	      if(VariableName[j]=="numS1"){
	      	BinaryFile.read((char*)(&numS1_v),Length[j]);
		numS1[l]=numS1_v;
	      } 
 
	      if(VariableName[j]=="S1"){
	      	BinaryFile.read((char*)(&S1_v),Length[j]);
		S1[l]=S1_v;
	      //cout <<" k = "<<k<<"; l = "<<l<<" ---- S1 : "<<S1[l]<< endl;
	      } 


	      if(VariableName[j]=="BS_ave"){
	      	BinaryFile.read((char*)(&BS_ave_v),Length[j]);
		BS_ave[l]=BS_ave_v;
	      } 

	      if(VariableName[j]=="BS_std"){
	      	BinaryFile.read((char*)(&BS_std_v),Length[j]);
		BS_std[l]=BS_std_v;
	      } 

	      if(VariableName[j]=="Sat"){
	      	BinaryFile.read((char*)(&Sat_v),Length[j]);
		Sat[l]=Sat_v;
	      } 
	     
	    }
	  }
	}// end loop over the variables

	tree->Fill(); 
       

      }// end loop over the events
 }//end of if (only S1)		  

 }//end of second else

  }// end of first else


  hfile->Write(); 
  hfile->Close();
  
  
  // Close the file --------- ----------
  BinaryFile.close();
  
 }

