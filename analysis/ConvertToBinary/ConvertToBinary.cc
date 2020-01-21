#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TMath.h"
#include "TChain.h"
#include "TStyle.h"

#include "vector"
#include "string"
#include "iostream"
#include "fstream"
#include "stdio.h"
#include "stdlib.h"


#include "genStruct.hh"

int main()
{	

	// load data
	//const char *DataFolder 	= "/Users/alexkish/data/xuerich/mc/compton";	// idark05x
	const char *DataFolder 	= "/Users/alexkish/data/xuerich/mc/Co57source";	// idark05x
		
	//const char *inFileName = "Xu_Compton_3.5deg_4e10-t2";
	//const char *inFileName = "Xu_Compton_4.25deg_4e10-t2";
	//const char *inFileName = "Xu_Compton_4.25deg1_4e10-t2";
	//const char *inFileName = "Xu_Compton_6.25deg_4e10-t2";
	//const char *inFileName = "Xu_Compton_8.5deg_4e10-t2";
	//const char *inFileName = "Xu_Compton_16.25deg_4e10-t2";
	//const char *inFileName = "Xu_Compton_34.5deg_4e10-t2";
	//const char *inFileName = "Xu_Compton_34.5deg1_4e10-t2";

	//const char *inFileName = "Xu_Compton_3.5deg-70cm_4e10-t2";
	//const char *inFileName = "Xu_Compton_4.25deg-70cm_4e10-t2";
	//const char *inFileName = "Xu_Compton_6.25deg-70cm_4e10-t2";
	//const char *inFileName = "Xu_Compton_8.5deg-70cm_4e10-t2";
	//const char *inFileName = "Xu_Compton_16.25deg-70cm_4e10-t2";
	//const char *inFileName = "Xu_Compton_34.5deg-70cm_4e10-t2";

	//const char *inFileName = "Xu_Compton_6.125deg-70cm-Vacuum_4e10-t2";
	//const char *inFileName = "Xu_Compton_6.25deg-70cm-Vacuum_4e10-t2";
	//const char *inFileName = "Xu_Compton_6.375deg-70cm-Vacuum_4e10-t2";

	//const char *inFileName = "Xu_Compton_8.375deg-70cm-Vacuum_4e10-t2";
	//const char *inFileName = "Xu_Compton_8.5deg-70cm-Vacuum_4e10-t2";
	//const char *inFileName = "Xu_Compton_8.625deg-70cm-Vacuum_4e10-t2";

	//const char *inFileName = "Xu_Compton_4.125deg-28cm-Vacuum_4e11-t2";
	//const char *inFileName = "Xu_Compton_4.25deg-28cm-Vacuum_4e11-t2";
	//const char *inFileName = "Xu_Compton_4.375deg-28cm-Vacuum_4e11-t2";
	//const char *inFileName = "Xu_Compton_5.25deg-28cm-Vacuum_4e11-t2";
	//const char *inFileName = "Xu_Compton_6.125deg-70cm-Vacuum_4e11-t2";
	//const char *inFileName = "Xu_Compton_6.25deg-70cm-Vacuum_4e11-t2";
	//const char *inFileName = "Xu_Compton_6.375deg-70cm-Vacuum_4e11-t2";
	//const char *inFileName = "Xu_Compton_8.5deg-70cm-Vacuum_4e11-t2";
	//const char *inFileName = "Xu_Compton_8.375deg-28cm-Vacuum_3.952e11-t2";
	//const char *inFileName = "Xu_Compton_8.5deg-28cm-Vacuum_4e11-t2";
	//const char *inFileName = "Xu_Compton_16.25deg-70cm-Vacuum_4e11-t2";
	const char *inFileName = "Xu_Compton_34.5deg-70cm-Vacuum_4e11-t2";
	
	//const char *inFileName = "Xu_Co57sourceRight_2e10-t2";
	//const char *inFileName = "XuV3-0deg_Co57source-0mm_7.88e10-t2";

	//const char *inFileName = "Xu0_Co57source-0mm_2e10-t2";
	//const char *inFileName = "Xu0_Co57source-Minus10mm_2e10-t2";
	//const char *inFileName = "Xu0_Co57source-Minus20mm_1e10-t2";
	//const char *inFileName = "Xu0_Co57source-Plus10mm_1e10-t2";
	
	
	Char_t  inFile[500];
	sprintf(inFile, "%s/%s.root", DataFolder, inFileName);

	// load tree
	TChain *t2	= new TChain("t2");
			t2	->Add(inFile);

 	//--------------------------------------------------------------------------------------å
	// GET NUMBER OF ENTRIES IN TTREE
	int k = t2->GetEntries();
	cout <<"File "<< inFileName <<" has "<< k <<" entries."<< endl;
	//--------------------------------------------------------------------------------------
	
	
	// define variables in the original tree
	int hit[8];
	float hit_Edep[8];
	int hitBefore[8];
	float hitBefore_Edep[8];
	int hitAfter[8];
	float hitAfter_Edep[8];

	int LXe_ns;
	float LXe_Edep;
	float LXe_tm[20];
	float LXe_EkinPre[20];
	float LXe_EkinPost[20];
	float LXe_x[20];
	float LXe_y[20];
	float LXe_z[20];

	float NaI_Edep;
	float NaI_tm;
	float NaI_x;
	float NaI_y;
	float NaI_z;

	float LXe_tmMin;
	float LXe_tmMax;
	float LXe_EkinMax;
	float LXe_EkinMin;
	float LXe_EkinDelta[20];
	float LXe_EkinDeltaTot;

	float NaI_tmMin;
	float NaI_tmMax;
	
	float primary_x;
	float primary_y;
	float primary_z;

	// define tree structure and pointers
	t2->SetBranchAddress("LXe_ns", 				&LXe_ns);
	t2->SetBranchAddress("LXe_Edep", 			&LXe_Edep);
	t2->SetBranchAddress("LXe_EkinPre", 		&LXe_EkinPre);
	t2->SetBranchAddress("LXe_EkinPost", 		&LXe_EkinPost);
	t2->SetBranchAddress("LXe_tm", 				&LXe_tm);
	t2->SetBranchAddress("LXe_x", 				&LXe_x);
	t2->SetBranchAddress("LXe_y", 				&LXe_y);
	t2->SetBranchAddress("LXe_z", 				&LXe_z);

	t2->SetBranchAddress("NaI_Edep", 			&NaI_Edep);
	t2->SetBranchAddress("NaI_tm", 				&NaI_tm);
	t2->SetBranchAddress("NaI_x", 				&NaI_x);
	t2->SetBranchAddress("NaI_y", 				&NaI_y);
	t2->SetBranchAddress("NaI_z", 				&NaI_z);

	t2->SetBranchAddress("hit", 				&hit);
	t2->SetBranchAddress("hit_Edep", 			&hit_Edep);

	// additional  branches
	t2->SetBranchAddress("hitBefore", 			&hitBefore);
	t2->SetBranchAddress("hitBefore_Edep",		&hitBefore);
	t2->SetBranchAddress("hitAfter", 			&hitAfter);
	t2->SetBranchAddress("hitAfter_Edep", 		&hitAfter);

	t2->SetBranchAddress("primary_x", 			&primary_x);
	t2->SetBranchAddress("primary_y", 			&primary_y);
	t2->SetBranchAddress("primary_z", 			&primary_z);

	t2->SetBranchAddress("LXe_tmMin", 			&LXe_tmMin);
	t2->SetBranchAddress("LXe_tmMax", 			&LXe_tmMax);
	t2->SetBranchAddress("NaI_tmMin", 			&NaI_tmMin);
	t2->SetBranchAddress("NaI_tmMax", 			&NaI_tmMax);
	t2->SetBranchAddress("LXe_EkinMax", 		&LXe_EkinMax);
	t2->SetBranchAddress("LXe_EkinMin", 		&LXe_EkinMin);
	t2->SetBranchAddress("LXe_EkinDelta", 		&LXe_EkinDelta);
	t2->SetBranchAddress("LXe_EkinDeltaTot", 	&LXe_EkinDeltaTot);


	// OUTPUT BINARY FILE
	Char_t  outFileName[500];
	sprintf(outFileName, "%s/%s.dat", DataFolder, inFileName);

	genStruct *outFile = new genStruct(outFileName);

	int *b_LXe_ns = new int;
	double *b_LXe_Edep = new double;
	double *b_LXe_EkinMin = new double;
	double *b_LXe_EkinMax = new double;
	double *b_LXe_EkinDeltaTot = new double;

	double *b_NaI_Edep = new double;
	double *b_NaI_tm = new double;
	double *b_NaI_x = new double;
	double *b_NaI_y = new double;
	double *b_NaI_z = new double;

	double *b_primary_x = new double;
	double *b_primary_y = new double;
	double *b_primary_z = new double;

	vector<double> *b_LXe_EkinPre = new vector<double>;
	vector<double> *b_LXe_EkinPost = new vector<double>;
	vector<double> *b_LXe_EkinDelta = new vector<double>;
	vector<double> *b_LXe_tm = new vector<double>;
	vector<double> *b_LXe_x = new vector<double>;
	vector<double> *b_LXe_y = new vector<double>;
	vector<double> *b_LXe_z = new vector<double>;

	vector<double> *b_hit = new vector<double>;
	vector<double> *b_hit_Edep = new vector<double>;
	vector<double> *b_hitBefore = new vector<double>;
	vector<double> *b_hitBefore_Edep = new vector<double>;
	vector<double> *b_hitAfter = new vector<double>;
	vector<double> *b_hitAfter_Edep = new vector<double>;

	outFile->AddField(b_LXe_ns,"b_LXe_ns","1,1");
	outFile->AddField(b_LXe_Edep,"b_LXe_Edep","1,1");
	outFile->AddField(b_LXe_EkinPre,"b_LXe_EkinPre","1,1");
	outFile->AddField(b_LXe_EkinPost,"b_LXe_EkinPost","1,1");
	outFile->AddField(b_LXe_EkinMin,"b_LXe_EkinMin","1,1");
	outFile->AddField(b_LXe_EkinMax,"b_LXe_EkinMax","1,1");
	outFile->AddField(b_LXe_EkinDelta,"b_LXe_EkinDelta","1,1");
	outFile->AddField(b_LXe_EkinDeltaTot,"b_LXe_EkinDeltaTot","1,1");
	outFile->AddField(b_LXe_tm,"b_LXe_tm","1,1");
	outFile->AddField(b_LXe_x,"b_LXe_x","1,1");
	outFile->AddField(b_LXe_y,"b_LXe_y","1,1");
	outFile->AddField(b_LXe_z,"b_LXe_z","1,1");

	outFile->AddField(b_NaI_Edep,"b_NaI_Edep","1,1");
	outFile->AddField(b_NaI_tm,"b_NaI_tm","1,1");
	outFile->AddField(b_NaI_x,"b_NaI_x","1,1");
	outFile->AddField(b_NaI_y,"b_NaI_y","1,1");
	outFile->AddField(b_NaI_z,"b_NaI_z","1,1");

	outFile->AddField(b_primary_x,"b_primary_x","1,1");
	outFile->AddField(b_primary_y,"b_primary_y","1,1");
	outFile->AddField(b_primary_z,"b_primary_z","1,1");

	outFile->AddField(b_hit,"b_hit","1,1");
	outFile->AddField(b_hit_Edep,"b_hit_Edep","1,1");
	outFile->AddField(b_hitBefore,"b_hitBefore","1,1");
	outFile->AddField(b_hitBefore_Edep,"b_hitBefore_Edep","1,1");
	outFile->AddField(b_hitAfter,"b_hitAfter","1,1");
	outFile->AddField(b_hitAfter_Edep,"b_hitAfter_Edep","1,1");

	outFile->MakeHeader();
	
	outFile->Write();


	// PROCESS DATA
	for(int i=0; i<k; i++){ // loop on events
	//for(int i=0; i<10; i++){ // loop on events
		t2->GetEntry(i);

		//cout <<"Event = "<< i <<";  LXe_ns = "<< LXe_ns << endl;
		//cout <<"i = "<< i <<"; LXe_Edep = " << LXe_Edep << endl;
		// second iteration
		*b_primary_x = primary_x;
		*b_primary_y = primary_y;
		*b_primary_z = primary_z;

		*b_LXe_ns = LXe_ns;
		*b_LXe_Edep = LXe_Edep;
		*b_LXe_EkinMin = LXe_EkinMin;
		*b_LXe_EkinMax = LXe_EkinMax;
		*b_LXe_EkinDeltaTot = LXe_EkinDeltaTot;

		*b_NaI_Edep = NaI_Edep;
		*b_NaI_tm = NaI_tm;
		*b_NaI_x = NaI_x;
		*b_NaI_y = NaI_y;
		*b_NaI_z = NaI_z;

		b_LXe_EkinPre->clear();
		b_LXe_EkinPost->clear();
		b_LXe_EkinDelta->clear();
		b_LXe_tm->clear();
		b_LXe_x->clear();
		b_LXe_y->clear();
		b_LXe_z->clear();
		
		for(int j=0; j<LXe_ns; j++)	{ // loop on scatters in each event
			b_LXe_EkinPre->push_back(LXe_EkinPre[j]);
			b_LXe_EkinPost->push_back(LXe_EkinPost[j]);
			b_LXe_EkinDelta->push_back(LXe_EkinDelta[j]);
			b_LXe_tm->push_back(LXe_tm[j]);
			b_LXe_x->push_back(LXe_x[j]);
			b_LXe_y->push_back(LXe_y[j]);
			b_LXe_z->push_back(LXe_z[j]);
		}

		b_hit->clear();
		b_hit_Edep->clear();
		b_hitBefore->clear();
		b_hitBefore_Edep->clear();
		b_hitAfter->clear();
		b_hitAfter_Edep->clear();
		for(int m=0; m<8; m++)	{
			b_hit->push_back(hit[m]);
			b_hit_Edep->push_back(hit_Edep[m]);
			b_hitBefore->push_back(hitBefore[m]);
			b_hitBefore_Edep->push_back(hitBefore_Edep[m]);
			b_hitAfter->push_back(hitAfter[m]);
			b_hitAfter_Edep->push_back(hitAfter_Edep[m]);		
		}
		outFile->Write();
	} // end loop on events


	delete outFile;
	
	delete b_LXe_ns;
	delete b_LXe_Edep;
	delete b_LXe_EkinMin;
	delete b_LXe_EkinMax;
	delete b_LXe_EkinDeltaTot;

	delete b_NaI_Edep;
	delete b_NaI_tm;
	delete b_NaI_x;
	delete b_NaI_y;
	delete b_NaI_z;

	delete b_primary_x;
	delete b_primary_y;
	delete b_primary_z;

	delete b_LXe_EkinPre;
	delete b_LXe_EkinPost;
	delete b_LXe_EkinDelta;
	delete b_LXe_tm;
	delete b_LXe_x;
	delete b_LXe_y;
	delete b_LXe_z;

	delete b_hit;
	delete b_hit_Edep;
	delete b_hitBefore;
	delete b_hitBefore_Edep;
	delete b_hitAfter;
	delete b_hitAfter_Edep;
	
	
	
	return 0;

}
