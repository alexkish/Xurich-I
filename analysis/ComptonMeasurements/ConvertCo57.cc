#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TMath.h"
#include "TRandom.h"
#include "TChain.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCut.h"
#include "TPie.h"
#include "TPieSlice.h"
#include "TPaletteAxis.h"
#include "TList.h"
#include "TColor.h"
#include "TTreePlayer.h"
#include "TEventList.h"

#include "vector"
#include "string"
#include "iostream"
#include "fstream"
#include "stdio.h"
#include "stdlib.h"


void ConvertCo57()
{	
	gStyle	->SetOptStat(0);
	gStyle	->SetOptFit(0);
	
	gStyle->SetStatBorderSize(0);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFillColor(10);
	gStyle->SetStatColor(10);
	gStyle->SetStatFont(42);
	//gStyle->SetMarkerColor(4);
	gStyle->SetMarkerStyle(8);
	gStyle->SetMarkerSize(0.6);
	
	// define color gradinet
	const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(30);

	// load data
	const char *DataFolder 	= "/Users/alexkish/data/xuerich/mc/Co57source";	// idark05x
		
	//const char *inFileName = "Xu_Compton_3.5deg_4e10";
	//const char *inFileName = "Xu_Compton_4.25deg_4e10";
	//const char *inFileName = "Xu_Compton_4.25deg1_4e10";
	//const char *inFileName = "Xu_Compton_6.25deg_4e10";
	//const char *inFileName = "Xu_Compton_6.25deg1_4e10";
	//const char *inFileName = "Xu_Compton_8.5deg_4e10";
	//const char *inFileName = "Xu_Compton_16.25deg_4e10";
	//const char *inFileName = "Xu_Compton_34.5deg_4e10";
	//const char *inFileName = "Xu_Compton_34.5deg1_4e10";

	//const char *inFileName = "Xu_Compton_3.5deg-70cm_4e10";
	//const char *inFileName = "Xu_Compton_4.25deg-70cm_4e10";
	//const char *inFileName = "Xu_Compton_6.25deg-70cm_4e10";
	//const char *inFileName = "Xu_Compton_8.5deg-70cm_4e10";
	//const char *inFileName = "Xu_Compton_16.25deg-70cm_4e10";
	//const char *inFileName = "Xu_Compton_34.5deg-70cm_4e10";
	//const char *inFileName = "Xu_Compton_4.25deg-Beam662keV_2e10";

	//const char *inFileName = "Xu_Co57sourceRight_2e10";

	//const char *inFileName = "Xu0_Co57source-0mm_2e10";
	//const char *inFileName = "Xu0_Co57source-0mm_1.16e11";
	//const char *inFileName = "Xu0_Co57source-Minus10mm_2e10";
	//const char *inFileName = "Xu0_Co57source-Minus20mm_1e10";
	//const char *inFileName = "Xu0_Co57source-Plus10mm_1e10";
	const char *inFileName = "XuV3-0deg_Co57source-0mm_7.88e10";

	Char_t  inFile[500];
	sprintf(inFile, "%s/%s.root", DataFolder, inFileName);

	TChain *t1	= new TChain("t1");
			t1	->Add(inFile);

	int eventid;
	int nsteps;
	float etot;
	float xp_pri;
	float yp_pri;
	float zp_pri;
	std::vector<int> *parentid 		= new vector<int>;
	std::vector<double> *pre_kinE 	= new vector<double>;
	std::vector<double> *post_kinE 	= new vector<double>;
	std::vector<float> *time 		= new vector<float>;
	std::vector<float> *ed 			= new vector<float>;
	std::vector<float> *xp 			= new vector<float>;
	std::vector<float> *yp 			= new vector<float>;
	std::vector<float> *zp 			= new vector<float>;

	int NaI_nsteps;
	float NaI_etot;
	std::vector<float> *NaI_time 	= new vector<float>;
	std::vector<float> *NaI_ed 		= new vector<float>;
	std::vector<float> *NaI_xp 		= new vector<float>;
	std::vector<float> *NaI_yp 		= new vector<float>;
	std::vector<float> *NaI_zp 		= new vector<float>;


	int DeadXe_nsteps;
	int GXe_nsteps;
	int Teflon_nsteps;
	int Steel_nsteps;
	int Lead_nsteps;
	int Al_nsteps;
	int Rest_nsteps;

	float DeadXe_etot;
	float GXe_etot;
	float Teflon_etot;
	float Steel_etot;
	float Lead_etot;
	float Al_etot;
	float Rest_etot;

	std::vector<float> *DeadXe_time = new vector<float>;
	std::vector<float> *GXe_time 	= new vector<float>;
	std::vector<float> *Teflon_time = new vector<float>;
	std::vector<float> *Steel_time 	= new vector<float>;
	std::vector<float> *Lead_time 	= new vector<float>;
	std::vector<float> *Al_time 	= new vector<float>;
	std::vector<float> *Rest_time 	= new vector<float>;

	std::vector<float> *DeadXe_ed 	= new vector<float>;
	std::vector<float> *GXe_ed 		= new vector<float>;
	std::vector<float> *Teflon_ed 	= new vector<float>;
	std::vector<float> *Steel_ed 	= new vector<float>;
	std::vector<float> *Lead_ed 	= new vector<float>;
	std::vector<float> *Al_ed 		= new vector<float>;
	std::vector<float> *Rest_ed 	= new vector<float>;

	t1->SetBranchStatus("*",0); // 1st - disable all branches to speed up processing
	t1->SetBranchAddress("eventid",			&eventid);	
	t1->SetBranchAddress("nsteps",			&nsteps);
	t1->SetBranchAddress("parentid",		&parentid);
	t1->SetBranchAddress("etot",			&etot);
	t1->SetBranchAddress("pre_kinE",		&pre_kinE);
	t1->SetBranchAddress("post_kinE",		&post_kinE);
	t1->SetBranchAddress("time",			&time);
	t1->SetBranchAddress("ed",				&ed);
	t1->SetBranchAddress("xp",				&xp);
	t1->SetBranchAddress("yp",				&yp);
	t1->SetBranchAddress("zp",				&zp);
	t1->SetBranchAddress("xp_pri",			&xp_pri);
	t1->SetBranchAddress("yp_pri",			&yp_pri);
	t1->SetBranchAddress("zp_pri",			&zp_pri);
	
	t1->SetBranchAddress("NaI_nsteps",		&NaI_nsteps);
	t1->SetBranchAddress("NaI_etot",		&NaI_etot);
	t1->SetBranchAddress("NaI_time",		&NaI_time);
	t1->SetBranchAddress("NaI_ed",			&NaI_ed);
	t1->SetBranchAddress("NaI_xp",			&NaI_xp);
	t1->SetBranchAddress("NaI_yp",			&NaI_yp);
	t1->SetBranchAddress("NaI_zp",			&NaI_zp);

	t1->SetBranchAddress("DeadXe_nsteps",	&DeadXe_nsteps);
	t1->SetBranchAddress("DeadXe_etot",		&DeadXe_etot);
	t1->SetBranchAddress("DeadXe_time",		&DeadXe_time);
	t1->SetBranchAddress("DeadXe_ed",		&DeadXe_ed);
	
	t1->SetBranchAddress("GXe_nsteps",		&GXe_nsteps);
	t1->SetBranchAddress("GXe_etot",		&GXe_etot);
	t1->SetBranchAddress("GXe_time",		&GXe_time);
	t1->SetBranchAddress("GXe_ed",			&GXe_ed);

	t1->SetBranchAddress("Teflon_nsteps",	&Teflon_nsteps);
	t1->SetBranchAddress("Teflon_etot",		&Teflon_etot);
	t1->SetBranchAddress("Teflon_time",		&Teflon_time);
	t1->SetBranchAddress("Teflon_ed",		&Teflon_ed);

	t1->SetBranchAddress("Steel_nsteps",	&Steel_nsteps);
	t1->SetBranchAddress("Steel_etot",		&Steel_etot);
	t1->SetBranchAddress("Steel_time",		&Steel_time);
	t1->SetBranchAddress("Steel_ed",		&Steel_ed);

	t1->SetBranchAddress("Lead_nsteps",		&Lead_nsteps);
	t1->SetBranchAddress("Lead_etot",		&Lead_etot);
	t1->SetBranchAddress("Lead_time",		&Lead_time);
	t1->SetBranchAddress("Lead_ed",			&Lead_ed);

	t1->SetBranchAddress("Al_nsteps",		&Al_nsteps);
	t1->SetBranchAddress("Al_etot",			&Al_etot);
	t1->SetBranchAddress("Al_time",			&Al_time);
	t1->SetBranchAddress("Al_ed",			&Al_ed);

	t1->SetBranchAddress("Rest_nsteps",		&Rest_nsteps);
	t1->SetBranchAddress("Rest_etot",		&Rest_etot);
	t1->SetBranchAddress("Rest_time",		&Rest_time);
	t1->SetBranchAddress("Rest_ed",			&Rest_ed);
 	
	//--------------------------------------------------------------------------------------
	int k = t1->GetEntries();
	cout <<"File "<< inFileName <<" has "<< k <<" entries."<< endl;
	//--------------------------------------------------------------------------------------
	
	// OUTPUT FILE
	int hit[8];
	float hit_Edep[8];
	int hitBefore[8];
	float hitBefore_Edep[8];
	int hitAfter[8];
	float hitAfter_Edep[8];
	int hitBetween[8];
	float hitBetween_Edep[8];

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

	float LXe_tmMin = 0;
	float LXe_tmMax = 0;
	float LXe_EkinMax;
	float LXe_EkinMin;
	float LXe_EkinDelta[20];
	float LXe_EkinDeltaTot;

	int NaI_ns;
	float NaI_tmMin = 0;
	float NaI_tmMax = 0;
	
	float primary_x = 0;
	float primary_y = 0;
	float primary_z = 0;

	int goodEvent;
	int eventsRemoved = 0;

	Char_t  outFileName[500];
	sprintf(outFileName, "%s/%s-t2.root", DataFolder, inFileName);

	TFile *outFile = new TFile(outFileName,"RECREATE");
	TTree *t2 = new TTree("t2","t2");

	t2->Branch("LXe_ns", 			&LXe_ns,			"LXe_ns/I");
	t2->Branch("LXe_Edep", 			&LXe_Edep,			"LXe_Edep/F");
	t2->Branch("LXe_EkinPre", 		&LXe_EkinPre,		"LXe_EkinPre[LXe_ns]/F");
	t2->Branch("LXe_EkinPost", 		&LXe_EkinPost,		"LXe_EkinPost[LXe_ns]/F");
	t2->Branch("LXe_tm", 			&LXe_tm,			"LXe_tm[LXe_ns]/F");
	t2->Branch("LXe_x", 			&LXe_x,				"LXe_x[LXe_ns]/F");
	t2->Branch("LXe_y", 			&LXe_y,				"LXe_y[LXe_ns]/F");
	t2->Branch("LXe_z", 			&LXe_z,				"LXe_z[LXe_ns]/F");

	t2->Branch("NaI_Edep", 			&NaI_Edep,			"NaI_Edep/F");
	t2->Branch("NaI_tm", 			&NaI_tm,			"NaI_tm/F");
	t2->Branch("NaI_x", 			&NaI_x,				"NaI_x/F");
	t2->Branch("NaI_y", 			&NaI_y,				"NaI_y/F");
	t2->Branch("NaI_z", 			&NaI_z,				"NaI_z/F");

	t2->Branch("hit", 				&hit,				"hit[8]/I");
	t2->Branch("hit_Edep", 			&hit_Edep,			"hit_Edep[8]/F");

	// additional  branches
	t2->Branch("hitBefore", 		&hitBefore,			"hitBefore[8]/I");
	t2->Branch("hitBefore_Edep",	&hitBefore_Edep,	"hitBefore_Edep[8]/F");
	t2->Branch("hitAfter", 			&hitAfter,			"hitAfter[8]/I");
	t2->Branch("hitAfter_Edep", 	&hitAfter_Edep,		"hitAfter_Edep[8]/F");
	t2->Branch("hitBetween", 		&hitBetween,		"hitBetween[8]/I");
	t2->Branch("hitBetween_Edep", 	&hitBetween_Edep,	"hitBetween_Edep[8]/F");

	t2->Branch("primary_x", 		&primary_x,			"primary_x/F");
	t2->Branch("primary_y", 		&primary_y,			"primary_y/F");
	t2->Branch("primary_z", 		&primary_z,			"primary_z/F");

	t2->Branch("LXe_tmMin", 		&LXe_tmMin,			"LXe_tmMin/F");
	t2->Branch("LXe_tmMax", 		&LXe_tmMax,			"LXe_tmMax/F");
	t2->Branch("NaI_tmMin", 		&NaI_tmMin,			"NaI_tmMin/F");
	t2->Branch("NaI_tmMax", 		&NaI_tmMax,			"NaI_tmMax/F");
	t2->Branch("LXe_EkinMax", 		&LXe_EkinMax,		"LXe_EkinMax/F");
	t2->Branch("LXe_EkinMin", 		&LXe_EkinMin,		"LXe_EkinMin/F");
	t2->Branch("LXe_EkinDelta", 	&LXe_EkinDelta,		"LXe_EkinDelta[LXe_ns]/F");
	t2->Branch("LXe_EkinDeltaTot", 	&LXe_EkinDeltaTot,	"LXe_EkinDeltaTot/F");

	//t2->Branch("goodEvent", 		&goodEvent,			"goodEvent/I");
	
	for(int i=0; i<k; i++){
	/*	if(i%100==0){
			cout << "\rProcessing event " << i << " of " << k;
			cout.flush();
		}	
	*/	t1->GetEntry(i);
		
		goodEvent = 0;

		hit_Edep[0] = 0;	hitBefore_Edep[0] = 0;	hitAfter_Edep[0] = 0;	hitBetween_Edep[0] = 0;
		hit_Edep[1] = 0;	hitBefore_Edep[1] = 0;	hitAfter_Edep[1] = 0;	hitBetween_Edep[1] = 0;
		hit_Edep[2] = 0;	hitBefore_Edep[2] = 0;	hitAfter_Edep[2] = 0;	hitBetween_Edep[2] = 0;
		hit_Edep[3] = 0;	hitBefore_Edep[3] = 0;	hitAfter_Edep[3] = 0;	hitBetween_Edep[3] = 0;
		hit_Edep[4] = 0;	hitBefore_Edep[4] = 0;	hitAfter_Edep[4] = 0;	hitBetween_Edep[4] = 0;
		hit_Edep[5] = 0;	hitBefore_Edep[5] = 0;	hitAfter_Edep[5] = 0;	hitBetween_Edep[5] = 0;
		hit_Edep[6] = 0;	hitBefore_Edep[6] = 0;	hitAfter_Edep[6] = 0;	hitBetween_Edep[6] = 0;
		hit_Edep[7] = 0;	hitBefore_Edep[7] = 0;	hitAfter_Edep[7] = 0;	hitBetween_Edep[7] = 0;

		// CHECK COINCIDENCE
		//if(NaI_etot>0 && etot>0 && nsteps>0){
		// do not check coincidence for Co57
		if(etot>0 && nsteps>0){

		LXe_Edep = etot;
		NaI_Edep = NaI_etot;

		LXe_ns = 0;
		NaI_ns = 0;
		LXe_EkinDeltaTot = 0;
		NaI_tm = 0;
		NaI_x = 0;
		NaI_y = 0;
		NaI_z = 0;

		primary_x = xp_pri;
		primary_y = yp_pri;
		primary_z = zp_pri;

		//----------------------------------------------------------	
		// LXe
		for(int j=0; j<nsteps; j++){
			//if((*parentid)[j]==0 && (*ed)[j]>0){
			if((*ed)[j]>0){
				// kinetic energy before and after each interaction
				LXe_EkinPre[LXe_ns]		= (*pre_kinE)[j];
				LXe_EkinPost[LXe_ns]	= (*post_kinE)[j];
				// positions of each interaction
				LXe_x[LXe_ns] = (*xp)[j];
				LXe_y[LXe_ns] = (*yp)[j];
				LXe_z[LXe_ns] = (*zp)[j];
				// kinetic energy before particle enters LXe
				LXe_EkinMax = (*pre_kinE)[0];
				// kinetic energy after all interactions in LXe
				if((*post_kinE)[j]<(*post_kinE)[j-1]){
				LXe_EkinMin = (*post_kinE)[j];
				}
				// time of each interaction
				LXe_tm[LXe_ns] = (*time)[j];
				// time of the 1st interaction in LXe
				LXe_tmMin = (*time)[0];
				// time of the last interaction in LXe
				if((*time)[j]>(*time)[j-1]){
				LXe_tmMax = (*time)[j];
				}
				// compute difference in kinetic energy before and after all steps, sum up
				LXe_EkinDelta[LXe_ns] = (*pre_kinE)[j]-(*post_kinE)[j];
				LXe_EkinDeltaTot += LXe_EkinDelta[LXe_ns];
				// count real number of scatters in LXe				
				LXe_ns++;
			}
		} // end loop on scatters


		//----------------------------------------------------------	
		// NaI
		NaI_tmMin = (*NaI_time)[0];
		// compute energy weighted time and positions
		for(int j_NaI=0; j_NaI<NaI_nsteps; j_NaI++){			
			if((*NaI_ed)[j_NaI]>0){
				NaI_tm	+= (*NaI_time)[j_NaI] 	* (*NaI_ed)[j_NaI];
				NaI_x	+= (*NaI_xp)[j_NaI] 	* (*NaI_ed)[j_NaI];
				NaI_y	+= (*NaI_yp)[j_NaI] 	* (*NaI_ed)[j_NaI];
				NaI_z	+= (*NaI_zp)[j_NaI] 	* (*NaI_ed)[j_NaI];
			}
		}
		NaI_tm	= NaI_tm / NaI_etot;
		NaI_x	= NaI_x	 / NaI_etot;
		NaI_y	= NaI_y  / NaI_etot;
		NaI_z	= NaI_z	 / NaI_etot;

		//----------------------------------------------------------	
		// DeadXe
		for(int j_DeadXe=0; j_DeadXe<DeadXe_nsteps; j_DeadXe++){			
			if((*DeadXe_ed)[j_DeadXe]>0 && (*DeadXe_time)[j_DeadXe]>LXe_tmMax && (*DeadXe_time)[j_DeadXe]<NaI_tmMin){
			hit_Edep[1] += (*DeadXe_ed)[j_DeadXe];
			}
			if((*DeadXe_ed)[j_DeadXe]>0 && (*DeadXe_time)[j_DeadXe]<LXe_tmMin){
			hitBefore_Edep[1] += (*DeadXe_ed)[j_DeadXe];
			}
			if((*DeadXe_ed)[j_DeadXe]>0 && (*DeadXe_time)[j_DeadXe]>NaI_tmMax){
			hitAfter_Edep[1] += (*DeadXe_ed)[j_DeadXe];
			}
			if((*DeadXe_ed)[j_DeadXe]>0 && (*DeadXe_time)[j_DeadXe]>LXe_tmMin && (*DeadXe_time)[j_DeadXe]<LXe_tmMax){
			hitBetween_Edep[1] += (*DeadXe_ed)[j_DeadXe];
			}
		}
		if(hit_Edep[1]>0) hit[1]=1;
		else hit[1]=0;
		if(hitBefore_Edep[1]>0) hitBefore[1]=1;
		else hitBefore[1]=0;
		if(hitAfter_Edep[1]>0) hitAfter[1]=1;
		else hitAfter[1]=0;
		if(hitBetween_Edep[1]>0) hitBetween[1]=1;
		else hitBetween[1]=0;

		//----------------------------------------------------------	
		// GXe
		for(int j_GXe=0; j_GXe<GXe_nsteps; j_GXe++){			
			if((*GXe_ed)[j_GXe]>0 && (*GXe_time)[j_GXe]>LXe_tmMax && (*GXe_time)[j_GXe]<NaI_tmMin){
			hit_Edep[2] += (*GXe_ed)[j_GXe];
			}
			if((*GXe_ed)[j_GXe]>0 && (*GXe_time)[j_GXe]<LXe_tmMin){
			hitBefore_Edep[2] += (*GXe_ed)[j_GXe];
			}
			if((*GXe_ed)[j_GXe]>0 && (*GXe_time)[j_GXe]>NaI_tmMax){
			hitAfter_Edep[2] += (*GXe_ed)[j_GXe];
			}
			if((*GXe_ed)[j_GXe]>0 && (*GXe_time)[j_GXe]>LXe_tmMin && (*GXe_time)[j_GXe]<LXe_tmMax){
			hitBetween_Edep[2] += (*GXe_ed)[j_GXe];
			}
		}
		if(hit_Edep[2]>0) hit[2]=1;
		else hit[2]=0;
		if(hitBefore_Edep[2]>0) hitBefore[2]=1;
		else hitBefore[2]=0;
		if(hitAfter_Edep[2]>0) hitAfter[2]=1;
		else hitAfter[2]=0;
		if(hitBetween_Edep[2]>0) hitBetween[2]=1;
		else hitBetween[2]=0;

		//----------------------------------------------------------	
		// Teflon
		for(int j_Teflon=0; j_Teflon<Teflon_nsteps; j_Teflon++){			
			if((*Teflon_ed)[j_Teflon]>0 && (*Teflon_time)[j_Teflon]>LXe_tmMax && (*Teflon_time)[j_Teflon]<NaI_tmMin){
			hit_Edep[3] += (*Teflon_ed)[j_Teflon];
			}
			if((*Teflon_ed)[j_Teflon]>0 && (*Teflon_time)[j_Teflon]<LXe_tmMin){
			hitBefore_Edep[3] += (*Teflon_ed)[j_Teflon];
			}
			if((*Teflon_ed)[j_Teflon]>0 && (*Teflon_time)[j_Teflon]>NaI_tmMax){
			hitAfter_Edep[3] += (*Teflon_ed)[j_Teflon];
			}
			if((*Teflon_ed)[j_Teflon]>0 && (*Teflon_time)[j_Teflon]>LXe_tmMin && (*Teflon_time)[j_Teflon]<LXe_tmMax){
			hitBetween_Edep[3] += (*Teflon_ed)[j_Teflon];
			}
		}
		if(hit_Edep[3]>0) hit[3]=1;
		else hit[3]=0;
		if(hitBefore_Edep[3]>0) hitBefore[3]=1;
		else hitBefore[3]=0;
		if(hitAfter_Edep[3]>0) hitAfter[3]=1;
		else hitAfter[3]=0;
		if(hitBetween_Edep[3]>0) hitBetween[3]=1;
		else hitBetween[3]=0;

		//----------------------------------------------------------	
		// Steel
		for(int j_Steel=0; j_Steel<Steel_nsteps; j_Steel++){			
			if((*Steel_ed)[j_Steel]>0 && (*Steel_time)[j_Steel]>LXe_tmMax && (*Steel_time)[j_Steel]<NaI_tmMin){
			hit_Edep[4] += (*Steel_ed)[j_Steel];
			}
			if((*Steel_ed)[j_Steel]>0 && (*Steel_time)[j_Steel]<LXe_tmMin){
			hitBefore_Edep[4] += (*Steel_ed)[j_Steel];
			}
			if((*Steel_ed)[j_Steel]>0 && (*Steel_time)[j_Steel]>NaI_tmMax){
			hitAfter_Edep[4] += (*Steel_ed)[j_Steel];
			}
			if((*Steel_ed)[j_Steel]>0 && (*Steel_time)[j_Steel]>LXe_tmMin && (*Steel_time)[j_Steel]<LXe_tmMax){
			hitBetween_Edep[4] += (*Steel_ed)[j_Steel];
			}
		}
		if(hit_Edep[4]>0) hit[4]=1;
		else hit[4]=0;
		if(hitBefore_Edep[4]>0) hitBefore[4]=1;
		else hitBefore[4]=0;
		if(hitAfter_Edep[4]>0) hitAfter[4]=1;
		else hitAfter[4]=0;
		if(hitBetween_Edep[4]>0) hitBetween[4]=1;
		else hitBetween[4]=0;

		//----------------------------------------------------------	
		// Lead
		for(int j_Lead=0; j_Lead<Lead_nsteps; j_Lead++){			
			if((*Lead_ed)[j_Lead]>0 && (*Lead_time)[j_Lead]>LXe_tmMax && (*Lead_time)[j_Lead]<NaI_tmMin){
			hit_Edep[5] += (*Lead_ed)[j_Lead];
			}
			if((*Lead_ed)[j_Lead]>0 && (*Lead_time)[j_Lead]<LXe_tmMin){
			hitBefore_Edep[5] += (*Lead_ed)[j_Lead];
			}
			if((*Lead_ed)[j_Lead]>0 && (*Lead_time)[j_Lead]>NaI_tmMax){
			hitAfter_Edep[5] += (*Lead_ed)[j_Lead];
			}
			if((*Lead_ed)[j_Lead]>0 && (*Lead_time)[j_Lead]>LXe_tmMin && (*Lead_time)[j_Lead]<LXe_tmMax){
			hitBetween_Edep[5] += (*Lead_ed)[j_Lead];
			}
		}
		if(hit_Edep[5]>0) hit[5]=1;
		else hit[5]=0;
		if(hitBefore_Edep[5]>0) hitBefore[5]=1;
		else hitBefore[5]=0;
		if(hitAfter_Edep[5]>0) hitAfter[5]=1;
		else hitAfter[5]=0;
		if(hitBetween_Edep[5]>0) hitBetween[5]=1;
		else hitBetween[5]=0;

		//----------------------------------------------------------	
		// Al
		for(int j_Al=0; j_Al<Al_nsteps; j_Al++){			
			if((*Al_ed)[j_Al]>0 && (*Al_time)[j_Al]>LXe_tmMax && (*Al_time)[j_Al]<NaI_tmMin){
			hit_Edep[6] += (*Al_ed)[j_Al];
			}
			if((*Al_ed)[j_Al]>0 && (*Al_time)[j_Al]<LXe_tmMin){
			hitBefore_Edep[6] += (*Al_ed)[j_Al];
			}
			if((*Al_ed)[j_Al]>0 && (*Al_time)[j_Al]>NaI_tmMax){
			hitAfter_Edep[6] += (*Al_ed)[j_Al];
			}
			if((*Al_ed)[j_Al]>0 && (*Al_time)[j_Al]>LXe_tmMin && (*Al_time)[j_Al]<LXe_tmMax){
			hitBetween_Edep[6] += (*Al_ed)[j_Al];
			}
		}
		if(hit_Edep[6]>0) hit[6]=1;
		else hit[6]=0;
		if(hitBefore_Edep[6]>0) hitBefore[6]=1;
		else hitBefore[6]=0;
		if(hitAfter_Edep[6]>0) hitAfter[6]=1;
		else hitAfter[6]=0;
		if(hitBetween_Edep[6]>0) hitBetween[6]=1;
		else hitBetween[6]=0;

		//----------------------------------------------------------	
		// Rest
		for(int j_Rest=0; j_Rest<Rest_nsteps; j_Rest++){			
			if((*Rest_ed)[j_Rest]>0 && (*Rest_time)[j_Rest]>LXe_tmMax && (*Rest_time)[j_Rest]<NaI_tmMin){
			hit_Edep[6] += (*Rest_ed)[j_Rest];
			}
			if((*Rest_ed)[j_Rest]>0 && (*Rest_time)[j_Rest]<LXe_tmMin){
			hitBefore_Edep[6] += (*Rest_ed)[j_Rest];
			}
			if((*Rest_ed)[j_Rest]>0 && (*Rest_time)[j_Rest]>NaI_tmMax){
			hitAfter_Edep[6] += (*Rest_ed)[j_Rest];
			}
			if((*Rest_ed)[j_Rest]>0 && (*Rest_time)[j_Rest]>LXe_tmMin && (*Rest_time)[j_Rest]<LXe_tmMax){
			hitBetween_Edep[6] += (*Rest_ed)[j_Rest];
			}
		}
		if(hit_Edep[7]>0) hit[7]=1;
		else hit[7]=0;
		if(hitBefore_Edep[7]>0) hitBefore[7]=1;
		else hitBefore[7]=0;
		if(hitAfter_Edep[7]>0) hitAfter[7]=1;
		else hitAfter[7]=0;
		if(hitBetween_Edep[7]>0) hitBetween[7]=1;
		else hitBetween[7]=0;


		//----------------------------------------------------------	
		// sum up
		hit_Edep[0] = hit_Edep[1] + hit_Edep[2] + hit_Edep[3] + hit_Edep[4] + hit_Edep[5] + hit_Edep[6] + hit_Edep[7]; 
		if(hit_Edep[0]>0) hit[0]=1;
		else hit[0]=0;
		// some up before
		hitBefore_Edep[0] = hitBefore_Edep[1] + hitBefore_Edep[2] + hitBefore_Edep[3] + hitBefore_Edep[4] + hitBefore_Edep[5] + hitBefore_Edep[6] + hitBefore_Edep[7]; 
		if(hitBefore_Edep[0]>0) hitBefore[0]=1;
		else hitBefore[0]=0;
		// some up after
		hitAfter_Edep[0] = hitAfter_Edep[1] + hitAfter_Edep[2] + hitAfter_Edep[3] + hitAfter_Edep[4] + hitAfter_Edep[5] + hitAfter_Edep[6] + hitAfter_Edep[7]; 
		if(hitAfter_Edep[0]>0) hitAfter[0]=1;
		else hitAfter[0]=0;
		// some up between
		hitBetween_Edep[0] = hitBetween_Edep[1] + hitBetween_Edep[2] + hitBetween_Edep[3] + hitBetween_Edep[4] + hitBetween_Edep[5] + hitBetween_Edep[6] + hitBetween_Edep[7]; 
		if(hitBetween_Edep[0]>0) hitBetween[0]=1;
		else hitBetween[0]=0;
		
		//if(fabs(etot-LXe_EkinDeltaTot)>1e-3){
		//if(fabs(etot-LXe_EkinDeltaTot)<1e-3){
			//cout <<"> Bad event:  "<< eventid << endl;
			//cout <<"> Bad event:  "<< eventid <<";   hitBetween = "<< hitBetween[0] <<";  hit = "<< hit[0] << endl;
			//goodEvent=0;
			//eventsRemoved ++;
			
		//}
		//else 
		t2->Fill();		
		} // END COINCIDENCE CONDITION		
	} // end loop on events
	cout <<">> Removed  "<< eventsRemoved <<"  events."<< endl;

	t2->Write();
	outFile->Close();
	

}
