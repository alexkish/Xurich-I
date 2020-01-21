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


void ConvertBeam()
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
	const char *DataFolder 	= "/Users/alexkish/data/xuerich/mc/beam";	// idark05x
		
	const char *inFileName = "XuerichHP-Xe_nBeam2.5MeV_1e5";
	//const char *inFileName = "XuerichHP-Ar_nBeam2.5MeV_1e5";
	
	
	Char_t  inFile[500];
	sprintf(inFile, "%s/%s.root", DataFolder, inFileName);

	TChain *t1	= new TChain("t1");
			t1	->Add(inFile);

	int eventid;
	int nsteps;
	float etotNR, etotER;
	float xp_pri;
	float yp_pri;
	float zp_pri;
	std::vector<int> *parentid 		= new vector<int>;
	//std::vector<double> *pre_kinE 	= new vector<double>;
	//std::vector<double> *post_kinE 	= new vector<double>;
	std::vector<float> *time 		= new vector<float>;
	std::vector<float> *ed 			= new vector<float>;
	std::vector<float> *xp 			= new vector<float>;
	std::vector<float> *yp 			= new vector<float>;
	std::vector<float> *zp 			= new vector<float>;
	std::vector<string> *type 		= new vector<string>;
	
	float Teflon_etot, DeadXe_etot, Rest_etot;
	

	t1->SetBranchStatus("*",0); // 1st - disable all branches to speed up processing
	t1->SetBranchAddress("eventid",			&eventid);	
	t1->SetBranchAddress("nsteps",			&nsteps);
	t1->SetBranchAddress("parentid",		&parentid);
	t1->SetBranchAddress("etotNR",			&etotNR);
	t1->SetBranchAddress("etotER",			&etotER);
	//t1->SetBranchAddress("pre_kinE",		&pre_kinE);
	//t1->SetBranchAddress("post_kinE",		&post_kinE);
	t1->SetBranchAddress("time",			&time);
	t1->SetBranchAddress("ed",				&ed);
	t1->SetBranchAddress("xp",				&xp);
	t1->SetBranchAddress("yp",				&yp);
	t1->SetBranchAddress("zp",				&zp);
	t1->SetBranchAddress("xp_pri",			&xp_pri);
	t1->SetBranchAddress("yp_pri",			&yp_pri);
	t1->SetBranchAddress("zp_pri",			&zp_pri);
	t1->SetBranchAddress("type",			&type);

	t1->SetBranchAddress("Teflon_etot",		&Teflon_etot);
	t1->SetBranchAddress("DeadXe_etot",		&DeadXe_etot);
	t1->SetBranchAddress("Rest_etot",		&Rest_etot);
	 	
	//--------------------------------------------------------------------------------------
	int k = t1->GetEntries();
	cout <<"File "<< inFileName <<" has "<< k <<" entries."<< endl;
	//--------------------------------------------------------------------------------------
	
	// OUTPUT FILE
	int LXe_nsNR, LXe_nsER;
	float LXe_EdepNR, LXe_EdepER;
	float LXe_tmNR[50], LXe_tmER[50];
	float LXe_xNR[50], LXe_xER[50];
	float LXe_yNR[50], LXe_yER[50];
	float LXe_zNR[50], LXe_zER[50];

	float LXe_edNR[50];
	float LXe_edER[50];

	float LXe_tmMinNR = 0;
	float LXe_tmMaxNR = 0;

	int hitBefore;
	float hitBefore_Edep;
	float etotOthers;
	

	string xeStr="Xe", nucName;

	Char_t  outFileName[500];
	sprintf(outFileName, "%s/%s-t2.root", DataFolder, inFileName);

	TFile *outFile = new TFile(outFileName,"RECREATE");
	TTree *t2 = new TTree("t2","t2");

	t2->Branch("LXe_nsNR", 			&LXe_nsNR,			"LXe_nsNR/I");
	t2->Branch("LXe_nsER", 			&LXe_nsER,			"LXe_nsER/I");
	t2->Branch("LXe_EdepNR", 		&LXe_EdepNR,		"LXe_EdepNR/F");
	t2->Branch("LXe_EdepER", 		&LXe_EdepER,		"LXe_EdepER/F");
	t2->Branch("LXe_edNR", 			&LXe_edNR,			"LXe_edNR[LXe_nsNR]/F");
	t2->Branch("LXe_tmNR", 			&LXe_tmNR,			"LXe_tmNR[LXe_nsNR]/F");
	t2->Branch("LXe_xNR", 			&LXe_xNR,			"LXe_xNR[LXe_nsNR]/F");
	t2->Branch("LXe_yNR", 			&LXe_yNR,			"LXe_yNR[LXe_nsNR]/F");
	t2->Branch("LXe_zNR", 			&LXe_zNR,			"LXe_zNR[LXe_nsNR]/F");
	t2->Branch("LXe_edER", 			&LXe_edER,			"LXe_edER[LXe_nsER]/F");
	t2->Branch("LXe_tmER", 			&LXe_tmER,			"LXe_tmER[LXe_nsER]/F");
	t2->Branch("LXe_xER", 			&LXe_xER,			"LXe_xER[LXe_nsER]/F");
	t2->Branch("LXe_yER", 			&LXe_yER,			"LXe_yER[LXe_nsER]/F");
	t2->Branch("LXe_zER", 			&LXe_zER,			"LXe_zER[LXe_nsER]/F");
	t2->Branch("hitBefore", 		&hitBefore,			"hitBefore/I");
	t2->Branch("etotOthers", 		&etotOthers,		"etotOthers/F");
	t2->Branch("LXe_tmMinNR", 		&LXe_tmMinNR,		"LXe_tmMinNR/F");
	t2->Branch("LXe_tmMaxNR", 		&LXe_tmMaxNR,		"LXe_tmMaxNR/F");

	
	
	// START
	for(int i=0; i<k; i++){
	//for(int i=0; i<10; i++){
	/*	if(i%100==0){
			cout << "\rProcessing event " << i << " of " << k;
			cout.flush();
		}	
	*/	t1->GetEntry(i);
		
		LXe_EdepNR = etotNR;
		LXe_EdepER = etotER;

		LXe_nsNR = 0;
		LXe_nsER = 0;

		hitBefore_Edep = 0;

		etotOthers = Teflon_etot + DeadXe_etot + Rest_etot;


		//----------------------------------------------------------	
		// LXe
		for(int j=0; j<nsteps; j++){
			//if((*parentid)[j]==0 && (*ed)[j]>0){
	
			nucName = ((*type)[j]).substr(0,2);
			//cout <<"nucName = "<< nucName << endl;
			

			// NRs
			//if(nucName=="Ar" && (*ed)[j]>0){ // for Ar
			if(nucName=="Xe" && (*ed)[j]>0){ // for Xe
				LXe_edNR[LXe_nsNR] = (*ed)[j];
				// positions of each interaction
				LXe_xNR[LXe_nsNR] = (*xp)[j];
				LXe_yNR[LXe_nsNR] = (*yp)[j];
				LXe_zNR[LXe_nsNR] = (*zp)[j];
				// time of each interaction
				LXe_tmNR[LXe_nsNR] = (*time)[j];
				// time of the 1st interaction in LXe
				LXe_tmMinNR = (*time)[0];
				// time of the last interaction in LXe
				if((*time)[j]>(*time)[j-1]){
				LXe_tmMaxNR = (*time)[j];
				}
				// number of NRs
				LXe_nsNR++;
				//cout <<"i = "<< i<<";   nsteps= "<< nsteps <<";   ed = "<< (*ed)[j] <<";   LXe_ns = "<< LXe_nsNR <<";   type = "<< (*type)[j] << endl;
			} 
			// ERs
			//else if(nucName!="Ar" && (*ed)[j]>0){
			if(nucName!="Xe" && (*ed)[j]>0){ // for Xe
				LXe_edER[LXe_nsER] = (*ed)[j];
				// positions of each interaction
				LXe_xER[LXe_nsER] = (*xp)[j];
				LXe_yER[LXe_nsER] = (*yp)[j];
				LXe_zER[LXe_nsER] = (*zp)[j];
				// time of each interaction
				LXe_tmER[LXe_nsER] = (*time)[j];
				// number of ERs
				LXe_nsER++;
				//cout <<"i = "<< i<<";   nsteps= "<< nsteps <<";   ed = "<< (*ed)[j] <<";   LXe_ns = "<< LXe_nsER <<";   type = "<< (*type)[j] << endl;
			} 
		
		} // end loop on scatters


		t2->Fill();		
		//}// // END COINCIDENCE CONDITION		
	} // end loop on events

	t2->Write();
	outFile->Close();
	

}
