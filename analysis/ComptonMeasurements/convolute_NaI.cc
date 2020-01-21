#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TObject.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TF1.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TAxis.h"


#include "vector"
#include "string"
#include "iostream"
#include "fstream"
#include "stdio.h"
#include "stdlib.h"

void convolute_NaI()
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


	//NaI light yield:
	//f = m*x, m=(40+-1) a.u/keV
	
	//NaI resolution:
	//f = a/sqrt(x), a=(83+-2)sqrt(keV), b=(0.2+-0.1)
	TF1 *fres_NaI = new TF1("fres_NaI", "83/sqrt(x)+0.2", 0, 1000);
		 fres_NaI->SetLineColor(2);
		 fres_NaI->SetLineWidth(2);

	TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 1000, 500);
			 c1->SetFillColor(10);
		fres_NaI	->Draw();
		fres_NaI	->GetXaxis()->SetTitle("Resolution (#sigma/#mu) [%]");
		fres_NaI	->GetYaxis()->SetTitle("Energy [keV]");
		fres_NaI	->GetXaxis()->CenterTitle(true);
		fres_NaI	->GetYaxis()->CenterTitle(true);
		fres_NaI	->GetXaxis()->SetLabelSize(0.03);
		fres_NaI	->GetYaxis()->SetLabelSize(0.03);
	c1->Update();
	

	////////////////////
	// LOAD DATA
	////////////////////
	const char *DataFolder 	= "/Users/alexkish/data/xuerich/mc/compton";
	const char *OutFolder 	= "./plots";
	
	Char_t  DataFile1[500];
	
	sprintf(DataFile1, "%s/Xu_Compton_4.25deg_5e9.root", DataFolder);
	//sprintf(DataFile2, "%s/Xu_Compton_4.25deg_5e9c.root", DataFolder);
	TChain *t1	= new TChain("t1");
			t1	->Add(DataFile1);


	// INPUT TREE
	float NaI_etot;

	// OUTPUT TREEE
	// raw energy in keV
	float NaI_etotRes;





}
