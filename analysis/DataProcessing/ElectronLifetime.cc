#include <iostream>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <ctime>
#include <math.h>

#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TDatime.h"

using namespace std;

void ElectronLifetime()
{
	// settings
	gStyle->SetStatBorderSize(0);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFillColor(10);
	gStyle->SetStatColor(10);
	gStyle->SetStatFont(42);
	gStyle->SetMarkerStyle(8);
	gStyle->SetMarkerColor(4);
	gStyle->SetOptStat(0);

	char dataset[200];
	float livetime;
	float linecut_a, linecut_b;
	float S1fit_a, S1fit_b;
	float S2fit_c, S2fit_d;
	float livetime_tot = 0;
	char date[200], time[200], yr[200], mon[200], day[200], hour[200], min[200];
	Int_t yr_i[100], day_i[100], mon_i[100], hour_i[100], min_i[100];

	float x[100], xer[100], y[100], yer[100];
	
	char *outFolder1 = "/Users/jaehzorn/analysis/xeData/CalibrationSources/Cs137/plots";

	ifstream logfile;

	char logfile_char[200];
	sprintf(logfile_char,"%s/logFile.txt",outFolder1);
	logfile.open(logfile_char);

	////////////////////////
	// Read data into buffer
	////////////////////////

	int nf = 0;
	
	while((logfile >> dataset >> livetime >> linecut_a >> linecut_b >> S1fit_a >> S1fit_b >> S2fit_c >> S2fit_d )){
	yr[0] 	= dataset[6];
	yr[1] 	= dataset[7];
	//yr[2] 	= dataset[17];
	
	mon[0] 	= dataset[8];
	mon[1] 	= dataset[9];
	mon[2] 	= dataset[17];
	
	day[0] 	= dataset[10];
	day[1] 	= dataset[11];
	day[2] 	= dataset[17];

	hour[0] = dataset[13];
	hour[1]	= dataset[14];
	hour[2]	= dataset[17];
	
	min[0]	= dataset[15];
	min[1]	= dataset[16];
	min[2]	= dataset[17];
	
	char temp_yr[100];
	sprintf(temp_yr,"20%s",yr);
	char temp_mon[100];
	sprintf(temp_mon,"%s",mon);
	char temp_day[100];
	sprintf(temp_day,"%s",day);
	char temp_hour[100];
	sprintf(temp_hour,"%s",hour);
	char temp_min[100];
	sprintf(temp_min,"%s",min);

	yr_i[nf] 	= atoi(temp_yr);
	mon_i[nf] 	= atoi(temp_mon);
	day_i[nf] 	= atoi(temp_day);
	hour_i[nf] 	= atoi(temp_hour);
	min_i[nf] 	= atoi(temp_min);
		
	TDatime da(yr_i[nf],mon_i[nf],day_i[nf],hour_i[nf],min_i[nf],00); 
	x[nf] = da.Convert();
	y[nf] = TMath::Abs(1/S2fit_d);
	
	//cout << yr <<" "<< mon <<" "<< day <<" "<< hour <<" "<< min << endl;
	
	nf++;

	}
	logfile.close();


	TCanvas *c_elifetime = new TCanvas("c_elifetime", "Electron Lifetime", 0, 0, 700, 500);
			 c_elifetime->SetFillColor(10);
			 c_elifetime->SetGrid();
			 c_elifetime->SetFillStyle(0);

    TGraph *gr = new TGraph(nf, x, y);
    		gr->SetMarkerStyle(22);
    		gr->SetMarkerSize(1.3);
    		gr->SetMarkerColor(4);
    		gr->SetLineColor(2);
    		gr->SetTitle("");
			gr->GetXaxis()->SetTimeFormat("%m/%d");
			gr->GetXaxis()->SetTimeDisplay(1);
			gr->GetXaxis()->SetTitle("Time");
    		gr->GetYaxis()->SetTitle("Electron lifetime [#mus]");
			gr->GetXaxis()->CenterTitle(true);
			gr->GetYaxis()->CenterTitle(true);
    		//gain_er->GetYaxis()->SetRangeUser(0.0,3.0);
			gr->Draw("AP");

	c_elifetime->Update();

	Char_t  name_elifetime[500];
	sprintf(name_elifetime, "%s/e_lifetime.png", outFolder1);
	c_elifetime->Print(name_elifetime);

}
