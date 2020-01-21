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

#include "vector"
#include "string"
#include "iostream"
#include "fstream"
#include "stdio.h"
#include "stdlib.h"


void readCompton()
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


/*	// 34deg, no Lead Pipe
	TFile *tfile1 = new TFile("./plots/34deg_noLeadPipe.root","READ");

	TH1F *file1_E1all_1 = (TH1F*)tfile1->Get("h1_E1all_1");	// energy in NaI, all events
	TH1F *file1_E1_1 	= (TH1F*)tfile1->Get("h1_E1_1");	// energy in NaI, all coincident events
	TH1F *file1_E1_2 	= (TH1F*)tfile1->Get("h1_E1_2");	// energy in NaI, true coincidence
	TH1F *file1_E2_1 	= (TH1F*)tfile1->Get("h1_E2_1");	// energy in LXe, all coincident events
	TH1F *file1_E2_2 	= (TH1F*)tfile1->Get("h1_E2_2");	// energy in LXe, all coincident events
	TH1F *file1_EE_1 	= (TH1F*)tfile1->Get("h2_EE_1");	// energy in NaI vs LXe energy, all coincident events
	TH1F *file1_EE_2 	= (TH1F*)tfile1->Get("h2_EE_2");	// energy in NaI vs LXe energy, true coincidence
	TH1F *file1_tof_1 	= (TH1F*)tfile1->Get("h1_tof_1");	// time of flight, all coincident events
	TH1F *file1_tof_2 	= (TH1F*)tfile1->Get("h1_tof_2");	// time of flight, true coincidence
	TH1F *file1_tofE_1 	= (TH1F*)tfile1->Get("h1_tofE_1");	// NaI energy VS time of flight
	TH1F *file1_tofE_2	= (TH1F*)tfile1->Get("h1_tofE_2");	// NaI energy VS time of flight
	// 34deg, with Lead Pipe 50cm
	TFile *tfile2 = new TFile("./plots/34deg_withLeadPipe.root","READ");

	TH1F *file2_E1all_1 = (TH1F*)tfile2->Get("h1_E1all_1");	// energy in NaI, all events
	TH1F *file2_E1_1 	= (TH1F*)tfile2->Get("h1_E1_1");	// energy in NaI, all coincident events
	TH1F *file2_E1_2 	= (TH1F*)tfile2->Get("h1_E1_2");	// energy in NaI, true coincidence
	TH1F *file2_E2_1 	= (TH1F*)tfile2->Get("h1_E2_1");	// energy in LXe, all coincident events
	TH1F *file2_E2_2 	= (TH1F*)tfile2->Get("h1_E2_2");	// energy in LXe, all coincident events
	TH1F *file2_EE_1 	= (TH1F*)tfile2->Get("h2_EE_1");	// energy in NaI vs LXe energy, all coincident events
	TH1F *file2_EE_2 	= (TH1F*)tfile2->Get("h2_EE_2");	// energy in NaI vs LXe energy, true coincidence
	TH1F *file2_tof_1 	= (TH1F*)tfile2->Get("h1_tof_1");	// time of flight, all coincident events
	TH1F *file2_tof_2 	= (TH1F*)tfile2->Get("h1_tof_2");	// time of flight, true coincidence
	TH1F *file2_tofE_1 	= (TH1F*)tfile2->Get("h1_tofE_1");	// NaI energy VS time of flight
	TH1F *file2_tofE_2	= (TH1F*)tfile2->Get("h1_tofE_2");	// NaI energy VS time of flight

	file2_E1all_1	->SetLineColor(2);
	file2_E1_1 		->SetLineStyle(7);
	file2_E1_2 		->SetLineStyle(7);
	file2_E2_1 		->SetLineStyle(7);
	file2_E2_2 		->SetLineStyle(7);
	file2_tof_1 	->SetLineStyle(7);
	file2_tof_2 	->SetLineStyle(7);
*/
	// 10deg, no Lead Pipe
	TFile *tfile1 = new TFile("./plots/10deg_noLeadPipe.root","READ");

	TH1F *file1_E1all_1 = (TH1F*)tfile1->Get("h1_E1all_1");	// energy in NaI, all events
	TH1F *file1_E1_1 	= (TH1F*)tfile1->Get("h1_E1_1");	// energy in NaI, all coincident events
	TH1F *file1_E1_2 	= (TH1F*)tfile1->Get("h1_E1_2");	// energy in NaI, true coincidence
	TH1F *file1_E2_1 	= (TH1F*)tfile1->Get("h1_E2_1");	// energy in LXe, all coincident events
	TH1F *file1_E2_2 	= (TH1F*)tfile1->Get("h1_E2_2");	// energy in LXe, all coincident events
	TH1F *file1_EE_1 	= (TH1F*)tfile1->Get("h2_EE_1");	// energy in NaI vs LXe energy, all coincident events
	TH1F *file1_EE_2 	= (TH1F*)tfile1->Get("h2_EE_2");	// energy in NaI vs LXe energy, true coincidence
	TH1F *file1_tof_1 	= (TH1F*)tfile1->Get("h1_tof_1");	// time of flight, all coincident events
	TH1F *file1_tof_2 	= (TH1F*)tfile1->Get("h1_tof_2");	// time of flight, true coincidence
	TH1F *file1_tofE_1 	= (TH1F*)tfile1->Get("h1_tofE_1");	// NaI energy VS time of flight
	TH1F *file1_tofE_2	= (TH1F*)tfile1->Get("h1_tofE_2");	// NaI energy VS time of flight
	// 34deg, with Lead Pipe 50cm
	TFile *tfile2 = new TFile("./plots/10deg_withLeadPipe50cm.root","READ");

	TH1F *file2_E1all_1 = (TH1F*)tfile2->Get("h1_E1all_1");	// energy in NaI, all events
	TH1F *file2_E1_1 	= (TH1F*)tfile2->Get("h1_E1_1");	// energy in NaI, all coincident events
	TH1F *file2_E1_2 	= (TH1F*)tfile2->Get("h1_E1_2");	// energy in NaI, true coincidence
	TH1F *file2_E2_1 	= (TH1F*)tfile2->Get("h1_E2_1");	// energy in LXe, all coincident events
	TH1F *file2_E2_2 	= (TH1F*)tfile2->Get("h1_E2_2");	// energy in LXe, all coincident events
	TH1F *file2_EE_1 	= (TH1F*)tfile2->Get("h2_EE_1");	// energy in NaI vs LXe energy, all coincident events
	TH1F *file2_EE_2 	= (TH1F*)tfile2->Get("h2_EE_2");	// energy in NaI vs LXe energy, true coincidence
	TH1F *file2_tof_1 	= (TH1F*)tfile2->Get("h1_tof_1");	// time of flight, all coincident events
	TH1F *file2_tof_2 	= (TH1F*)tfile2->Get("h1_tof_2");	// time of flight, true coincidence
	TH1F *file2_tofE_1 	= (TH1F*)tfile2->Get("h1_tofE_1");	// NaI energy VS time of flight
	TH1F *file2_tofE_2	= (TH1F*)tfile2->Get("h1_tofE_2");	// NaI energy VS time of flight
	// 34deg, with Lead Pipe 75cm
	TFile *tfile3 = new TFile("./plots/10deg_withLeadPipe75cm.root","READ");

	TH1F *file3_E1all_1 = (TH1F*)tfile3->Get("h1_E1all_1");	// energy in NaI, all events
	TH1F *file3_E1_1 	= (TH1F*)tfile3->Get("h1_E1_1");	// energy in NaI, all coincident events
	TH1F *file3_E1_2 	= (TH1F*)tfile3->Get("h1_E1_2");	// energy in NaI, true coincidence
	TH1F *file3_E2_1 	= (TH1F*)tfile3->Get("h1_E2_1");	// energy in LXe, all coincident events
	TH1F *file3_E2_2 	= (TH1F*)tfile3->Get("h1_E2_2");	// energy in LXe, all coincident events
	TH1F *file3_EE_1 	= (TH1F*)tfile3->Get("h2_EE_1");	// energy in NaI vs LXe energy, all coincident events
	TH1F *file3_EE_2 	= (TH1F*)tfile3->Get("h2_EE_2");	// energy in NaI vs LXe energy, true coincidence
	TH1F *file3_tof_1 	= (TH1F*)tfile3->Get("h1_tof_1");	// time of flight, all coincident events
	TH1F *file3_tof_2 	= (TH1F*)tfile3->Get("h1_tof_2");	// time of flight, true coincidence
	TH1F *file3_tofE_1 	= (TH1F*)tfile3->Get("h1_tofE_1");	// NaI energy VS time of flight
	TH1F *file3_tofE_2	= (TH1F*)tfile3->Get("h1_tofE_2");	// NaI energy VS time of flight

	file1_E1all_1	->SetLineColor(1);
	file2_E1all_1	->SetLineColor(4);
	file3_E1all_1	->SetLineColor(2);

	file1_E1_1 		->SetLineColor(1);
	file2_E1_1 		->SetLineColor(4);
	file3_E1_1 		->SetLineColor(2);

	file1_E1_2 		->SetLineColor(1);
	file2_E1_2 		->SetLineColor(4);
	file3_E1_2 		->SetLineColor(2);

	file1_E2_1 		->SetLineColor(1);
	file2_E2_1 		->SetLineColor(4);
	file3_E2_1 		->SetLineColor(2);

	file1_E2_2 		->SetLineColor(1);
	file2_E2_2 		->SetLineColor(4);
	file3_E2_2 		->SetLineColor(2);

	file2_tof_1 	->SetLineStyle(7);
	file2_tof_2 	->SetLineStyle(7);

	file3_tof_1 	->SetLineStyle(10);
	file3_tof_2 	->SetLineStyle(10);
	
	// energy in NaI, all events
	TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 700, 500);
			 c1->SetFillColor(10);
			 c1->SetLogy();
		file1_E1all_1	->Draw();
		file2_E1all_1	->Draw("same");
		file3_E1all_1	->Draw("same");
		file1_E1all_1	->SetMinimum(100);
		file1_E1all_1	->GetXaxis()->SetTitle("Energy deposited in NaI [keV]");
		file1_E1all_1	->GetYaxis()->SetTitle("Entries");
		file1_E1all_1	->GetXaxis()->CenterTitle(true);
		file1_E1all_1	->GetYaxis()->CenterTitle(true);
		file1_E1all_1	->GetXaxis()->SetLabelSize(0.03);
		file1_E1all_1	->GetYaxis()->SetLabelSize(0.03);
		TLegend *leg1 = new TLegend(0.1408046,0.7330508,0.7140805,0.8834746,NULL,"brNDC");
				 leg1 ->SetBorderSize(0);
				 leg1 ->SetFillColor(10);
				 leg1->AddEntry(file1_E1all_1,	"standard NaI shield", 	"L");
				 leg1->AddEntry(file2_E1all_1,	"extended NaI shield (50cm channel)", 	"L");
				 leg1->AddEntry(file3_E1all_1,	"extended NaI shield (75cm channel)", 	"L");
				 leg1 ->Draw();
	c1->Update();

/*	// energy in NaI, all coincident events
	TCanvas *c2a = new TCanvas("c2a", "NaI, all coincident events", 0, 0, 700, 500);
			 c2a->SetFillColor(10);
			 c2a->SetLogy();
		file1_E1_1	->Draw();
		file2_E1_1	->Draw("same");
		file3_E1_1	->Draw("same");
		
		file1_E1_1	->Rebin(2);
		file2_E1_1	->Rebin(2);
		file3_E1_1	->Rebin(2);

		file1_E1_1	->Scale(0.5);
		file2_E1_1	->Scale(0.5);
		file3_E1_1	->Scale(0.5);
		
		file1_E1_1	->GetXaxis()->SetTitle("Energy deposited in NaI [keV]");
		file1_E1_1	->GetYaxis()->SetTitle("Entries");
		file1_E1_1	->GetXaxis()->CenterTitle(true);
		file1_E1_1	->GetYaxis()->CenterTitle(true);
		file1_E1_1	->GetXaxis()->SetLabelSize(0.03);
		file1_E1_1	->GetYaxis()->SetLabelSize(0.03);
		TLegend *leg2a = new TLegend(0.1408046,0.7330508,0.7140805,0.8834746,NULL,"brNDC");
				 leg2a ->SetBorderSize(0);
				 leg2a ->SetFillColor(10);
				 leg2a->AddEntry(file1_E1_1,	"standard NaI shield", 	"L");
				 leg2a->AddEntry(file2_E1_1,	"extended NaI shield (50cm channel)", 	"L");
				 leg2a->AddEntry(file3_E1_1,	"extended NaI shield (75cm channel)", 	"L");
				 leg2a ->Draw();
	c2a->Update();

	// energy in NaI, true coincident events
	TCanvas *c2b = new TCanvas("c2b", "NaI, true coincidence", 0, 0, 700, 500);
			 c2b->SetFillColor(10);
			 c2b->SetLogy();
		file1_E1_2	->Draw();
		file2_E1_2	->Draw("same");
		file3_E1_2	->Draw("same");
		
/*		file1_E1_2	->Rebin(2);
		file2_E1_2	->Rebin(2);
		file3_E1_2	->Rebin(2);

		file1_E1_2	->Scale(0.5);
		file2_E1_2	->Scale(0.5);
		file3_E1_2	->Scale(0.5);
*/		
		file1_E1_2	->GetXaxis()->SetTitle("Energy deposited in NaI [keV]");
		file1_E1_2	->GetYaxis()->SetTitle("Entries");
		file1_E1_2	->GetXaxis()->CenterTitle(true);
		file1_E1_2	->GetYaxis()->CenterTitle(true);
		file1_E1_2	->GetXaxis()->SetLabelSize(0.03);
		file1_E1_2	->GetYaxis()->SetLabelSize(0.03);
		TLegend *leg2b = new TLegend(0.1408046,0.7330508,0.7140805,0.8834746,NULL,"brNDC");
				 leg2b ->SetBorderSize(0);
				 leg2b ->SetFillColor(10);
				 leg2b->AddEntry(file1_E1_2,	"standard NaI shield", 	"L");
				 leg2b->AddEntry(file2_E1_2,	"extended NaI shield (50cm channel)", 	"L");
				 leg2b->AddEntry(file3_E1_2,	"extended NaI shield (75cm channel)", 	"L");
				 leg2b ->Draw();
	c2b->Update();

	// energy in LXe, all coincident events
	TCanvas *c3a = new TCanvas("c3a", "LXe, all coincident events", 0, 0, 700, 500);
			 c3a->SetFillColor(10);
			 c3a->SetLogy();
		file1_E2_1	->Draw();
		file2_E2_1	->Draw("same");
		file3_E2_1	->Draw("same");
		
		file1_E2_1	->Rebin(2);
		file2_E2_1	->Rebin(2);
		file3_E2_1	->Rebin(2);

		file1_E2_1	->Scale(0.5);
		file2_E2_1	->Scale(0.5);
		file3_E2_1	->Scale(0.5);
		
		file1_E2_1	->GetXaxis()->SetTitle("Energy deposited in LXe [keV]");
		file1_E2_1	->GetYaxis()->SetTitle("Entries");
		file1_E2_1	->GetXaxis()->CenterTitle(true);
		file1_E2_1	->GetYaxis()->CenterTitle(true);
		file1_E2_1	->GetXaxis()->SetLabelSize(0.03);
		file1_E2_1	->GetYaxis()->SetLabelSize(0.03);
		TLegend *leg3a = new TLegend(0.1408046,0.7330508,0.7140805,0.8834746,NULL,"brNDC");
				 leg3a ->SetBorderSize(0);
				 leg3a ->SetFillColor(10);
				 leg3a->AddEntry(file1_E2_1,	"standard NaI shield", 	"L");
				 leg3a->AddEntry(file2_E2_1,	"extended NaI shield (50cm channel)", 	"L");
				 leg3a->AddEntry(file3_E2_1,	"extended NaI shield (75cm channel)", 	"L");
				 leg3a ->Draw();
	c3a->Update();


	// energy in LXe, coincident events
	TCanvas *c3b = new TCanvas("c3b", "LXe, true coincidence", 0, 0, 700, 500);
			 c3b->SetFillColor(10);
			 c3b->SetLogy();
		file1_E2_2	->Draw();
		file2_E2_2	->Draw("same");
		file3_E2_2	->Draw("same");
		
		file1_E2_2	->Rebin(2);
		file2_E2_2	->Rebin(2);
		file3_E2_2	->Rebin(2);

		file1_E2_2	->Scale(0.5);
		file2_E2_2	->Scale(0.5);
		file3_E2_2	->Scale(0.5);
		
		file1_E2_2	->GetXaxis()->SetTitle("Energy deposited in LXe [keV]");
		file1_E2_2	->GetYaxis()->SetTitle("Entries");
		file1_E2_2	->GetXaxis()->CenterTitle(true);
		file1_E2_2	->GetYaxis()->CenterTitle(true);
		file1_E2_2	->GetXaxis()->SetLabelSize(0.03);
		file1_E2_2	->GetYaxis()->SetLabelSize(0.03);
		TLegend *leg3b = new TLegend(0.1408046,0.7330508,0.7140805,0.8834746,NULL,"brNDC");
				 leg3b ->SetBorderSize(0);
				 leg3b ->SetFillColor(10);
				 leg3b->AddEntry(file1_E2_2,	"standard NaI shield", 	"L");
				 leg3b->AddEntry(file2_E2_2,	"extended NaI shield (50cm channel)", 	"L");
				 leg3b->AddEntry(file3_E2_2,	"extended NaI shield (75cm channel)", 	"L");
				 leg3b ->Draw();
	c3b->Update();

	const double livetime = 226350; // sec
	double n_NaI_allEvents_1 = file1_E1_1->GetEntries();
	double n_NaI_allEvents_2 = file2_E1_1->GetEntries();
	double n_NaI_allEvents_3 = file3_E1_1->GetEntries();
	cout <<"NaI, all events:           "<< n_NaI_allEvents_1/livetime <<" Hz"<< endl;
	cout <<"NaI, all events, pipe50cm: "<< n_NaI_allEvents_2/livetime <<" Hz"<< endl;
	cout <<"NaI, all events, pipe75cm: "<< n_NaI_allEvents_3/livetime <<" Hz"<< endl;
*/
}
