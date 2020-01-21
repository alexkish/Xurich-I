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


void readProCompton()
{	
	gStyle	->SetOptStat(1);
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
	const char *DataFolder 	= "/Users/alexkish/data/xuerich/mc/compton";	// idark05x
	const char *OutFolder 	= "./plots";
	
	Char_t  DataFile1[500];
	Char_t  DataFile2[500];
	
	//sprintf(DataFile1, "%s/Xu_Compton_4.25deg-70cm_4e10-t2.root", DataFolder);
	//sprintf(DataFile1, "%s/Xu_Compton_4.25deg-70cm-Air_4e10-t2.root", DataFolder);
	//sprintf(DataFile1, "%s/Xu_Compton_4.25deg-70cm-Vacuum_4e10-t2.root", DataFolder);

	sprintf(DataFile1, "%s/Xu_Compton_6.25deg-70cm-Vacuum_4e11-t2.root", DataFolder);

	TChain *t2	= new TChain("t2");
			t2	->Add(DataFile1);


	////////////////////////////////////////////////////////////////
	// ALIASES, CUTS
	t2->SetAlias("TOF", "(NaI_tm-LXe_tm)*1e9");	

	TCut cutTrue1 = "LXe_Edep<40 && LXe_ns==1 && hit[0]==0 && hitBefore[0]==0 && hitAfter[0]==0";
	TCut cutTrue2 = "LXe_Edep>40 && LXe_ns==1 && hit[0]==0 && hitBefore[0]==0 && hitAfter[0]==0";
		
	////////////////////////////////////////////////////////////////
	// HISTOGRAMS
	////////////////////////////////////////////////////////////////

	int nbins_tof 	= 400;
	double min_tof 	= 0.0;
	double max_tof 	= 6.0;

	TH1F *hTOF_1 = new TH1F("hTOF_1", "", nbins_tof, min_tof, max_tof);
		  hTOF_1 ->SetLineColor(1);
		  hTOF_1 ->SetLineWidth(2);
	TH1F *hTOF_2 = new TH1F("hTOF_2", "", nbins_tof, min_tof, max_tof);
		  hTOF_2 ->SetLineColor(2);
		  hTOF_2 ->SetLineWidth(2);
	TH1F *hTOF_3 = new TH1F("hTOF_3", "", nbins_tof, min_tof, max_tof);
		  hTOF_3 ->SetLineColor(3);
		  hTOF_3 ->SetLineWidth(2);
	TH1F *hTOF_4 = new TH1F("hTOF_4", "", nbins_tof, min_tof, max_tof);
		  hTOF_4 ->SetLineColor(4);
		  hTOF_4 ->SetLineWidth(2);

	int nbinsE_NaI 	= 1400;
	double minE_NaI = 0;
	double maxE_NaI = 700;

	TH1F *hE_NaI_1 	= new TH1F("hE_NaI_1", "", nbinsE_NaI, minE_NaI, maxE_NaI);
		  hE_NaI_1 	->SetLineColor(1);
		  hE_NaI_1 	->SetLineWidth(2);
	TH1F *hE_NaI_2 	= new TH1F("hE_NaI_2", "", nbinsE_NaI, minE_NaI, maxE_NaI);
		  hE_NaI_2 	->SetLineColor(2);
		  hE_NaI_2 	->SetLineWidth(2);
	TH1F *hE_NaI_3 	= new TH1F("hE_NaI_3", "", nbinsE_NaI, minE_NaI, maxE_NaI);
		  hE_NaI_3 	->SetLineColor(3);
		  hE_NaI_3 	->SetLineWidth(2);
	TH1F *hE_NaI_4 	= new TH1F("hE_NaI_4", "", nbinsE_NaI, minE_NaI, maxE_NaI);
		  hE_NaI_4 	->SetLineColor(4);
		  hE_NaI_4 	->SetLineWidth(2);

	TH2F *hETOF_NaI_1 	= new TH2F("hETOF_NaI_1", "", 1000, min_tof, max_tof, 1000, minE_NaI, maxE_NaI);
		  hETOF_NaI_1 	->SetMarkerColor(1);
	TH2F *hETOF_NaI_2 	= new TH2F("hETOF_NaI_2", "", 1000, min_tof, max_tof, 1000, minE_NaI, maxE_NaI);
		  hETOF_NaI_2 	->SetMarkerColor(2);
	TH2F *hETOF_NaI_3 	= new TH2F("hETOF_NaI_3", "", 1000, min_tof, max_tof, 1000, minE_NaI, maxE_NaI);
		  hETOF_NaI_3 	->SetMarkerColor(3);
	TH2F *hETOF_NaI_4 	= new TH2F("hETOF_NaI_4", "", 1000, min_tof, max_tof, 1000, minE_NaI, maxE_NaI);
		  hETOF_NaI_4 	->SetMarkerColor(4);

	int nbinsE_LXe 	= 1200;
	double minE_LXe = 0;
	double maxE_LXe = 700;

	TH1F *hE_LXe_1 	= new TH1F("hE_LXe_1", "", nbinsE_LXe, minE_LXe, maxE_LXe);
		  hE_LXe_1 	->SetLineColor(1);
		  hE_LXe_1 	->SetLineWidth(2);
	TH1F *hE_LXe_2 	= new TH1F("hE_LXe_2", "", nbinsE_LXe, minE_LXe, maxE_LXe);
		  hE_LXe_2 	->SetLineColor(2);
		  hE_LXe_2 	->SetLineWidth(2);
	TH1F *hE_LXe_3	= new TH1F("hE_LXe_3", "", nbinsE_LXe, minE_LXe, maxE_LXe);
		  hE_LXe_3	->SetLineColor(3);
		  hE_LXe_3	->SetLineWidth(2);
	TH1F *hE_LXe_4	= new TH1F("hE_LXe_4", "", nbinsE_LXe, minE_LXe, maxE_LXe);
		  hE_LXe_4	->SetLineColor(4);
		  hE_LXe_4	->SetLineWidth(2);

	TH2F *hE_NaI_LXe_1 	= new TH2F("hE_NaI_LXe_1", "", 1000, minE_LXe, maxE_LXe, 1000, minE_NaI, maxE_NaI);
		  hE_NaI_LXe_1 	->SetMarkerColor(1);
	TH2F *hE_NaI_LXe_2 	= new TH2F("hE_NaI_LXe_2", "", 1000, minE_LXe, maxE_LXe, 1000, minE_NaI, maxE_NaI);
		  hE_NaI_LXe_2 	->SetMarkerColor(2);
	TH2F *hE_NaI_LXe_3	= new TH2F("hE_NaI_LXe_3", "", 1000, minE_LXe, maxE_LXe, 1000, minE_NaI, maxE_NaI);
		  hE_NaI_LXe_3	->SetMarkerColor(3);
	TH2F *hE_NaI_LXe_4	= new TH2F("hE_NaI_LXe_4", "", 1000, minE_LXe, maxE_LXe, 1000, minE_NaI, maxE_NaI);
		  hE_NaI_LXe_4	->SetMarkerColor(4);


	////////////////////////////////////////////////////////////////
	// PLOTS
	////////////////////////////////////////////////////////////////
	
	// energy in NaI VS time of flight
	TCanvas *cETOF_NaI = new TCanvas("cETOF_NaI", "cETOF_NaI", 0, 0, 700, 700);
			 cETOF_NaI->SetFillColor(10);
	/*	t2	->Draw("NaI_Edep:TOF>>hETOF_NaI_1", "LXe_ns==1", "");
		t2	->Draw("NaI_Edep:TOF>>hETOF_NaI_2", "LXe_ns==1 && hit[0]==0", "same");
		t2	->Draw("NaI_Edep:TOF>>hETOF_NaI_3", "LXe_ns==1 && hit[0]==0 && hitBefore[0]==0", "same");
		t2	->Draw("NaI_Edep:TOF>>hETOF_NaI_4", "LXe_ns==1 && hit[0]==0 && hitBefore[0]==0 && hitAfter[0]==0", "same");
	*/
		t2	->Draw("NaI_Edep:TOF>>hETOF_NaI_1", "LXe_ns==1 && hit[0]==0 && hitBefore[0]==0 && hitAfter[0]==0", "");
		t2	->Draw("NaI_Edep:TOF>>hETOF_NaI_2", "LXe_Edep>40 && LXe_ns==1 && hit[0]==0 && hitBefore[0]==0 && hitAfter[0]==0", "same");
		hETOF_NaI_1	->GetXaxis()->SetTitle("time of flight [ns]");
		hETOF_NaI_1	->GetYaxis()->SetTitle("Energy deposited in NaI [keV]");
		hETOF_NaI_1	->GetXaxis()->CenterTitle(true);
		hETOF_NaI_1	->GetYaxis()->CenterTitle(true);
		hETOF_NaI_1	->GetXaxis()->SetTitleOffset(1.20);
		hETOF_NaI_1	->GetYaxis()->SetTitleOffset(1.25);
		hETOF_NaI_1	->GetXaxis()->SetLabelSize(0.03);
		hETOF_NaI_1	->GetYaxis()->SetLabelSize(0.03);
/*	TLegend *legETOF_NaI = new TLegend(0.2974138,0.7366071,0.9971264,0.8839286,NULL,"brNDC");
 			 legETOF_NaI->SetBorderSize(0);
			 legETOF_NaI->SetTextFont(62);
			 legETOF_NaI->SetFillStyle(0);
	TLegendEntry *entryETOF_NaI=legETOF_NaI->AddEntry("hETOF_NaI_all","all coincident events","");
				  entryETOF_NaI->SetTextColor(1);
  				  entryETOF_NaI=legETOF_NaI->AddEntry("hETOF_NaI_true","true coincidence","");
				  entryETOF_NaI->SetTextColor(2);
	legETOF_NaI->Draw();
	l1->Draw("same");
	l2->Draw("same");
	l3->Draw("same");
	l4->Draw("same");	
*/	cETOF_NaI->Update();
	//char c1_tofE_png[256];
	//sprintf(c1_tofE_png, "%s/%s_tofE.png", OutFolder, FileName);
	//c1_tofE->SaveAs(c1_tofE_png);
	

	// time of flight
	TCanvas *c1_tof = new TCanvas("c1_tof", "c1_tof", 0, 0, 1000, 500);
			 c1_tof->SetFillColor(10);
			 c1_tof->SetLogy();
		//t2	->Draw("TOF>>hTOF_1", cutTrue1, "");
		//t2	->Draw("TOF>>hTOF_2", cutTrue2, "same");
		t2	->Draw("LXe_tmMin*1e9>>hTOF_1", cutTrue1, "");
		t2	->Draw("LXe_tmMin*1e9>>hTOF_2", cutTrue2, "same");
		hTOF_1	->GetXaxis()->SetTitle("Time of flight [ns]");
		hTOF_1	->GetYaxis()->SetTitle("Entries");
		hTOF_1	->GetXaxis()->CenterTitle(true);
		hTOF_1	->GetYaxis()->CenterTitle(true);
		hTOF_1	->GetXaxis()->SetLabelSize(0.03);
		hTOF_1	->GetYaxis()->SetLabelSize(0.03);
	c1_tof->Update();


/*	// energy deposited in NaI
	TCanvas *cE_NaI = new TCanvas("cE_NaI", "cE_NaI", 0, 0, 1000, 500);
			 cE_NaI->SetFillColor(10);
			 //c1_E->SetLogy();
		t1	->Draw("NaI_etot>>hE_NaI_all", "", "");
		t1	->Draw("NaI_etot>>hE_NaI_true", cutTrue, "same");
		hE_NaI_all	->GetXaxis()->SetTitle("Energy deposited in NaI [keV]");
		hE_NaI_all	->GetYaxis()->SetTitle("Entries");
		hE_NaI_all	->GetXaxis()->CenterTitle(true);
		hE_NaI_all	->GetYaxis()->CenterTitle(true);
		hE_NaI_all	->GetXaxis()->SetLabelSize(0.03);
		hE_NaI_all	->GetYaxis()->SetLabelSize(0.03);
	TLegend *legE_NaI = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
 			 legE_NaI->SetBorderSize(0);
			 legE_NaI->SetTextFont(62);
			 legE_NaI->SetFillStyle(0);
	TLegendEntry *entryE_NaI=legE_NaI->AddEntry("hE_NaI_all","all coincident events","");
				  entryE_NaI->SetTextColor(1);
  				  entryE_NaI=legE_NaI->AddEntry("hE_NaI_true","true coincidence","");
				  entryE_NaI->SetTextColor(2);
	legE_NaI->Draw();
	cE_NaI->Update();
	//char c1_E_png[256];
	//sprintf(c1_E_png, "%s/%s_E_NaI.png", OutFolder, FileName);
	//c1_E->SaveAs(c1_E_png);
*/

	// energy deposited in LXe
	TCanvas *cE_LXe = new TCanvas("cE_LXe", "cE_LXe", 0, 0, 1000, 500);
			 cE_LXe->SetFillColor(10);
			 cE_LXe->SetLogy();
		t2	->Draw("LXe_Edep>>hE_LXe_1", "LXe_ns==1", "");
		t2	->Draw("LXe_Edep>>hE_LXe_2", "LXe_ns==1 && hit[0]==0", "same");
		t2	->Draw("LXe_Edep>>hE_LXe_3", "LXe_ns==1 && hit[0]==0 && hitBefore[0]==0", "same");
		t2	->Draw("LXe_Edep>>hE_LXe_4", "LXe_ns==1 && hit[0]==0 && hitBefore[0]==0 && hitAfter[0]==0", "same");
		hE_LXe_1	->GetXaxis()->SetTitle("Energy deposited in LXe [keV]");
		hE_LXe_1	->GetYaxis()->SetTitle("Entries");
		hE_LXe_1	->GetXaxis()->CenterTitle(true);
		hE_LXe_1	->GetYaxis()->CenterTitle(true);
		hE_LXe_1	->GetXaxis()->SetLabelSize(0.03);
		hE_LXe_1	->GetYaxis()->SetLabelSize(0.03);
/*	TLegend *legE_LXe = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
 			 legE_LXe->SetBorderSize(0);
			 legE_LXe->SetTextFont(62);
			 legE_LXe->SetFillStyle(0);
	TLegendEntry *entryE_LXe=legE_LXe->AddEntry("hE_LXe_all","all coincident events","");
				  entryE_LXe->SetTextColor(1);
  				  entryE_LXe=legE_LXe->AddEntry("hE_LXe_true","true coincidence","");
				  entryE_LXe->SetTextColor(2);
  				  //entryE_LXe=legE_LXe->AddEntry("hE_LXe_trueCut","true coincidence with a cut","");
				  //entryE_LXe->SetTextColor(4);
	legE_LXe->Draw();
*/	cE_LXe->Update();
	//char c2_E_png[256];
	//sprintf(c2_E_png, "%s/%s_E_LXe.png", OutFolder, FileName);
	//c2_E->SaveAs(c2_E_png);


	// energy in NaI VS energy in LXe 
	TCanvas *cE_NaI_LXe = new TCanvas("cE_NaI_LXe", "cETOF_NaI", 0, 0, 700, 700);
			 cE_NaI_LXe->SetFillColor(10);
		t2	->Draw("NaI_Edep:LXe_Edep>>hE_NaI_LXe_1", "LXe_ns==1", "");
		t2	->Draw("NaI_Edep:LXe_Edep>>hE_NaI_LXe_2", "LXe_ns==1 && hit[0]==0", "same");
		t2	->Draw("NaI_Edep:LXe_Edep>>hE_NaI_LXe_3", "LXe_ns==1 && hit[0]==0 && hitBefore[0]==0", "same");
		t2	->Draw("NaI_Edep:LXe_Edep>>hE_NaI_LXe_4", "LXe_ns==1 && hit[0]==0 && hitBefore[0]==0 && hitAfter[0]==0", "same");
	
		//t2	->Draw("NaI_Edep:LXe_Edep>>hE_NaI_LXe_4", "LXe_ns==1 && hit[0]==0 && hitBefore[0]==0 && hitAfter[0]==0", "");
		//t2	->Draw("NaI_Edep:LXe_Edep>>hE_NaI_LXe_4", "LXe_Edep>40 && LXe_ns==1 && hit[0]==0 && hitBefore[0]==0 && hitAfter[0]==0", "same");
		hE_NaI_LXe_1	->GetXaxis()->SetTitle("Energy deposited in LXe [keV]");
		hE_NaI_LXe_1	->GetYaxis()->SetTitle("Energy deposited in NaI [keV]");
		hE_NaI_LXe_1	->GetXaxis()->CenterTitle(true);
		hE_NaI_LXe_1	->GetYaxis()->CenterTitle(true);
		hE_NaI_LXe_1	->GetXaxis()->SetTitleOffset(1.25);
		hE_NaI_LXe_1	->GetYaxis()->SetTitleOffset(1.25);
		hE_NaI_LXe_1	->GetXaxis()->SetLabelSize(0.03);
		hE_NaI_LXe_1	->GetYaxis()->SetLabelSize(0.03);
/*	TLegend *legE_NaI_LXe = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
 			 legE_NaI_LXe->SetBorderSize(0);
			 legE_NaI_LXe->SetTextFont(62);
			 legE_NaI_LXe->SetFillStyle(0);
	TLegendEntry *entryE_NaI_LXe=legE_NaI_LXe->AddEntry("hE_NaI_LXe_all","all coincident events","");
				  entryE_NaI_LXe->SetTextColor(1);
  				  entryE_NaI_LXe=legE_NaI_LXe->AddEntry("hE_NaI_LXe_true","true coincidence","");
				  entryE_NaI_LXe->SetTextColor(2);
  				  entryE_NaI_LXe=legE_NaI_LXe->AddEntry("hE_NaI_LXe_trueCut","true coincidence with a cut","");
				  entryE_NaI_LXe->SetTextColor(4);
	legE_NaI_LXe->Draw();
*/	cE_NaI_LXe->Update();
	//char c1_EE_png[256];
	//sprintf(c1_EE_png, "%s/%s_EE.png", OutFolder, FileName);
	//c1_EE->SaveAs(c1_EE_png);
	

	int nbinsEkin 	= 700;
	double minEkin = 0;
	double maxEkin = 700;

	TH1F *hEkin_1 	= new TH1F("hEkin_1", "", nbinsEkin, minEkin, maxEkin);
		  hEkin_1 	->SetLineColor(1);
		  hEkin_1 	->SetLineWidth(2);
	TH1F *hEkin_2 	= new TH1F("hEkin_2", "", nbinsEkin, minEkin, maxEkin);
		  hEkin_2 	->SetLineColor(2);
		  hEkin_2 	->SetLineWidth(2);
	
	// initial kinetic energy of the gamma
	TCanvas *cEkin = new TCanvas("cEkin", "cEkin", 0, 0, 700, 700);
			 cEkin->SetFillColor(10);
			 cEkin->SetLogy();
		t2	->Draw("LXe_EkinMax>>hEkin_1", cutTrue1, "");
		t2	->Draw("LXe_EkinMax>>hEkin_2", cutTrue2, "same");
		hEkin_1	->GetXaxis()->SetTitle("Kinetic energy of the gamma [keV]");
		hEkin_1	->GetYaxis()->SetTitle("Entries");
		hEkin_1	->GetXaxis()->CenterTitle(true);
		hEkin_1	->GetYaxis()->CenterTitle(true);
		hEkin_1	->GetXaxis()->SetTitleOffset(1.25);
		hEkin_1	->GetYaxis()->SetTitleOffset(1.25);
		hEkin_1	->GetXaxis()->SetLabelSize(0.03);
		hEkin_1	->GetYaxis()->SetLabelSize(0.03);
	cEkin->Update();


/*	// position of interactions and primaries
	int nbinsY 	= 200;
	double minY = 0;
	double maxY = 50;

	TH1F *hY_1 	= new TH1F("hY_1", "", nbinsY, minY, maxY);
		  hY_1 	->SetLineColor(1);
		  hY_1 	->SetLineWidth(2);
	TH1F *hY_2 	= new TH1F("hY_2", "", nbinsY, minY, maxY);
		  hY_2 	->SetLineColor(2);
		  hY_2 	->SetLineWidth(2);

	TH1F *hR_1 	= new TH1F("hR_1", "", nbinsY, minY, maxY);
		  hR_1 	->SetLineColor(1);
		  hR_1 	->SetLineWidth(2);
	TH1F *hR_2 	= new TH1F("hR_2", "", nbinsY, minY, maxY);
		  hR_2 	->SetLineColor(2);
		  hR_2 	->SetLineWidth(2);
	
	TCanvas *cY = new TCanvas("cY", "cY", 0, 0, 700, 700);
			 cY->SetFillColor(10);
			 cY->SetLogy();
		//t2	->Draw("LXe_y>>hY_1", cutTrue1, "");
		//t2	->Draw("LXe_y>>hY_2", cutTrue2, "same");
		t2	->Draw("primary_x>>hY_1", cutTrue1, "");
		t2	->Draw("primary_x>>hY_2", cutTrue2, "same");
		hY_1	->GetXaxis()->SetTitle("Y [mm]");
		hY_1	->GetYaxis()->SetTitle("Entries");
		hY_1	->GetXaxis()->CenterTitle(true);
		hY_1	->GetYaxis()->CenterTitle(true);
		hY_1	->GetXaxis()->SetTitleOffset(1.25);
		hY_1	->GetYaxis()->SetTitleOffset(1.25);
		hY_1	->GetXaxis()->SetLabelSize(0.03);
		hY_1	->GetYaxis()->SetLabelSize(0.03);
	cY->Update();

	TCanvas *cR = new TCanvas("cR", "cR", 0, 0, 700, 700);
			 cR->SetFillColor(10);
			 cR->SetLogy();
		//t2	->Draw("sqrt(LXe_y*LXe_y+LXe_z*LXe_z)>>hY_1", cutTrue1, "");
		//t2	->Draw("sqrt(LXe_y*LXe_y+LXe_z*LXe_z)>>hY_2", cutTrue2, "same");
		t2	->Draw("sqrt(primary_x*primary_x+primary_z*primary_z)>>hY_1", cutTrue1, "");
		t2	->Draw("sqrt(primary_x*primary_x+primary_z*primary_z)>>hY_2", cutTrue2, "same");
		hR_1	->GetXaxis()->SetTitle("R [mm]");
		hR_1	->GetYaxis()->SetTitle("Entries");
		hR_1	->GetXaxis()->CenterTitle(true);
		hR_1	->GetYaxis()->CenterTitle(true);
		hR_1	->GetXaxis()->SetTitleOffset(1.25);
		hR_1	->GetYaxis()->SetTitleOffset(1.25);
		hR_1	->GetXaxis()->SetLabelSize(0.03);
		hR_1	->GetYaxis()->SetLabelSize(0.03);
	cR->Update();
*/
/*	// energy deposited in NaI, all events that do not deposit energy in the LXe
	TCanvas *c1_Eall = new TCanvas("c1_Eall", "c1_Eall", 0, 0, 1000, 500);
			 c1_Eall->SetFillColor(10);
			 c1_Eall->SetLogy();
		t1	->Draw("NaI_etot>>h1_E1all_1", cutNo_LXe &&"NaI_etot>0", "");
		h1_E1all_1	->GetXaxis()->SetTitle("Energy deposited in NaI [keV]");
		h1_E1all_1	->GetYaxis()->SetTitle("Entries");
		h1_E1all_1	->GetXaxis()->CenterTitle(true);
		h1_E1all_1	->GetYaxis()->CenterTitle(true);
		h1_E1all_1	->GetXaxis()->SetLabelSize(0.03);
		h1_E1all_1	->GetYaxis()->SetLabelSize(0.03);
	c1_Eall->Update();
	//char c1_Eall_png[256];
	//sprintf(c1_Eall_png, "%s/%s_EallNaI.png", OutFolder, FileName);
	//c1_Eall->SaveAs(c1_Eall_png);
*/


/*	////////////////////////////////////////////////////////////////
	// SCAN
	////////////////////////////////////////////////////////////////
	char LogFile[500];
	sprintf(LogFile, "./TEST.txt");
	((TTreePlayer*)(t1->GetPlayer()))->SetScanRedirect(true);
	((TTreePlayer*)(t1->GetPlayer()))->SetScanFileName(LogFile);
	t1->SetScanField(0);
	//t1->Scan("eventid:nsteps:ed:time:etot:type:edproc:xp:yp:zp:NaI_nsteps:NaI_ed:NaI_time:NaI_etot:Teflon_ed:Teflon_time:Teflon_etot:DeadXe_ed:DeadXe_time:DeadXe_etot:GXe_ed:DeadXe_time:GXe_etot:Steel_ed:Steel_time:Steel_etot:Al_ed:Al_time:Al_etot:Lead_ed:Lead_time:Lead_etot:Rest_ed:Rest_time:Rest_etot", "", "");
	//t1->Scan("eventid:nsteps:etot", "", "");
	t1->Scan("eventid:nsteps:ed:etot:pre_kinE:post_kinE:type:parentid:trackid:NaI_etot:zp:sqrt(xp*xp+yp*yp):(Teflon_etot+Steel_etot+DeadXe_etot+GXe_etot+Lead_etot+Al_etot)","etot>200 && NaI_etot>0 && ed>0 && parentid==0 && (NaI_time-time)>0");
*/


/*	// save histograms
	char OutFile_char[256];
	sprintf(OutFile_char, "%s/%s.root", OutFolder, FileName);
	TFile *OutFile = new TFile(OutFile_char,"RECREATE");

	h1_E1all_1	->Write();
	h1_E1_1		->Write();
	h1_E1_2		->Write();
	h1_E2_1		->Write();
	h1_E2_2		->Write();
	h2_EE_1		->Write();
	h2_EE_2		->Write();
	h1_tof_1	->Write();
	h1_tof_2	->Write();
	h2_tofE_1	->Write();
	h2_tofE_2	->Write();

	OutFile->Close();
*/

}
