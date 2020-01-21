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

#include "vector"
#include "string"
#include "iostream"
#include "fstream"
#include "stdio.h"
#include "stdlib.h"


void plotCompton()
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
	const char *DataFolder 	= "/Users/alexkish/data/xuerich/mc/compton";	// idark05x
	const char *OutFolder 	= "./plots";
	const char *FileName 	= "10deg_withLeadPipe75cm_Aperture";
	
	Char_t  DataFile1[500];
	Char_t  DataFile2[500];
	Char_t  DataFile3[500];
	Char_t  DataFile4[500];
	Char_t  DataFile5[500];
	Char_t  DataFile6[500];
	
/*	// 43 deg
	sprintf(DataFile1, "%s/Xu_Compton_43deg_1.6e9.root", DataFolder);
	sprintf(DataFile2, "%s/Xu_Compton_43deg_1.6e9_1.root", DataFolder);
	sprintf(DataFile3, "%s/Xu_Compton_43deg_1.6e9_2.root", DataFolder);
	TChain *t1	= new TChain("t1");
			t1	->Add(DataFile1);
			t1	->Add(DataFile2);
			t1	->Add(DataFile3);
*/
/*	// 34 deg, standard
	sprintf(DataFile1, "%s/Xu_Compton_34deg_3.2e9.root", DataFolder);
	sprintf(DataFile2, "%s/Xu_Compton_34deg_3.2e9_1.root", DataFolder);
	sprintf(DataFile3, "%s/Xu_Compton_34deg_3.2e9_2.root", DataFolder);
	sprintf(DataFile4, "%s/Xu_Compton_34deg_3.2e9_3.root", DataFolder);
	sprintf(DataFile5, "%s/Xu_Compton_34deg_3.2e9_4.root", DataFolder);
	sprintf(DataFile6, "%s/Xu_Compton_34deg_3.2e9_5.root", DataFolder);
	TChain *t1	= new TChain("t1");
			t1	->Add(DataFile1);
			t1	->Add(DataFile2);
			t1	->Add(DataFile3);
			t1	->Add(DataFile4);
			t1	->Add(DataFile5);
			t1	->Add(DataFile6);
*/
/*	// 34 deg, extended shield for NaI (50cm)
	sprintf(DataFile1, "%s/Xu_Compton_34deg_withLeadPipe_3.2e9.root", DataFolder);
	sprintf(DataFile2, "%s/Xu_Compton_34deg_withLeadPipe_3.2e9_1.root", DataFolder);
	sprintf(DataFile3, "%s/Xu_Compton_34deg_withLeadPipe_3.2e9_2.root", DataFolder);
	sprintf(DataFile4, "%s/Xu_Compton_34deg_withLeadPipe_3.2e9_3.root", DataFolder);
	sprintf(DataFile5, "%s/Xu_Compton_34deg_withLeadPipe_3.2e9_4.root", DataFolder);
	TChain *t1	= new TChain("t1");
			t1	->Add(DataFile1);
			t1	->Add(DataFile2);
			t1	->Add(DataFile3);
			t1	->Add(DataFile4);
			t1	->Add(DataFile5);
*/


/*	// 10 deg, standard
	sprintf(DataFile1, "%s/Xu_Compton_10deg_3.12e9.root", DataFolder);
	sprintf(DataFile2, "%s/Xu_Compton_10deg_3.12e9_1.root", DataFolder);
	sprintf(DataFile3, "%s/Xu_Compton_10deg_3.12e9_2.root", DataFolder);
	sprintf(DataFile4, "%s/Xu_Compton_10deg_3.12e9_3.root", DataFolder);
	sprintf(DataFile5, "%s/Xu_Compton_10deg_3.12e9_4.root", DataFolder);
	sprintf(DataFile6, "%s/Xu_Compton_10deg_3.12e9_5.root", DataFolder);
	TChain *t1	= new TChain("t1");
			t1	->Add(DataFile1);
			t1	->Add(DataFile2);
			t1	->Add(DataFile3);
			t1	->Add(DataFile4);
			t1	->Add(DataFile5);
			t1	->Add(DataFile6);
*/
/*	// 10 deg, lead pipe 50cm
	sprintf(DataFile1, "%s/Xu_Compton_10deg_Pipe50cm_3.12e9.root", DataFolder);
	sprintf(DataFile2, "%s/Xu_Compton_10deg_Pipe50cm_3.12e9_1.root", DataFolder);
	sprintf(DataFile3, "%s/Xu_Compton_10deg_Pipe50cm_3.12e9_2.root", DataFolder);
	sprintf(DataFile4, "%s/Xu_Compton_10deg_Pipe50cm_3.12e9_3.root", DataFolder);
	sprintf(DataFile5, "%s/Xu_Compton_10deg_Pipe50cm_3.12e9_4.root", DataFolder);
	sprintf(DataFile6, "%s/Xu_Compton_10deg_Pipe50cm_3.12e9_5.root", DataFolder);
	TChain *t1	= new TChain("t1");
			t1	->Add(DataFile1);
			t1	->Add(DataFile2);
			t1	->Add(DataFile3);
			t1	->Add(DataFile4);
			t1	->Add(DataFile5);
			t1	->Add(DataFile6);
*/
/*	// 10 deg, lead pipe 75cm
	sprintf(DataFile1, "%s/Xu_Compton_10deg_Pipe75cm_3.12e9.root", DataFolder);
	sprintf(DataFile2, "%s/Xu_Compton_10deg_Pipe75cm_3.12e9_1.root", DataFolder);
	sprintf(DataFile3, "%s/Xu_Compton_10deg_Pipe75cm_3.12e9_2.root", DataFolder);
	sprintf(DataFile4, "%s/Xu_Compton_10deg_Pipe75cm_3.12e9_3.root", DataFolder);
	sprintf(DataFile5, "%s/Xu_Compton_10deg_Pipe75cm_3.12e9_4.root", DataFolder);
	TChain *t1	= new TChain("t1");
			t1	->Add(DataFile1);
			t1	->Add(DataFile2);
			t1	->Add(DataFile3);
			t1	->Add(DataFile4);
			t1	->Add(DataFile5);
*/
	// 10 deg, lead pipe 75cm, aperture 5cm
	sprintf(DataFile1, "%s/Xu_Compton_10deg_Pipe75cm_Aperture_3.12e9.root", DataFolder);
	sprintf(DataFile2, "%s/Xu_Compton_10deg_Pipe75cm_Aperture_3.12e9_1.root", DataFolder);
	sprintf(DataFile3, "%s/Xu_Compton_10deg_Pipe75cm_Aperture_3.12e9_2.root", DataFolder);
	sprintf(DataFile4, "%s/Xu_Compton_10deg_Pipe75cm_Aperture_3.12e9_3.root", DataFolder);
	sprintf(DataFile5, "%s/Xu_Compton_10deg_Pipe75cm_Aperture_3.12e9_4.root", DataFolder);
	TChain *t1	= new TChain("t1");
			t1	->Add(DataFile1);
			t1	->Add(DataFile2);
			t1	->Add(DataFile3);
			t1	->Add(DataFile4);
			t1	->Add(DataFile5);



	////////////////////////////////////////////////////////////////
	// ALIASES, CUTS
	t1->SetAlias("TOF", "(Alt$(NaI_time,0)-Alt$(time,0))*1e9");

	TCut cut_coincidence 	= "etot>0 && NaI_etot>0";

	TCut cut_true_tof		= "Alt$(NaI_time,0)-Alt$(time,0)>0";

	//TCut cutNo_LXe		= "etot==0 || (etot>0 && (Alt$(NaI_time,0)-Alt$(time,0))>0)";
	TCut cutNo_LXe			= "etot==0";
	TCut cutNo_NaI			= "NaI_etot==0 || (NaI_etot>0 && (Alt$(NaI_time,0)-Alt$(time,0))>0)";

/*	TCut cutNo_DeadXe		= "DeadXe_etot==0 || (DeadXe_etot>0 && (Alt$(DeadXe_time,0)-Alt$(NaI_time,0))>0)";
	TCut cutNo_GXe			= "GXe_etot==0 || (GXe_etot>0 && (Alt$(GXe_time,0)-Alt$(NaI_time,0))>0)";
	TCut cutNo_Teflon		= "Teflon_etot==0 || (Teflon_etot>0 && (Alt$(Teflon_time,0)-Alt$(NaI_time,0))>0)";
	TCut cutNo_Steel		= "Steel_etot==0 || (Steel_etot>0 && (Alt$(Steel_time,0)-Alt$(NaI_time,0))>0)";
	TCut cutNo_Al			= "Al_etot==0 || (Al_etot>0 && (Alt$(Al_time,0)-Alt$(NaI_time,0))>0)";
	TCut cutNo_Lead			= "Lead_etot==0 || (Lead_etot>0 && (Alt$(Lead_time,0)-Alt$(NaI_time,0))>0)";
	TCut cutNo_Rest			= "Rest_etot==0 || (Rest_etot>0 && (Alt$(Rest_time,0)-Alt$(NaI_time,0))>0)";
	TCut cutTrue			= cutNo_DeadXe && cutNo_GXe && cutNo_Teflon && cutNo_Steel && cutNo_Al && cutNo_Lead && cutNo_Rest;
*/
	TCut cutNo_DeadXe		= "(Alt$(DeadXe_etot,0))==0";
	TCut cutNo_GXe			= "(Alt$(GXe_etot,0))==0";
	TCut cutNo_Teflon		= "(Alt$(Teflon_etot,0))==0";
	TCut cutNo_Steel		= "(Alt$(Steel_etot,0))==0";
	TCut cutNo_Al			= "(Alt$(Al_etot,0))==0";
	TCut cutNo_Lead			= "(Alt$(Lead_etot,0))==0";
	TCut cutNo_Rest			= "(Alt$(Rest_etot,0))==0";
	TCut cutTrue			= cutNo_DeadXe && cutNo_GXe && cutNo_Teflon && cutNo_Steel && cutNo_Al && cutNo_Lead && cutNo_Rest;

	TCut cutYs_DeadXe		= "DeadXe_etot>0 && (Alt$(DeadXe_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_GXe			= "GXe_etot>0 && (Alt$(GXe_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Teflon		= "Teflon_etot>0 && (Alt$(Teflon_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Steel		= "Steel_etot>0 && (Alt$(Steel_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Al			= "Al_etot>0 && (Alt$(Al_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Lead			= "Lead_etot>0 && (Alt$(Lead_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Rest			= "Rest_etot>0 && (Alt$(Rest_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutFalse			= cutYs_DeadXe || cutYs_GXe || cutYs_Teflon || cutYs_Steel || cutYs_Al || cutYs_Lead || cutYs_Rest;

	
	////////////////////////////////////////////////////////////////
/*	// SCAN
	char LogFile[500];
	sprintf(LogFile, "./plots/EventswithHighEtot.txt");
	((TTreePlayer*)(t1->GetPlayer()))->SetScanRedirect(true);
	((TTreePlayer*)(t1->GetPlayer()))->SetScanFileName(LogFile);
	t1->SetScanField(0);
	t1->Scan("eventid:nsteps:ed:time:etot:type:edproc:xp:yp:zp:NaI_nsteps:NaI_ed:NaI_time:NaI_etot:Teflon_ed:Teflon_time:Teflon_etot:DeadXe_ed:DeadXe_time:DeadXe_etot:GXe_ed:DeadXe_time:GXe_etot:Steel_ed:Steel_time:Steel_etot:Al_ed:Al_time:Al_etot:Lead_ed:Lead_time:Lead_etot:Rest_ed:Rest_time:Rest_etot", cut_coincidence && cutTrue &&"etot>100", "");
*/
	
	////////////////////////////////////////////////////////////////
	// HISTOGRAMS and PLOTS
	int nbins_tof 	= 100;
	//double min_tof 	= 0;
	//double max_tof 	= 7;
	double min_tof 	= 2.5;
	double max_tof 	= 6.0;

	TH1F *h1_tof_1 = new TH1F("h1_tof_1", "", nbins_tof, min_tof, max_tof);
		  h1_tof_1 ->SetLineColor(4);
		  h1_tof_1 ->SetLineWidth(2);
	TH1F *h1_tof_2 = new TH1F("h1_tof_2", "", nbins_tof, min_tof, max_tof);
		  h1_tof_2 ->SetLineColor(2);
		  h1_tof_2 ->SetLineWidth(2);

	int nbins_E 	= 700;
	double min_E 	= 0;
	double max_E 	= 700;

	TH1F *h1_E1_1 = new TH1F("h1_E1_1", "", nbins_E, min_E, max_E);
		  h1_E1_1 ->SetLineColor(1);
		  h1_E1_1 ->SetLineWidth(2);
	TH1F *h1_E1all_1 = new TH1F("h1_E1all_1", "", nbins_E, min_E, max_E);
		  h1_E1all_1 ->SetLineColor(1);
		  h1_E1all_1 ->SetLineWidth(2);
	TH1F *h1_E1_2 = new TH1F("h1_E1_2", "", nbins_E, min_E, max_E);
		  h1_E1_2 ->SetLineColor(2);
		  h1_E1_2 ->SetLineWidth(2);

	TH1F *h1_E2_1 = new TH1F("h1_E2_1", "", nbins_E, min_E, max_E);
		  h1_E2_1 ->SetLineColor(1);
		  h1_E2_1 ->SetLineWidth(2);
	TH1F *h1_E2_2 = new TH1F("h1_E2_2", "", nbins_E, min_E, max_E);
		  h1_E2_2 ->SetLineColor(2);
		  h1_E2_2 ->SetLineWidth(2);
	TH1F *h1_E2_3 = new TH1F("h1_E2_3", "", nbins_E, min_E, max_E);
		  h1_E2_3 ->SetLineColor(3);
		  h1_E2_3 ->SetLineWidth(2);

	TH2F *h2_tofE_1 = new TH2F("h2_tofE_1", "", nbins_tof, min_tof, max_tof, nbins_E, min_E, max_E);
		  h2_tofE_1 ->SetMarkerColor(1);
		  h2_tofE_1 ->SetLineColor(1);
		  h2_tofE_1 ->SetLineWidth(2);
	TH2F *h2_tofE_2 = new TH2F("h2_tofE_2", "", nbins_tof, min_tof, max_tof, nbins_E, min_E, max_E);
		  h2_tofE_2 ->SetMarkerColor(2);
		  h2_tofE_2 ->SetLineColor(2);
		  h2_tofE_2 ->SetLineWidth(2);
	TH2F *h2_tofE_3 = new TH2F("h2_tofE_3", "", nbins_tof, min_tof, max_tof, nbins_E, min_E, max_E);
		  h2_tofE_3 ->SetMarkerColor(3);
		  h2_tofE_3 ->SetLineColor(3);
		  h2_tofE_3 ->SetLineWidth(2);

	TH2F *h2_EE_1 = new TH2F("h2_EE_1", "", nbins_E, min_E, max_E, nbins_E, min_E, max_E);
		  h2_EE_1 ->SetMarkerColor(1);
		  h2_EE_1 ->SetLineColor(1);
		  h2_EE_1 ->SetLineWidth(2);
	TH2F *h2_EE_2 = new TH2F("h2_EE_2", "", nbins_E, min_E, max_E, nbins_E, min_E, max_E);
		  h2_EE_2 ->SetMarkerColor(2);
		  h2_EE_2 ->SetLineColor(2);
		  h2_EE_2 ->SetLineWidth(2);
	TH2F *h2_EE_3 = new TH2F("h2_EE_3", "", nbins_E, min_E, max_E, nbins_E, min_E, max_E);
		  h2_EE_3 ->SetMarkerColor(3);
		  h2_EE_3 ->SetLineColor(3);
		  h2_EE_3 ->SetLineWidth(2);

	

/*	// time of flight
	TCanvas *c1_tof = new TCanvas("c1_tof", "c1_tof", 0, 0, 1000, 500);
			 c1_tof->SetFillColor(10);
			 c1_tof->SetLogy();
		t1	->Draw("TOF>>h1_tof_1", cut_coincidence, "");
		t1	->Draw("TOF>>h1_tof_2", cut_coincidence && cutTrue, "same");
		h1_tof_1	->GetXaxis()->SetTitle("Time of flight [ns]");
		h1_tof_1	->GetYaxis()->SetTitle("Entries");
		h1_tof_1	->GetXaxis()->CenterTitle(true);
		h1_tof_1	->GetYaxis()->CenterTitle(true);
		h1_tof_1	->GetXaxis()->SetLabelSize(0.03);
		h1_tof_1	->GetYaxis()->SetLabelSize(0.03);
	TLegend *leg1_tof = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
 			 leg1_tof->SetBorderSize(0);
			 leg1_tof->SetTextFont(62);
			 leg1_tof->SetFillStyle(0);
	TLegendEntry *entry1_tof=leg1_tof->AddEntry("h1_tof_1","all coincident events","");
				  entry1_tof->SetTextColor(1);
  				  entry1_tof=leg1_tof->AddEntry("h1_tof_2","true coincidence","");
				  entry1_tof->SetTextColor(2);
	leg1_tof->Draw();
	c1_tof->Update();
	char c1_tof_png[256];
	sprintf(c1_tof_png, "%s/%s_TOF.png", OutFolder, FileName);
	c1_tof->SaveAs(c1_tof_png);


	// energy deposited in NaI
	TCanvas *c1_E = new TCanvas("c1_E", "c1_E", 0, 0, 1000, 500);
			 c1_E->SetFillColor(10);
			 //c1_E->SetLogy();
		t1	->Draw("NaI_etot>>h1_E1_1", cut_coincidence, "");
		t1	->Draw("NaI_etot>>h1_E1_2", cut_coincidence && cutTrue, "same");
		h1_E1_1	->GetXaxis()->SetTitle("Energy deposited in NaI [keV]");
		h1_E1_1	->GetYaxis()->SetTitle("Entries");
		h1_E1_1	->GetXaxis()->CenterTitle(true);
		h1_E1_1	->GetYaxis()->CenterTitle(true);
		h1_E1_1	->GetXaxis()->SetLabelSize(0.03);
		h1_E1_1	->GetYaxis()->SetLabelSize(0.03);
	TLegend *leg1_E = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
 			 leg1_E->SetBorderSize(0);
			 leg1_E->SetTextFont(62);
			 leg1_E->SetFillStyle(0);
	TLegendEntry *entry1_E=leg1_E->AddEntry("h1_E1_1","all coincident events","");
				  entry1_E->SetTextColor(1);
  				  entry1_E=leg1_E->AddEntry("h1_E1_2","true coincidence","");
				  entry1_E->SetTextColor(2);
	leg1_E->Draw();
	c1_E->Update();
	char c1_E_png[256];
	sprintf(c1_E_png, "%s/%s_E_NaI.png", OutFolder, FileName);
	c1_E->SaveAs(c1_E_png);
*/


	// energy deposited in LXe
	TCanvas *c2_E = new TCanvas("c2_E", "c2_E", 0, 0, 1000, 500);
			 c2_E->SetFillColor(10);
			 //c2_E->SetLogy();
		t1	->Draw("etot>>h1_E2_1", cut_coincidence, "");
		t1	->Draw("etot>>h1_E2_2", cut_coincidence && cutTrue, "same");
		t1	->Draw("etot>>h1_E2_3", cut_coincidence && cutTrue &&"etot>100", "same");
		h1_E2_1	->GetXaxis()->SetTitle("Energy deposited in LXe [keV]");
		h1_E2_1	->GetYaxis()->SetTitle("Entries");
		h1_E2_1	->GetXaxis()->CenterTitle(true);
		h1_E2_1	->GetYaxis()->CenterTitle(true);
		h1_E2_1	->GetXaxis()->SetLabelSize(0.03);
		h1_E2_1	->GetYaxis()->SetLabelSize(0.03);
	TLegend *leg2_E = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
 			 leg2_E->SetBorderSize(0);
			 leg2_E->SetTextFont(62);
			 leg2_E->SetFillStyle(0);
	TLegendEntry *entry2_E=leg2_E->AddEntry("h1_E2_1","all coincident events","");
				  entry2_E->SetTextColor(1);
  				  entry2_E=leg2_E->AddEntry("h1_E2_2","true coincidence","");
				  entry2_E->SetTextColor(2);
  				  entry2_E=leg2_E->AddEntry("h1_E2_3","true coincidence, E>100keV","");
				  entry2_E->SetTextColor(3);
	leg2_E->Draw();
	c2_E->Update();
	char c2_E_png[256];
	sprintf(c2_E_png, "%s/%s_E_LXe.png", OutFolder, FileName);
	//c2_E->SaveAs(c2_E_png);


	// energy in NaI VS energy in LXe 
	TCanvas *c1_EE = new TCanvas("c1_EE", "c1_EE", 0, 0, 700, 700);
			 c1_EE->SetFillColor(10);
		t1	->Draw("NaI_etot:etot>>h2_EE_1", cut_coincidence, "");
		t1	->Draw("NaI_etot:etot>>h2_EE_2", cut_coincidence && cutTrue, "same");
		t1	->Draw("NaI_etot:etot>>h2_EE_3", cut_coincidence && cutTrue &&"etot>100", "same");
		h2_EE_1	->GetXaxis()->SetTitle("Energy deposited in LXe [keV]");
		h2_EE_1	->GetYaxis()->SetTitle("Energy deposited in NaI [keV]");
		h2_EE_1	->GetXaxis()->CenterTitle(true);
		h2_EE_1	->GetYaxis()->CenterTitle(true);
		h2_EE_1	->GetXaxis()->SetTitleOffset(1.25);
		h2_EE_1	->GetYaxis()->SetTitleOffset(1.25);
		h2_EE_1	->GetXaxis()->SetLabelSize(0.03);
		h2_EE_1	->GetYaxis()->SetLabelSize(0.03);
	TLegend *leg1_EE = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
 			 leg1_EE->SetBorderSize(0);
			 leg1_EE->SetTextFont(62);
			 leg1_EE->SetFillStyle(0);
	TLegendEntry *entry1_EE=leg1_EE->AddEntry("h2_EE_1","all coincident events","");
				  entry1_EE->SetTextColor(1);
  				  entry1_EE=leg1_EE->AddEntry("h2_EE_2","true coincidence","");
				  entry1_EE->SetTextColor(2);
  				  entry1_EE=leg1_EE->AddEntry("h2_EE_3","true coincidence, E>100keV","");
				  entry1_EE->SetTextColor(3);
	leg1_EE->Draw();
	c1_EE->Update();
	char c1_EE_png[256];
	sprintf(c1_EE_png, "%s/%s_EE.png", OutFolder, FileName);
	//c1_EE->SaveAs(c1_EE_png);

	// energy in NaI VS time of flight
	TCanvas *c1_tofE = new TCanvas("c1_tofE", "c1_tofE", 0, 0, 700, 700);
			 c1_tofE->SetFillColor(10);
		t1	->Draw("NaI_etot:TOF>>h2_tofE_1", cut_coincidence, "");
		t1	->Draw("NaI_etot:TOF>>h2_tofE_2", cut_coincidence && cutTrue, "same");
		t1	->Draw("NaI_etot:TOF>>h2_tofE_3", cut_coincidence && cutTrue &&"etot>100", "same");
		h2_tofE_1	->GetXaxis()->SetTitle("time of flight [ns]");
		h2_tofE_1	->GetYaxis()->SetTitle("Energy deposited in NaI [keV]");
		h2_tofE_1	->GetXaxis()->CenterTitle(true);
		h2_tofE_1	->GetYaxis()->CenterTitle(true);
		h2_tofE_1	->GetXaxis()->SetTitleOffset(1.25);
		h2_tofE_1	->GetYaxis()->SetTitleOffset(1.25);
		h2_tofE_1	->GetXaxis()->SetLabelSize(0.03);
		h2_tofE_1	->GetYaxis()->SetLabelSize(0.03);
	TLegend *leg1_tofE = new TLegend(0.2974138,0.7366071,0.9971264,0.8839286,NULL,"brNDC");
 			 leg1_tofE->SetBorderSize(0);
			 leg1_tofE->SetTextFont(62);
			 leg1_tofE->SetFillStyle(0);
	TLegendEntry *entry1_tofE=leg1_tofE->AddEntry("h2_tofE_1","all coincident events","");
				  entry1_tofE->SetTextColor(1);
  				  entry1_tofE=leg1_tofE->AddEntry("h2_tofE_2","true coincidence","");
				  entry1_tofE->SetTextColor(2);
  				  entry1_tofE=leg1_tofE->AddEntry("h2_tofE_3","true coincidence, E>100keV","");
				  entry1_tofE->SetTextColor(3);
	leg1_tofE->Draw();
	c1_tofE->Update();
	char c1_tofE_png[256];
	sprintf(c1_tofE_png, "%s/%s_tofE.png", OutFolder, FileName);
	//c1_tofE->SaveAs(c1_tofE_png);


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
	char c1_Eall_png[256];
	sprintf(c1_Eall_png, "%s/%s_EallNaI.png", OutFolder, FileName);
	c1_Eall->SaveAs(c1_Eall_png);


	// save histograms
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
