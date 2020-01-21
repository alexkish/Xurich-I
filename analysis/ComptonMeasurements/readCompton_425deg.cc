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


void readCompton_425deg()
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
	
	//sprintf(DataFile1, "%s/Xu_Compton_4.25deg_5e9.root", DataFolder);
	//sprintf(DataFile2, "%s/Xu_Compton_4.25deg_5e9c.root", DataFolder);
	//sprintf(DataFile1, "%s/Xu_Compton_4.25deg_4e10.root", DataFolder);
	sprintf(DataFile1, "%s/Xu_Compton_6.25deg_4e10.root", DataFolder);
	TChain *t1	= new TChain("t1");
			t1	->Add(DataFile1);
			//t1	->Add(DataFile2);


	////////////////////////////////////////////////////////////////
	// ALIASES, CUTS
	t1->SetAlias("TOF", "(Alt$(NaI_time,0)-Alt$(time,0))*1e9");

	TCut cut_coincidence 	= "etot>0 && NaI_etot>0";

	TCut cut_true_tof		= "Alt$(NaI_time,0)-Alt$(time,0)>0";

	//TCut cutNo_LXe		= "etot==0 || (etot>0 && (Alt$(NaI_time,0)-Alt$(time,0))>0)";
	TCut cutNo_LXe			= "etot==0";
	TCut cutNo_NaI			= "NaI_etot==0 || (NaI_etot>0 && (Alt$(NaI_time,0)-Alt$(time,0))>0)";

	TCut cutNo_DeadXe		= "(Alt$(DeadXe_etot,0))==0";
	TCut cutNo_GXe			= "(Alt$(GXe_etot,0))==0";
	TCut cutNo_Teflon		= "(Alt$(Teflon_etot,0))==0";
	TCut cutNo_Steel		= "(Alt$(Steel_etot,0))==0";
	TCut cutNo_Al			= "(Alt$(Al_etot,0))==0";
	TCut cutNo_Lead			= "(Alt$(Lead_etot,0))==0";
	TCut cutNo_Rest			= "(Alt$(Rest_etot,0))==0";
	TCut cutTrue			= cutNo_DeadXe && cutNo_GXe && cutNo_Teflon && cutNo_Steel && cutNo_Al && cutNo_Lead && cutNo_Rest;

/*	TCut cutNo_DeadXe		= "DeadXe_etot==0 && (Alt$(NaI_time,0)-Alt$(DeadXe_time,0))>0";
	TCut cutNo_GXe			= "GXe_etot==0 && (Alt$(NaI_time,0)-Alt$(GXe_time,0))>0";
	TCut cutNo_Teflon		= "Teflon_etot==0 && (Alt$(NaI_time,0)-Alt$(Teflon_time,0))>0";
	TCut cutNo_Steel		= "Steel_etot==0 && (Alt$(NaI_time,0)-Alt$(Steel_time,0))>0";
	TCut cutNo_Al			= "Al_etot==0 && (Alt$(NaI_time,0)-Alt$(Al_time,0))>0";
	TCut cutNo_Lead			= "Lead_etot==0 && (Alt$(NaI_time,0)-Alt$(Lead_time,0))>0";
	TCut cutNo_Rest			= "Rest_etot==0 && (Alt$(NaI_time,0)-Alt$(Rest_time,0))>0";
	TCut cutTrue			= cutNo_DeadXe && cutNo_GXe && cutNo_Teflon && cutNo_Steel && cutNo_Al && cutNo_Lead && cutNo_Rest;
*/
/*	TCut cutNo_DeadXe		= "DeadXe_etot==0 && (Alt$(NaI_time,0)-Alt$(DeadXe_time,0))>0 && (Alt$(time,0)-Alt$(DeadXe_time,0))>0";
	TCut cutNo_GXe			= "GXe_etot==0 && (Alt$(NaI_time,0)-Alt$(GXe_time,0))>0 && (Alt$(time,0)-Alt$(GXe_time,0))>0";
	TCut cutNo_Teflon		= "Teflon_etot==0 && (Alt$(NaI_time,0)-Alt$(Teflon_time,0))>0 && (Alt$(time,0)-Alt$(Teflon_time,0))>0";
	TCut cutNo_Steel		= "Steel_etot==0 && (Alt$(NaI_time,0)-Alt$(Steel_time,0))>0 && (Alt$(time,0)-Alt$(Steel_time,0))>0";
	TCut cutNo_Al			= "Al_etot==0 && (Alt$(NaI_time,0)-Alt$(Al_time,0))>0 && (Alt$(time,0)-Alt$(Al_time,0))>0";
	TCut cutNo_Lead			= "Lead_etot==0 && (Alt$(NaI_time,0)-Alt$(Lead_time,0))>0 && (Alt$(time,0)-Alt$(Lead_time,0))>0";
	TCut cutNo_Rest			= "Rest_etot==0 && (Alt$(NaI_time,0)-Alt$(Rest_time,0))>0 && (Alt$(time,0)-Alt$(Rest_time,0))>0";
	TCut cutTrue			= cutNo_DeadXe && cutNo_GXe && cutNo_Teflon && cutNo_Steel && cutNo_Al && cutNo_Lead && cutNo_Rest;
*/
	TCut cutYs_DeadXe		= "DeadXe_etot>0 && (Alt$(DeadXe_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_GXe			= "GXe_etot>0 && (Alt$(GXe_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Teflon		= "Teflon_etot>0 && (Alt$(Teflon_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Steel		= "Steel_etot>0 && (Alt$(Steel_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Al			= "Al_etot>0 && (Alt$(Al_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Lead			= "Lead_etot>0 && (Alt$(Lead_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutYs_Rest			= "Rest_etot>0 && (Alt$(Rest_time,0)-Alt$(NaI_time,0))<0)";
	TCut cutFalse			= cutYs_DeadXe || cutYs_GXe || cutYs_Teflon || cutYs_Steel || cutYs_Al || cutYs_Lead || cutYs_Rest;

	
	////////////////////////////////////////////////////////////////
	// APPLY CUTS AND CREATE EVENT LIST
	////////////////////////////////////////////////////////////////
/*	//t1->Draw(">>ll", cut_coincidence);
	t1->Draw(">>ll", cut_coincidence && cutTrue);
	TEventList *list = (TEventList *)gDirectory->Get("ll");
	// SAVE LIST
	TFile *efile = new TFile("./list/EventList_4.25deg_4e10_CoincidenceTrue.root","recreate");
	list->Write();
	efile->Close();
*/
/*	////////////////////////////////////////////////////////////////
	// LOAD EVENT LIST WITH CUTS
	////////////////////////////////////////////////////////////////
	TFile *efile = new TFile("./list/EventList_4.25deg_4e10_CoincidenceAll.root","READ");
	//TFile *efile = new TFile("./list/EventList_4.25deg_1e10_CoincidenceTrue.root","READ");
	//TFile *efile = new TFile("./list/EventList_4.25deg_1e10_CoincidenceAll.root","READ");
	TEventList *list = (TEventList*)efile->Get("ll");
	t1->SetEventList(list);
*/

/*	////////////////////////////////////////////////////////////////
	// SCAN
	////////////////////////////////////////////////////////////////
	char LogFile[500];
	sprintf(LogFile, "./list/EventsShort_CoincidenceTrue.txt");
	((TTreePlayer*)(t1->GetPlayer()))->SetScanRedirect(true);
	((TTreePlayer*)(t1->GetPlayer()))->SetScanFileName(LogFile);
	t1->SetScanField(0);
	//t1->Scan("eventid:nsteps:ed:time:etot:type:edproc:xp:yp:zp:NaI_nsteps:NaI_ed:NaI_time:NaI_etot:Teflon_ed:Teflon_time:Teflon_etot:DeadXe_ed:DeadXe_time:DeadXe_etot:GXe_ed:DeadXe_time:GXe_etot:Steel_ed:Steel_time:Steel_etot:Al_ed:Al_time:Al_etot:Lead_ed:Lead_time:Lead_etot:Rest_ed:Rest_time:Rest_etot", "", "");
	t1->Scan("eventid:nsteps:etot", "", "");
*/
	
	////////////////////////////////////////////////////////////////
	// HISTOGRAMS
	////////////////////////////////////////////////////////////////
	int nbins_tof 	= 200;
	//double min_tof = 0;
	//double max_tof = 7;
	double min_tof 	= 3.0;
	double max_tof 	= 6.0;

	TH1F *hTOF_all 	= new TH1F("hTOF_all", "", nbins_tof, min_tof, max_tof);
		  hTOF_all 	->SetLineColor(1);
		  hTOF_all 	->SetLineWidth(2);
	TH1F *hTOF_true = new TH1F("hTOF_true", "", nbins_tof, min_tof, max_tof);
		  hTOF_true ->SetLineColor(2);
		  hTOF_true ->SetLineWidth(2);
	TH1F *hTOF_trueCut = new TH1F("hTOF_trueCut", "", nbins_tof, min_tof, max_tof);
		  hTOF_trueCut ->SetLineColor(4);
		  hTOF_trueCut ->SetLineWidth(2);

	int nbinsE_NaI 	= 1400;
	double minE_NaI = 0;
	double maxE_NaI = 700;

	TH1F *hE_NaI_all 	= new TH1F("hE_NaI_all", "", nbinsE_NaI, minE_NaI, maxE_NaI);
		  hE_NaI_all 	->SetLineColor(1);
		  hE_NaI_all 	->SetLineWidth(2);
	TH1F *hE_NaI_true 	= new TH1F("hE_NaI_true", "", nbinsE_NaI, minE_NaI, maxE_NaI);
		  hE_NaI_true 	->SetLineColor(2);
		  hE_NaI_true 	->SetLineWidth(2);

	TH2F *hETOF_NaI_all 	= new TH2F("hETOF_NaI_all", "", 1000, min_tof, max_tof, 1000, minE_NaI, maxE_NaI);
		  hETOF_NaI_all 	->SetMarkerColor(1);
	TH2F *hETOF_NaI_true 	= new TH2F("hETOF_NaI_true", "", 1000, min_tof, max_tof, 1000, minE_NaI, maxE_NaI);
		  hETOF_NaI_true 	->SetMarkerColor(2);

	int nbinsE_LXe 	= 1200;
	double minE_LXe = 0;
	double maxE_LXe = 700;

	TH1F *hE_LXe_all 		= new TH1F("hE_LXe_all", "", nbinsE_LXe, minE_LXe, maxE_LXe);
		  hE_LXe_all 		->SetLineColor(1);
		  hE_LXe_all 		->SetLineWidth(2);
	TH1F *hE_LXe_true 		= new TH1F("hE_LXe_true", "", nbinsE_LXe, minE_LXe, maxE_LXe);
		  hE_LXe_true 		->SetLineColor(2);
		  hE_LXe_true 		->SetLineWidth(2);
	TH1F *hE_LXe_trueCut	= new TH1F("hE_LXe_trueCut", "", nbinsE_LXe, minE_LXe, maxE_LXe);
		  hE_LXe_trueCut	->SetLineColor(4);
		  hE_LXe_trueCut	->SetLineWidth(2);

	TH2F *hE_NaI_LXe_all 	= new TH2F("hE_NaI_LXe_all", "", 1000, minE_LXe, maxE_LXe, 1000, minE_NaI, maxE_NaI);
		  hE_NaI_LXe_all 	->SetMarkerColor(1);
	TH2F *hE_NaI_LXe_true 	= new TH2F("hE_NaI_LXe_true", "", 1000, minE_LXe, maxE_LXe, 1000, minE_NaI, maxE_NaI);
		  hE_NaI_LXe_true 	->SetMarkerColor(2);
	TH2F *hE_NaI_LXe_trueCut= new TH2F("hE_NaI_LXe_trueCut", "", 1000, minE_LXe, maxE_LXe, 1000, minE_NaI, maxE_NaI);
		  hE_NaI_LXe_trueCut->SetMarkerColor(4);


	double tof1		= 3.65;
	double tof2		= 4.10;
	double NaIe1	= 655;
	double NaIe2	= 665;
	
	TLine *l1 = new TLine(tof1, NaIe1, tof1, NaIe2);
		   l1->SetLineColor(4);
		   l1->SetLineWidth(2);
	TLine *l2 = new TLine(tof2, NaIe1, tof2, NaIe2);
		   l2->SetLineColor(4);
		   l2->SetLineWidth(2);
	TLine *l3 = new TLine(tof1, NaIe1, tof2, NaIe1);
		   l3->SetLineColor(4);
		   l3->SetLineWidth(2);
	TLine *l4 = new TLine(tof1, NaIe2, tof2, NaIe2);
		   l4->SetLineColor(4);
		   l4->SetLineWidth(2);
	
	//TCut cutTrueCut = Form("NaI_etot>%g && NaI_etot<%g && TOF>%g && TOF<%g", NaIe1, NaIe2, tof1, tof2);
	//TCut cutTrueCut = "eventid==82544468";
	TCut cutTrueCut = "pre_kinE==0";


	////////////////////////////////////////////////////////////////
	// PLOTS
	////////////////////////////////////////////////////////////////
	
	// energy in NaI VS time of flight
	TCanvas *cETOF_NaI = new TCanvas("cETOF_NaI", "cETOF_NaI", 0, 0, 700, 700);
			 cETOF_NaI->SetFillColor(10);
		t1	->Draw("NaI_etot:TOF>>hETOF_NaI_all", "", "");
		t1	->Draw("NaI_etot:TOF>>hETOF_NaI_true", cutTrue, "sames");
		hETOF_NaI_all	->GetXaxis()->SetTitle("time of flight [ns]");
		hETOF_NaI_all	->GetYaxis()->SetTitle("Energy deposited in NaI [keV]");
		hETOF_NaI_all	->GetXaxis()->CenterTitle(true);
		hETOF_NaI_all	->GetYaxis()->CenterTitle(true);
		hETOF_NaI_all	->GetXaxis()->SetTitleOffset(1.20);
		hETOF_NaI_all	->GetYaxis()->SetTitleOffset(1.25);
		hETOF_NaI_all	->GetXaxis()->SetLabelSize(0.03);
		hETOF_NaI_all	->GetYaxis()->SetLabelSize(0.03);
	TLegend *legETOF_NaI = new TLegend(0.2974138,0.7366071,0.9971264,0.8839286,NULL,"brNDC");
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
	cETOF_NaI->Update();
	//char c1_tofE_png[256];
	//sprintf(c1_tofE_png, "%s/%s_tofE.png", OutFolder, FileName);
	//c1_tofE->SaveAs(c1_tofE_png);
	

/*	// time of flight
	TCanvas *c1_tof = new TCanvas("c1_tof", "c1_tof", 0, 0, 1000, 500);
			 c1_tof->SetFillColor(10);
			 c1_tof->SetLogy();
		t1	->Draw("TOF>>hTOF_all", "", "");
		t1	->Draw("TOF>>hTOF_true", cutTrue, "same");
		t1	->Draw("TOF>>hTOF_trueCut", cutTrue && cutTrueCut, "same");
		hTOF_all	->GetXaxis()->SetTitle("Time of flight [ns]");
		hTOF_all	->GetYaxis()->SetTitle("Entries");
		hTOF_all	->GetXaxis()->CenterTitle(true);
		hTOF_all	->GetYaxis()->CenterTitle(true);
		hTOF_all	->GetXaxis()->SetLabelSize(0.03);
		hTOF_all	->GetYaxis()->SetLabelSize(0.03);
	TLegend *leg1_tof = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
 			 leg1_tof->SetBorderSize(0);
			 leg1_tof->SetTextFont(62);
			 leg1_tof->SetFillStyle(0);
	TLegendEntry *entry1_tof=leg1_tof->AddEntry("hTOF_all",	"all coincident events","");
				  entry1_tof->SetTextColor(1);
  				  entry1_tof=leg1_tof->AddEntry("hTOF_true","true coincidence","");
				  entry1_tof->SetTextColor(2);
  				  entry1_tof=leg1_tof->AddEntry("hTOF_trueCut","true coincidence with a cut","");
				  entry1_tof->SetTextColor(4);
	leg1_tof->Draw();
	c1_tof->Update();
	//char c1_tof_png[256];
	//sprintf(c1_tof_png, "%s/%s_TOF.png", OutFolder, FileName);
	//c1_tof->SaveAs(c1_tof_png);
*/

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
			 //c2_E->SetLogy();
		t1	->Draw("etot>>hE_LXe_all", 		"", "");
		t1	->Draw("etot>>hE_LXe_true", 	cutTrue, "sames");
		//t1	->Draw("etot>>hE_LXe_trueCut", 	cutTrue && cutTrueCut, "sames");
		hE_LXe_all	->GetXaxis()->SetTitle("Energy deposited in LXe [keV]");
		hE_LXe_all	->GetYaxis()->SetTitle("Entries");
		hE_LXe_all	->GetXaxis()->CenterTitle(true);
		hE_LXe_all	->GetYaxis()->CenterTitle(true);
		hE_LXe_all	->GetXaxis()->SetLabelSize(0.03);
		hE_LXe_all	->GetYaxis()->SetLabelSize(0.03);
	TLegend *legE_LXe = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
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
	cE_LXe->Update();
	//char c2_E_png[256];
	//sprintf(c2_E_png, "%s/%s_E_LXe.png", OutFolder, FileName);
	//c2_E->SaveAs(c2_E_png);


	// energy in NaI VS energy in LXe 
	TCanvas *cE_NaI_LXe = new TCanvas("cE_NaI_LXe", "cETOF_NaI", 0, 0, 700, 700);
			 cE_NaI_LXe->SetFillColor(10);
		t1	->Draw("NaI_etot:etot>>hE_NaI_LXe_all", "", "");
		t1	->Draw("NaI_etot:etot>>hE_NaI_LXe_true", cutTrue, "sames");
		//t1	->Draw("NaI_etot:etot>>hE_NaI_LXe_trueCut", cutTrue && cutTrueCut, "sames");
		t1	->Draw("NaI_etot:etot>>hE_NaI_LXe_trueCut", cutTrueCut, "sames");
		hE_NaI_LXe_all	->GetXaxis()->SetTitle("Energy deposited in LXe [keV]");
		hE_NaI_LXe_all	->GetYaxis()->SetTitle("Energy deposited in NaI [keV]");
		hE_NaI_LXe_all	->GetXaxis()->CenterTitle(true);
		hE_NaI_LXe_all	->GetYaxis()->CenterTitle(true);
		hE_NaI_LXe_all	->GetXaxis()->SetTitleOffset(1.25);
		hE_NaI_LXe_all	->GetYaxis()->SetTitleOffset(1.25);
		hE_NaI_LXe_all	->GetXaxis()->SetLabelSize(0.03);
		hE_NaI_LXe_all	->GetYaxis()->SetLabelSize(0.03);
	TLegend *legE_NaI_LXe = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
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
	cE_NaI_LXe->Update();
	//char c1_EE_png[256];
	//sprintf(c1_EE_png, "%s/%s_EE.png", OutFolder, FileName);
	//c1_EE->SaveAs(c1_EE_png);
	

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
