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


	TF1 *f1_NaIres = new TF1("f1_NaIres","[0]/sqrt(x)+[1]", 0, 3000);
		 f1_NaIres->SetLineColor(2);
		 f1_NaIres->SetLineWidth(2);
	 
	f1_NaIres->SetParameter(0, 46);
	f1_NaIres->SetParameter(1, 1.3);


	// XENON100 energy resolution (combined scale)
 	TF1 *f1_Xe100res = new TF1("f1_Xe100res","[0]/sqrt(x)+[1]", 0, 3000);
		 f1_Xe100res->SetLineStyle(9);
		 f1_Xe100res->SetLineColor(4);
		 f1_Xe100res->SetLineWidth(2);
	 
	f1_Xe100res->SetParameter(0, 48.5);
	f1_Xe100res->SetParameter(1, 0.9);

	TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 1000, 500);
			 c1->SetFillColor(10);
			 //c1->SetLogy();
		f1_NaIres	->Draw();
		f1_Xe100res	->Draw("same");
		f1_NaIres	->GetXaxis()->SetTitle("Energy [keV]");
		f1_NaIres	->GetYaxis()->SetTitle("Resolution [%]");
		f1_NaIres	->GetXaxis()->CenterTitle(true);
		f1_NaIres	->GetYaxis()->CenterTitle(true);
		f1_NaIres	->GetXaxis()->SetLabelSize(0.03);
		f1_NaIres	->GetYaxis()->SetLabelSize(0.03);
	TLegend *leg1 = new TLegend(0.6124498,0.6716102,0.8895582,0.8813559,NULL,"brNDC");
 			 leg1->SetBorderSize(0);
			 leg1->SetTextFont(62);
			 leg1->SetFillStyle(0);
	TLegendEntry *entry1=leg1->AddEntry("f1_NaIres","NaI resolution","");
				  entry1->SetTextColor(2);
  				  entry1=leg1->AddEntry("f1_Xe100res","Xe100 resolution (combined scale)","");
				  entry1->SetTextColor(4);
	leg1->Draw();
	
	f1_NaIres->SetMaximum(10);
	f1_NaIres->SetMinimum(0);

	f1_NaIres	->SetTitle(0);
	f1_Xe100res	->SetTitle(0);
	
	c1->Update();




/*	TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 1000, 500);
			 c2->SetFillColor(10);
			 //c2->SetLogy();
		f1_Xe100res	->Draw();
		f1_Xe100res	->GetXaxis()->SetTitle("Energy [keV]");
		f1_Xe100res	->GetYaxis()->SetTitle("Resolution [#sigma/#mu]");
		f1_Xe100res	->GetXaxis()->CenterTitle(true);
		f1_Xe100res	->GetYaxis()->CenterTitle(true);
		f1_Xe100res	->GetXaxis()->SetLabelSize(0.03);
		f1_Xe100res	->GetYaxis()->SetLabelSize(0.03);
	c2->Update();
*/

}