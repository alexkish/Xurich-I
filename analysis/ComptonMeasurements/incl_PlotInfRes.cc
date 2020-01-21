	// energy in NaI, all events
	TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 700, 500);
			 c1->SetFillColor(10);
			 c1->SetLogy();
		file1_E1all_1	->Draw();
		file2_E1all_1	->Draw("same");
		file3_E1all_1	->Draw("same");

		file1_E1all_1	->Rebin(4);
		file2_E1all_1	->Rebin(4);
		file3_E1all_1	->Rebin(4);

		file1_E1all_1	->Scale(0.25);
		file2_E1all_1	->Scale(0.25);
		file3_E1all_1	->Scale(0.25);

		file1_E1all_1	->SetMinimum(10);

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
	double int_noPipe 	= file1_E1all_1->Integral();
	double int_Pipe50cm = file2_E1all_1->Integral();
	double int_Pipe75cm = file3_E1all_1->Integral();
	double reduction_Pipe50cm = (1-int_Pipe50cm/int_noPipe)*100;
	double reduction_Pipe75cm = (1-int_Pipe75cm/int_noPipe)*100;
	cout <<"int_noPipe   = "<< int_noPipe 	<< endl;
	cout <<"int_Pipe50cm = "<< int_Pipe50cm <<"; reduction = "<< reduction_Pipe50cm << endl;
	cout <<"int_Pipe75cm = "<< int_Pipe75cm <<"; reduction = "<< reduction_Pipe75cm << endl;

	// energy in NaI, all coincident events
	TCanvas *c2a = new TCanvas("c2a", "NaI, all coincident events", 0, 0, 700, 500);
			 c2a->SetFillColor(10);
			 //c2a->SetLogy();
		file1_E1_1	->Draw();
		file2_E1_1	->Draw("same");
		file3_E1_1	->Draw("same");
		
		file1_E1_1	->Rebin(2);
		file2_E1_1	->Rebin(2);
		file3_E1_1	->Rebin(2);

		file1_E1_1	->Scale(0.5);
		file2_E1_1	->Scale(0.5);
		file3_E1_1	->Scale(0.5);
		
		file1_E1_1->GetXaxis()->SetRangeUser(600, 680);
		
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
			 //c2b->SetLogy();
		file1_E1_2	->Draw();
		file2_E1_2	->Draw("same");
		file3_E1_2	->Draw("same");
		
		file1_E1_2	->Rebin(2);
		file2_E1_2	->Rebin(2);
		file3_E1_2	->Rebin(2);

		file1_E1_2	->Scale(0.5);
		file2_E1_2	->Scale(0.5);
		file3_E1_2	->Scale(0.5);

		file1_E1_2->GetXaxis()->SetRangeUser(600, 680);
		
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
			 //c3a->SetLogy();
		file1_E2_1	->Draw();
		file2_E2_1	->Draw("same");
		file3_E2_1	->Draw("same");
		
		file1_E2_1	->Rebin(2);
		file2_E2_1	->Rebin(2);
		file3_E2_1	->Rebin(2);

		file1_E2_1	->Scale(0.5);
		file2_E2_1	->Scale(0.5);
		file3_E2_1	->Scale(0.5);

		file1_E2_1->GetXaxis()->SetRangeUser(0, 100);
		file1_E2_1->SetMaximum(110);
		
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
			 //c3b->SetLogy();
		file1_E2_2	->Draw();
		file2_E2_2	->Draw("same");
		file3_E2_2	->Draw("same");
		
		file1_E2_2	->Rebin(2);
		file2_E2_2	->Rebin(2);
		file3_E2_2	->Rebin(2);

		file1_E2_2	->Scale(0.5);
		file2_E2_2	->Scale(0.5);
		file3_E2_2	->Scale(0.5);

		file1_E2_2->GetXaxis()->SetRangeUser(0, 100);
		file1_E2_2->SetMaximum(110);
		
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

	//const double livetime = 133333; // sec
	//double n_NaI_allEvents_1 = file1_E1_1->GetEntries();
	//double n_NaI_allEvents_2 = file2_E1_1->GetEntries();
	//double n_NaI_allEvents_3 = file3_E1_1->GetEntries();
	//cout <<"NaI, all events:           "<< n_NaI_allEvents_1/livetime <<" Hz"<< endl;
	//cout <<"NaI, all events, pipe50cm: "<< n_NaI_allEvents_2/livetime <<" Hz"<< endl;
	//cout <<"NaI, all events, pipe75cm: "<< n_NaI_allEvents_3/livetime <<" Hz"<< endl;
