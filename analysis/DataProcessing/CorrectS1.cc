const char *OutFolder_CorrectS1	= "/Users/jaehzorn/analysis/Xuerich/plots/Co57/CorrectS1";

double dt_step 	= 1.60;
double dt_l 	= 0;
double dt_r 	= dt_l + dt_step;

int n = 10;
double x[n];
double xer[n];
double y[n];
double yer[n];
double y_ly[n];
double yer_ly[n];
double y3p[n];
double y3m[n];

for (int i=0; i<n; i++){

	Char_t  Hname[500];
	sprintf(Hname, "histo_%d", i);
	cout << Hname <<" : "<< dt_l <<" to "<< dt_r <<" us"<< endl;

	// Energy Histogram
	TH1F *h_S1tot = new TH1F("h_S1tot","", nbins_phe, min_phe, max_phe);
		  h_S1tot->SetLineColor(1);
		  h_S1tot->SetLineWidth(2);
	
	// Gaussian Fit
	TF1 *gausFit = new TF1("gausFit","gaus", min_phe, max_phe);
		 gausFit->SetLineColor(2);
		 gausFit->SetLineWidth(2);
	
	// Drift Time Cut
	TCut cut_dt			= Form("dt > %g && dt < %g", dt_l, dt_r);
	
	// PLOT
	TCanvas *c_S1tot = new TCanvas("c_S1tot", "c_S1tot", 0, 0, 700, 500);
			 c_S1tot->SetFillColor(10);
			 c_S1tot->SetGrid();
	T1->Draw("S1tot>>h_S1tot", cutData && cut_dt,"");
		//h_S1tot	->Scale(rescaleF_mHz);
		h_S1tot->GetXaxis()->SetTitle("S1 [phe]");
		h_S1tot->GetXaxis()->CenterTitle(true);
		h_S1tot->GetYaxis()->SetTitle("rate [mHz]");
		h_S1tot->GetYaxis()->CenterTitle(true);
		
		h_S1tot->Fit(gausFit,"R");
		gausFit->Draw("SAME");
		
	c_S1tot->Update();
	
	
	double xmean = dt_l + (dt_r - dt_l)/2;
	x[i] = xmean;
	xer[i] = dt_step/2;
	
	double ymean =  gausFit->GetParameter(1);
	y[i] 		= ymean;
	y_ly[i] 	= ymean / Co57peak;
	double ysigma = gausFit->GetParameter(2);
	yer[i] = ysigma;
	yer_ly[i] = ysigma / Co57peak;
	// +-sigma
	y3p[i] 	= (ymean+ysigma);
	y3m[i] 	= (ymean-ysigma);
	
	dt_l += dt_step;
	dt_r += dt_step;

	Char_t  Sname[500];
	sprintf(Sname, "%s/S1slice_%d.png", OutFolder_CorrectS1, i+1);
	c_S1tot->Print(Sname);

}


	TGraph		 *gr_polfit = new TGraph(n, x, y);
				  gr_polfit->SetMarkerStyle(8);
				  gr_polfit->SetMarkerSize(0.5);
				  gr_polfit->SetMarkerColor(1);
				  gr_polfit->SetLineColor(1);
	TGraph		 *gr_polfit_3p = new TGraph(n, x, y3p);
				  gr_polfit_3p->SetMarkerStyle(8);
				  gr_polfit_3p->SetMarkerSize(0.5);
				  gr_polfit_3p->SetMarkerColor(1);
				  gr_polfit_3p->SetLineColor(1);
	TGraph		 *gr_polfit_3m = new TGraph(n, x, y3m);
				  gr_polfit_3m->SetMarkerStyle(8);
				  gr_polfit_3m->SetMarkerSize(0.3);
				  gr_polfit_3m->SetMarkerColor(1);
				  gr_polfit_3m->SetLineColor(1);

	TGraphErrors *gr_polfit_ly = new TGraphErrors(n, x, y_ly, xer, yer_ly);
				  gr_polfit_ly->SetMarkerStyle(8);
				  gr_polfit_ly->SetMarkerSize(0.3);
				  gr_polfit_ly->SetMarkerColor(1);
				  gr_polfit_ly->SetLineColor(1);


	// plot polfit in phe
	TCanvas *c_polfit = new TCanvas("c_polfit", "c_polfit", 0, 0, 700, 500);
			 c_polfit->SetFillColor(10);
			 c_polfit->SetGrid();
			 c_polfit->SetFillStyle(0);
		gr_polfit		->	Draw("AP");
		gr_polfit_3p	->	Draw("P");
		gr_polfit_3m	->	Draw("P");
		gr_polfit	->	SetTitle(0);
		gr_polfit	->	GetXaxis()	->	SetTitle("drift time [us]");
		gr_polfit	->	GetYaxis()	->	SetTitle("S1 [phe]");
		gr_polfit	-> 	GetXaxis() 	-> 	CenterTitle(true);
		gr_polfit	-> 	GetYaxis() 	-> 	CenterTitle(true);
		gr_polfit	-> 	GetYaxis() 	-> 	SetTitleOffset(1.30);
		gr_polfit		->	SetMinimum(min_phe);
		gr_polfit		->	SetMaximum(max_phe);


	TF1 *f1_polfit5 = new TF1("f1_polfit5","pol3",0.8, 16.0);
		 f1_polfit5->SetLineColor(2);
		 f1_polfit5->SetLineWidth(2);
		gr_polfit	->Fit(f1_polfit5,"R");
		f1_polfit5	->Draw("SAME");

	TF1 *f1_polfit5_3p = new TF1("f1_polfit5_3p","pol3",0.8, 16.0);
		 f1_polfit5_3p->SetLineColor(2);
		 f1_polfit5_3p->SetLineWidth(1);
		gr_polfit_3p	->Fit(f1_polfit5_3p,"R");
		f1_polfit5_3p	->Draw("SAME");

	TF1 *f1_polfit5_3m = new TF1("f1_polfit5_3m","pol3",0.8, 16.0);
		 f1_polfit5_3m->SetLineColor(2);
		 f1_polfit5_3m->SetLineWidth(1);
		gr_polfit_3m	->Fit(f1_polfit5_3m,"R");
		f1_polfit5_3m	->Draw("SAME");

   TPaveStats *pts_polfit = new TPaveStats(0.1879085,0.6747573,0.5245098,0.881068,"brNDC");
			   pts_polfit->SetName("stats");
			   pts_polfit->SetBorderSize(0);
			   pts_polfit->SetFillColor(10);
			   pts_polfit->SetTextAlign(12);
			   pts_polfit->SetTextFont(42);
			   pts_polfit->SetOptStat(0);
			   pts_polfit->SetOptFit(111);
			   pts_polfit->Draw();
		   gr_polfit->GetListOfFunctions()->Add(pts_polfit);

	double S1fit_0 = f1_polfit5->GetParameter(0);
	double S1fit_1 = f1_polfit5->GetParameter(1);
	double S1fit_2 = f1_polfit5->GetParameter(2);
	double S1fit_3 = f1_polfit5->GetParameter(3);
	//double S1fit_4 = f1_polfit5->GetParameter(4);
	//double S1fit_5 = f1_polfit5->GetParameter(5);
		
	c_polfit->Update();

	Char_t  name_polfit[500];
	sprintf(name_polfit, "%s/S1slices_polFit_phe.png", OutFolder_CorrectS1);
	c_polfit->Print(name_polfit);


	// plot polfit in LY units (phe/keVee)
	TCanvas *c_polfit_ly = new TCanvas("c_polfit_ly", "c_polfit_ly", 0, 0, 700, 500);
			 c_polfit_ly->SetFillColor(10);
			 c_polfit_ly->SetGrid();
			 c_polfit_ly->SetFillStyle(0);
		gr_polfit_ly	->	Draw("AP");
		gr_polfit_ly	->	SetTitle(0);
		gr_polfit_ly	->	GetXaxis()	->	SetTitle("drift time [us]");
		gr_polfit_ly	->	GetYaxis()	->	SetTitle("Light yield [phe/keVee]");
		gr_polfit_ly	-> 	GetXaxis() 	-> 	CenterTitle(true);
		gr_polfit_ly	-> 	GetYaxis() 	-> 	CenterTitle(true);
		gr_polfit_ly	-> 	GetYaxis() 	-> 	SetTitleOffset(1.30);
	TF1 *f1_polfit5_ly = new TF1("f1_polfit5_ly","pol6",0.80, 16.0);
		 f1_polfit5_ly->SetLineColor(2);
		 f1_polfit5_ly->SetLineWidth(2);
		gr_polfit_ly->Fit(f1_polfit5_ly,"R");
		f1_polfit5_ly->Draw("SAME");

   TPaveStats *pts_polfit_ly = new TPaveStats(0.5604575,0.1334951,0.8251634,0.4514563,"brNDC");
			   pts_polfit_ly->SetName("stats");
			   pts_polfit_ly->SetBorderSize(0);
			   pts_polfit_ly->SetFillColor(10);
			   pts_polfit_ly->SetTextAlign(12);
			   pts_polfit_ly->SetTextFont(42);
			   pts_polfit_ly->SetOptStat(0);
			   pts_polfit_ly->SetOptFit(111);
			   pts_polfit_ly->Draw();
		   gr_polfit_ly->GetListOfFunctions()->Add(pts_polfit_ly);

	double S1fit_ly_0 = f1_polfit5_ly->GetParameter(0);
	double S1fit_ly_1 = f1_polfit5_ly->GetParameter(1);
	double S1fit_ly_2 = f1_polfit5_ly->GetParameter(2);
	double S1fit_ly_3 = f1_polfit5_ly->GetParameter(3);
	double S1fit_ly_4 = f1_polfit5_ly->GetParameter(4);
	//double S1fit_ly_5 = f1_polfit5_ly->GetParameter(5);
		
	c_polfit_ly->Update();

	Char_t  name_polfit_ly[500];
	sprintf(name_polfit_ly, "%s/S1slices_polFit_ly.png", OutFolder_CorrectS1);
	c_polfit_ly->Print(name_polfit_ly);



	// compute S1 correction
	Char_t  polfitX_char[500];
	//sprintf(polfitX_char, "%f + %f*(x) + %f*(x)*(x) + %f*(x)*(x)*(x) + %f*(x)*(x)*(x)*(x) + %f*(x)*(x)*(x)*(x)*(x)", S1fit_0,S1fit_1,S1fit_2,S1fit_3,S1fit_4,S1fit_5);
	sprintf(polfitX_char, "%f + %f*(x) + %f*(x)*(x) + %f*(x)*(x)*(x)", S1fit_0,S1fit_1,S1fit_2,S1fit_3);
	TF1 *S1polfitX = new TF1("S1polfitX", polfitX_char, min_dt, max_dt);
		 double dtS1_min 	= S1polfitX->GetXmin();
		 double dtS1_max 	= S1polfitX->GetXmax();
		 double dtS1_mean 	= dtS1_min + (dtS1_max - dtS1_min)/2;
		 double S1_mean 	= S1fit_0 + S1fit_1*dtS1_mean + S1fit_2*dtS1_mean*dtS1_mean + S1fit_3*dtS1_mean*dtS1_mean*dtS1_mean;

	Char_t  S1cor_phe_char[500];
	//sprintf(S1cor_phe_char, "S1tot * %f/(%f + %f*dt + %f*dt*dt + %f*dt*dt*dt + %f*dt*dt*dt*dt + %f*dt*dt*dt*dt*dt)", S1_mean, S1fit_0,S1fit_1,S1fit_2,S1fit_3,S1fit_4,S1fit_5);
	sprintf(S1cor_phe_char, "S1tot * %f/(%f + %f*dt + %f*dt*dt + %f*dt*dt*dt)", S1_mean, S1fit_0,S1fit_1,S1fit_2,S1fit_3);

	T1->SetAlias("S1totCor", S1cor_phe_char);

