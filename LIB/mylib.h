
void lets_pause (){
	TTimer * timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
	timer->TurnOn();
	timer->Reset();
	std::cout << "q/Q to quit, other to continuee: ";
	char kkey;
	std::cin.get(kkey);
	if( kkey == 'q' || kkey == 'Q') 			gSystem->Exit(0);
	timer->TurnOff();
	delete timer;
}

int autorangemin( TH1F* myhist)
{
       for (Int_t i=0;i<myhist->GetNbinsX();i++) if (myhist->GetBinContent(i)>1) return i;
	return -1;
}
int autorangemax( TH1F* myhist)
{
       for (Int_t i=myhist->GetNbinsX();i>0;i--) if (myhist->GetBinContent(i)>1) return i;
	return -1;
}

void setlogscale( TCanvas* c1, TMultiGraph *mg, Double_t xlow, Double_t xup)
{
	//c1->SetGrid();
	c1->SetLogy();
	c1->SetLogx();
gPad-> SetLogx();
gPad-> SetLogy();
	TStyle *plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");
	   plain->SetCanvasBorderMode(0);
	   plain->SetPadBorderMode(0);
	   plain->SetPadColor(0);
	   plain->SetCanvasColor(0);
	   plain->SetStatColor(0);
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0000); //no mostrar cuadro fit (0011, etc)

//Lines needed to have the value in Log Scale
  //xlow   = 1000;
  //xup    = 1800;
  Int_t    nbinsx = Int_t(xup - xlow);
//cout << " ho visto xmin" << endl;
  Int_t    step_normal     = 100;
  Int_t    factor_large    = 1;
  //  Int_t    factor_large    = ;
  Double_t  fraction_normal = 0.015;
  Double_t  fraction_large  = 0.035;

  Double_t uxmin = c1->GetUxmin();
  Double_t uxmax = c1->GetUxmax();
  Double_t uymin = c1->GetUymin();
  Double_t uymax = c1->GetUymax();

  // x-axis ticks
//----------------------------------------------------------------------------
  if (c1->GetLogx()) {

    Int_t ixlow      = Int_t(xlow);
    Int_t ixup       = Int_t(xup);
    Int_t step_large = factor_large * step_normal;

    mg->GetXaxis()->SetNdivisions(0);

    Double_t ymin = (c1->GetLogy()) ? pow(10, uymin) : uymin;
    Double_t ymax = (c1->GetLogy()) ? pow(10, uymax) : uymax;

    Double_t bottom_tick_normal = fraction_normal * (uymax - uymin);
    Double_t bottom_tick_large  = fraction_large  * (uymax - uymin);

    Double_t top_tick_normal = fraction_normal * (uymax - uymin);
    Double_t top_tick_large  = fraction_large  * (uymax - uymin);

    if (c1->GetLogy()) {

      bottom_tick_normal = pow(10, uymin + bottom_tick_normal) - pow(10, uymin);
      bottom_tick_large  = pow(10, uymin + bottom_tick_large)  - pow(10, uymin);

      TLine tick;

      //for (Int_t i=xlow; i<=xup; i+=100) {
      for (Int_t i=xlow+100; i<xup; i+=100) {

    Double_t xx = i;

    tick.DrawLine(xx, ymin, xx, ymin + ((i-ixlow+50)%step_large == 0 ? bottom_tick_large : bottom_tick_normal));

      }

      // x-axis labels
//--------------------------------------------------------------------------

      Double_t ylatex = ymin - bottom_tick_normal;

      if ((ixup-ixlow)%step_large >= step_normal) {

    TLatex* latex = new TLatex(xup, ylatex, Form("%.0f", xup));

    latex->SetTextAlign(   23);
    latex->SetTextFont (   42);
    latex->SetTextSize (0.035);
    latex->Draw("same");
      }

      for (Int_t i=ixlow+100; i<ixup; i+=step_large) {

    Double_t xx = i;
    TLatex* latex = new TLatex(xx, ylatex, Form("%.0f", xx));

    latex->SetTextAlign(   23);
    latex->SetTextFont (   42);
    latex->SetTextSize (0.035);
    latex->Draw("same");
      }
    }
  }//end initial if
}

bool lets_check_overflow(string histname, string histfile)
{ 

 	TFile data1(histfile.c_str(), "READ");//,"RECREATE");


 	TH1F* hist = (TH1F*)data1.Get(histname.c_str());

         TAxis *xaxis = hist->GetXaxis();
         float vh[4097];
         float vx[4097];
         float adcfirst=0;
	float adclast=0;
         int nloop=0;
         float vhmax=0;
         float bmx=0;

	float second_max=0;
	float h_second_max=0;
	float bin_mean;


         float ntot=0;
         for (int ic=1; ic<=4096; ic++){
	   //         for (int ic=1; ic<=2048; ic++){
           vh[ic] = hist->GetBinContent(ic);//altura del bin
           vx[ic] = xaxis->GetBinCenter(ic);//valor x del bin
           ntot+=vh[ic];

	   if (vh[ic]>vhmax)
	   {
		vhmax=vh[ic];
		bmx=vx[ic];
	   }

	 }
	float trigger=2;//vhmax*0.005;
	float ini_second_max=4096;

	 for (int ic=1; ic<=4096; ic++){
	  
	   if (vh[ic]>trigger && nloop==0){ 
	     nloop++;
	// First bin with contents
	     adcfirst = vx[ic]; //primer bin con señal
	     float aux=2*bmx-adcfirst;
	     ini_second_max=(int)aux;

	     float aux2=2.05*bmx-adcfirst;
	     h_second_max=vh[(int)aux2];
	     second_max=(int)aux2;

		//cout << "ini second max " << ini_second_max << "bmx y adc" << bmx << " " << adcfirst<< endl;
	   }
      // Maximum of the pedestal
	   if (ic>ini_second_max) {

		if(vh[ic]>h_second_max) 
		{
			h_second_max=vh[ic];
			second_max=vx[ic];
		}

	     //}
	   }
	 }
	data1.Close();

	if (bmx==4095 || second_max==4095)
	{
		cout << "   ------------------------  OVERFLOW -in-"<< histname<< "----------------------" << endl;
		return true;
	}
	else
	{	return false;
	}

}



void print_pdf(int i, int tanda, string voltage_str, float voltage, float maxvolt, string histfile, string backupfile)
{		  
	string canvasname="SPE HI-RES V"+voltage_str;
	 TFile data1(histfile.c_str(), "READ");//,"RECREATE");
	TH1F* hist;

	TCanvas *c1 = new TCanvas(canvasname.c_str(),canvasname.c_str());
	c1->SetTitle(canvasname.c_str());    	
	c1->Divide(3,2);


	for( int ch = 1; ch <= 5; ch++ ){
		string histn="pmt"+std::to_string(ch+tanda*5)+"_hi"+voltage_str;
		cout << "Abro histograma " << histn<< endl;
		 hist = (TH1F*)data1.Get(histn.c_str());
		hist->SetTitle(histn.c_str());
		hist->GetXaxis()->SetRange(autorangemin(hist)*0.9,autorangemax(hist)*1.1);
		c1->cd(ch); hist -> Draw("HIST"); gPad->SetLogy();//pm2->Draw();		
		c1->cd(ch)->Modified(); c1->cd(ch)->Update();
	}
	canvasname=canvasname+".png";
	if(i==0) c1->Print(Form("%s %s %s","results/",backupfile.c_str(),"/HIRES.pdf("),"pdf"); 
	else c1->Print(Form("%s %s %s","results/",backupfile.c_str(),"/HIRES.pdf"),"pdf");
	if (voltage==maxvolt) c1->Print(Form("%s %s %s","results/",backupfile.c_str(),"/HIRES.pdf)"),"pdf");

	cout << "Vamos a los low res "<< endl;

	//lets_pause();
	//c1->Print(canvasname.c_str());
	c1->Close();


	canvasname="SPE LOW-RES V"+voltage_str;
	TCanvas *c = new TCanvas(canvasname.c_str(),canvasname.c_str());
	c->SetTitle(canvasname.c_str());    		
	c->Divide(3,2);


	for( int ch = 1; ch <=5; ch++ ){
		string histn="pmt"+std::to_string(ch+5*tanda);
		histn=histn+"_low"+voltage_str;
		cout << "Abro histograma " << histn<< endl;
		 hist = (TH1F*)data1.Get(histn.c_str());
		cout << "minimobin " << autorangemin(hist)<< endl;
		cout << "maxbin " << autorangemax(hist)<< endl;
		hist->SetTitle(histn.c_str());
		hist->GetXaxis()->SetRange(autorangemin(hist)*0.9,autorangemax(hist)*1.1);
		c->cd(ch); hist -> Draw("HIST");gPad->SetLogy();//pm2->Draw();		
		c->cd(ch)->Modified(); c->cd(ch)->Update();
	}
	canvasname=canvasname+".png";
	if(i==0) c->Print(Form("%s %s %s","results/",backupfile.c_str(),"/LOWRES.pdf("),"pdf"); 
	else c->Print(Form("%s %s %s","results/",backupfile.c_str(),"/LOWRES.pdf"),"pdf");
	if (voltage==maxvolt) c->Print(Form("%s %s %s","results/",backupfile.c_str(),"/LOWRES.pdf)"),"pdf");
	
	//c->Print(canvasname.c_str());
	c->Close();
	data1.Close();

}

void print_fitter_pervoltage(int i, int tanda, string voltage_str, float voltage, float maxvolt, std::vector<string> plots[5],string fecha, int pausa, string backupfile, string mode)
{
	cout << "-------------START-print_fitter_pervoltage---------------"<< endl;
	string canvasname="SPE FIT - Double gaussian method - "+voltage_str+"V - "+fecha;
	TCanvas *c1 = new TCanvas("c1",canvasname.c_str(),1200,800);
	c1->SetTitle(canvasname.c_str());    	
	c1->Divide(3,2);
	TImage* _img[5];

	for( int ch = 0; ch < 5; ch++ ){

		cout << "VEAMOS MI VECTR "<<i<<" "<< plots[ch][i]<< endl;
 
		_img[ch]= TImage::Open(plots[ch][i].c_str());
    		c1->cd(ch+1);
		_img[ch]->Draw("xxx");
   		_img[ch]->SetEditable(kTRUE);
		c1->cd(ch+1)->Modified(); c1->cd(ch+1)->Update();

	}

	canvasname="results/"+backupfile+"/"+mode+voltage_str+".png";
	c1->Print(canvasname.c_str(),"png");
	//if(i==0) c1->Print("results/SPEperV.pdf("); 
	//else c1->Print("results/SPEperV.pdf");
	//if (voltage==maxvolt){ c1->Print("results/SPEperV.pdf)");cout  << "ARCHIVO CERRADO"  << endl;}
	if(pausa==1) lets_pause();
	c1->Delete();
	//lets_pause();
	cout << "--------------print_fitter_pervoltage-END--------------"<< endl;

}



TH1F* quickrebin(TH1F* hist,string histn)
{
	TH1F* hnew;
	cout << autorangemax(hist)<< " " << autorangemin(hist) << endl;

	if (autorangemax(hist)-autorangemin(hist)>3000)
	{
		cout << histn <<" rebinned 1:32" << endl;
		hnew = dynamic_cast<TH1F*>(hist->Rebin(32,histn.c_str()));

	}
	else
	{
		if (autorangemax(hist)-autorangemin(hist)>2000)
		{
			cout << histn <<" rebinned 1:16" << endl;
			hnew = dynamic_cast<TH1F*>(hist->Rebin(16,histn.c_str()));

		}
		else
		{
			if (autorangemax(hist)-autorangemin(hist)>1000)
			{
				cout << histn <<" rebinned 1:8" << endl;
				hnew = dynamic_cast<TH1F*>(hist->Rebin(8,histn.c_str()));

			}
			else
			{
				if (autorangemax(hist)-autorangemin(hist)>500)
				{
					cout << histn <<" rebinned 1:4" << endl;
					hnew = dynamic_cast<TH1F*>(hist->Rebin(4,histn.c_str()));
				}
				else
				{

					if (autorangemax(hist)-autorangemin(hist)>250)
					{
						cout << histn <<" rebinned 1:2" << endl;
						hnew = dynamic_cast<TH1F*>(hist->Rebin(2,histn.c_str()));
					}
					else
					{
						hnew = hist;
		
					}
				}
			}
		}
	}

	return hnew;

}




/*
void fitandprintgains(string round, int tanda,std::vector<Double_t>  xx[5],std::vector<Double_t>  yy[5],std::vector<Double_t>  xxe[5],std::vector<Double_t>  yye[5],string fecha_str, ofstream &ofile2)
{


	TCanvas *c3 = new TCanvas("c3","c3",1200, 800);
	   //c3->Divide(3,2);

	TMultiGraph * mg1 = new TMultiGraph("mg1","");


	std::vector<TGraphErrors*> _graph;
	std::vector<TF1*> _myfit;
	std::vector<bool> _fit_successful;


	std::vector<Double_t> A, errorA, B, errorB, vv, vv9 ,chi2, ndf, covarianzaAB, fA, fB, fB9, errorV, errorV9, errorVtilde, errorVtilde9;

	    Double_t gainn=1.e7;
		Double_t gain9=1.e9; //<-- nominal gain

	

	for (int i=0;i<5;i++)
	{
	   cout << "Creo TGraph " << i << endl;

           _graph.push_back( new TGraphErrors(xx[i].size(),&(xx[i][0]),&(yy[i][0]),&(xxe[i][0]),&(yye[i][0])));
	   _fit_successful.push_back(false);
		
	   _graph[i]->SetName("gr1");
		int auxiliar=i+1+tanda*5;
	   _graph[i]->SetTitle(Form("%s%i","PMT",auxiliar));
	   _graph[i]->SetMarkerStyle(20+i);
	   _graph[i]->SetMarkerColor(2+i);
	   if (i==0) _graph[i]->SetDrawOption("AP");
	   else _graph[i]->SetDrawOption("P");
		 _graph[i]->Draw();
	//lets_pause();
	   
	   _graph[i]->SetFillStyle(0);

		cout << "Creo myfit " << i << endl;

	    _myfit.push_back(new TF1(Form("%s %i","myfit_",i),"pow(10,[0])*pow(x,[1])", 1100, 1800));  //<-- RANGE!
	    _myfit[i]->SetParName(0,"A");
	    _myfit[i]->SetParName(1,"B");
	    _myfit[i]->SetLineStyle(4);
	    _myfit[i]->SetLineWidth(3);
	    _myfit[i]->SetLineColor(i+2);


	   cout << "Ajusto! " << endl;


	    TFitResultPtr r = _graph[i]->Fit(Form("%s %i","myfit_",i),"QESV");
	    TString estado=gMinuit->fCstatu;
	    //float covarianzaAB;
		cout << "resultado ajuste: " << int(r) << " " << estado << endl;
	     //lets_pause();
	    if (int(r)==0  && estado=="SUCCESSFUL"){
			cout << _fit_successful[i] << " "<< i << endl;
			_fit_successful[i]=true;
			cout << _fit_successful[i] << " "<< i << endl << endl;
		 	//TMatrixDSym cov = r->GetCovarianceMatrix();
			//cov.Print();
			//cout << cov.Print();
			//cout << "numero columnas" << cov.GetNrows() << endl;//" " << 	cov.GetSub(1,2)<< endl;
			covarianzaAB.push_back(r->GetCovarianceMatrix()[0][1]);
			cout << " eerror ferist par " << covarianzaAB[i] <<endl;
		}
		else
		{
			covarianzaAB.push_back(0.0);
		}
		//lets_pause();

	    gStyle->SetOptFit(0000);

	    A.push_back(_myfit[i]->GetParameter(0));
	    errorA.push_back(_myfit[i]->GetParError(0));
	    B.push_back(_myfit[i]->GetParameter(1));
	    errorB.push_back(_myfit[i]->GetParError(1));
	    chi2.push_back(_myfit[i]->GetChisquare());
	    ndf.push_back(_myfit[i]->GetNDF());

	    fA.push_back(-1/B[i]);
	    fB.push_back(-(TMath::Log10(gainn)-A[i])/(B[i]*B[i]));
	    fB9.push_back(-(TMath::Log10(gain9)-A[i])/(B[i]*B[i]));

	     errorVtilde.push_back(TMath::Sqrt(fA[i]*fA[i]*errorA[i]*errorA[i]+fB[i]*fB[i]*errorB[i]*errorB[i]+2*fA[i]*fB[i]*covarianzaAB[i]));
	     errorVtilde9.push_back(TMath::Sqrt(fA[i]*fA[i]*errorA[i]*errorA[i]+fB9[i]*fB9[i]*errorB[i]*errorB[i]+2*fA[i]*fB9[i]*covarianzaAB[i]));

		//TMatrixDSym cov = _myfit[i]->GetCovarianceMatrix();
	    //->GetChisquare();
	    //--------------------------------------------------
	    //             HV for nominal gain
	    //--------------------------------------------------
	    //G=10^A*V^B
	    vv.push_back(pow(gainn/(pow(10,A[i])),1/B[i]));
	    vv9.push_back(pow(gain9/(pow(10,A[i])),1/B[i]));
	    errorV.push_back(errorVtilde[i]*vv[i]/TMath::Log10(TMath::E()));
	    errorV9.push_back(errorVtilde9[i]*vv9[i]/TMath::Log10(TMath::E()));

	    cout<<" "<<endl;
	    cout<<"HV(gain="<<gainn<<")"<<"="<<vv[i]<<" "<<"V"<<endl;

	    mg1-> Add(_graph[i]);
	    cout << "Acabo fit " <<endl;
	    //ofile2 << "Pos,Date,Analysis,log_10A,err,14k,Err,V_G=1e7,V_G=1e9" << endl;
	    ofile2 << auxiliar<< "," << fecha_str<< ",jose_autogaussian,"<<Form("%.6lf",A[i])<< ","<<Form("%.6lf",errorA[i])<<","<<Form("%.6lf",B[i])<<","<<Form("%.6lf",errorB[i])<<","<<covarianzaAB[i]<<","<<Form("%.0lf",vv[i])<<","<<Form("%.6lf",errorV[i])<<","<<Form("%.0lf",vv9[i])<<","<<Form("%.6lf",errorV9[i]) <<","<<fA[i]<<","<<fB[i]<<","<<fB9[i]<<","<<errorVtilde[i]<<","<<errorVtilde9[i]<<endl;     

	}

	mg1->Draw("AP");

	float x1=0.02;
	float y1=0.5;
	float x2=x1+0.4;//anchura
	float y2=0.97; //altura
	  
       TLegend *leg = new TLegend(x1,y1,x2,y2);
       //leg->SetHeader(Form("%s %i %s %i %s","PMTs from ",1+tanda*5," to ",5+tanda*5," @ Room Temperature (G=10^{A}V^{B})"));   //<-- change
	leg->SetHeader(Form("%s %i %s %i %s %s %s","PMTs from ",1+tanda*5," to ",5+tanda*5," @ RoomT ",fecha_str.c_str()," (G=10^{A}V^{B})"));   //<-- change

	for (int i=0;i<_graph.size();i++)
	{
	  int auxx3=i+1+tanda*5;
	   
           //leg2->AddEntry(_graph[i],Form("%s %i %s %s  %s%.2lf %s %.2lf","PMT ",PMT_pos[i]," - ",serialnum[i].c_str(), " #chi^{2}/ndf = ",chi2[i],"/",ndf[i]));
		cout << "Fit successful?? " << _fit_successful[i] << endl;
	   if (_fit_successful[i])
	   {
		   leg->AddEntry(_graph[i],Form("%s %i %s %s %s %.0lf\n %s %.0lf\n","PMT",auxx3,"-",Serialnum(round,i+5*tanda+1).c_str()," - #chi^{2}/ndf = ",chi2[i],"/",ndf[i]));       //<-- change
		   leg->AddEntry(_myfit[i],Form("%s %.2lf\n %s %.2lf\n %s %.2lf\n %s %.2lf\n","Fit : A=",A[i],"#pm",errorA[i],"; B=",B[i],"#pm",errorB[i]),"L");
		   leg->AddEntry(_graph[i],Form("%s %.1lf\n %s %.1lf\n %s %.1lf\n %s%.1lf\n %s ","HV (Gain=10^{7}) = ",vv[i],"#pm",errorV[i]," V - HV (Gain=10^{9}) = ",vv9[i],"#pm",errorV9[i],"V"),"");
	   }
	}

    //  leg->AddEntry(myfit,"Fit: G=10^{A}V^{B}: A= -25.5231 #pm 0.0046; B= 10.6243 #pm 0.0014","L");
	    leg->SetTextSize(0.034);
	    leg->SetBorderSize(0.0);
	    leg->Draw();

	setlogscale(c3,mg1,1000,1900);
lets_pause();
	cout << "lets log scale! " << endl;
	c3->Print("results/Superfit_exp.png");



}*/



#include"WaveUtils.h"
// MAL: Added optional parameter to provide pedestal value as input,
// rather than relying on bins 0:30
float get_S1full(TH1F *h, float &S1_width, int &binpeak, float ped=0)
{
	// This functions returns the integral of the S1 over the baseline of a waveform in a range of 2xFWHM (in ADCxbinwidth units)

	//
	WaveUtils wu;

	int peak_pos=0;
	double charge=0;
	short amp=h->GetBinContent(1);
	//cout << "bin width" << h->GetBinWidth(1) << endl;
	//cout << "bin numdth" << h->GetSize()+1 << endl;
	//lets_pause();
	//range: 57200,57400 with a binning of 4
	//range: 57200,57400 with a binning of 8
	if (h->GetSize()==262148) h->GetXaxis()->SetRange(57200,57400);
	if (h->GetSize()==131074) h->GetXaxis()->SetRange(28500,28700);
	if (h->GetSize()==1000) h->GetXaxis()->SetRange(60,150);
	binpeak = h->GetMinimumBin();
	//cout << "maximum bin : "<< binmax;
	h->GetXaxis()->SetRange(0,0);
	//h->Draw("HIST");
	
	if (!ped) ped=wu.FindPedestal(*h, 0, 30);
	float semiamp=0.5*(ped-(float)h->GetBinContent(binpeak));
	//cout << " semiamp" << semiamp<< endl;	
	int S1_end=0;
	int S1_start=0;
	S1_width=0;
	for(int i = binpeak; i <  h->GetSize()+1; i++)
	{
		if(h->GetBinContent(i)>ped-semiamp){S1_width=i; S1_end=i; break;}
		
	}
	for(int i = binpeak; i > 0; i--)
	{
		if(h->GetBinContent(i)>ped-semiamp){S1_width=S1_width-i; S1_start=i; break;}
		
	}



	float q_S1=0;


	//cout << "interpolating start 1. " << S1_start<< " " <<ped-h->GetBinContent(S1_start) << endl;//small
	//cout << "interpolating start 2. " << S1_start+1<< " " <<ped-h->GetBinContent(S1_start+1) << endl; //big
	//cout << "semiamp. " << semiamp << endl; //small


	//TAxis *axis = h->GetXaxis();
	//cout<< "binwidth " <<  axis->GetBinWidth(S1_start)<< endl;

	float x1, y1, x2, y2;
	x1 = -0.5;
	x2 = 0.5;
	y1 =ped-h->GetBinContent(S1_start);
	y2 =ped-h->GetBinContent(S1_start+1);
	float m,n;
	m=(y2-y1)/(x2-x1);
	n=y1-m*x1;

	float start_x=(semiamp-n)/m;
	//cout << "m n  " << m<< " " << n << endl;
	//cout << "interpolated xstart " << start_x << endl;

	//start_x=start_x-0.5;
	//if (start_x>0)q_S1-=(ped-h->GetBinContent(S1_start+1))*start_x;
	//else q_S1-=(ped-h->GetBinContent(S1_start))*start_x;
	//cout << "result " << q_S1<< endl;
	//cout << endl << endl;


	//cout << "interpolating end 1. " << S1_end-1<< " " <<ped-h->GetBinContent(S1_end-1) << endl;
	//cout << "interpolating end 2. " << S1_end<< " " <<ped-h->GetBinContent(S1_end) << endl; 
	//cout << "interpolating semiam " << semiamp << endl; //small


	x1 = -0.5;
	x2 = 0.5;
	y1 =ped-h->GetBinContent(S1_end-1);
	y2 =ped-h->GetBinContent(S1_end);

	m=(y2-y1)/(x2-x1);
	n=y1-m*x1;

	float end_x=(semiamp-n)/m;
	//cout << "m n  " << m<< " " << n << endl;

	//cout << "interpolated xend " << end_x << endl;
	//if (end_x>0)q_S1+=(ped-h->GetBinContent(S1_end))*end_x;
	//else q_S1+=(ped-h->GetBinContent(S1_end-1))*end_x;


	start_x=start_x+(S1_start+0.5);
	end_x=(S1_end-0.5)+end_x;
	//cout << "FWHM window " << start_x << " " << end_x << endl;
	S1_width = (end_x) - (start_x);

	start_x=start_x-0.5*S1_width;
	end_x=end_x+0.5*S1_width;



	//cout << "result q_S1" << q_S1<< " width " << S1_width<< endl;
	//cout << "integration window " << start_x << " " << end_x << endl;

	//cout << "integration window " << (int)start_x << " " << (int)end_x << endl;
	S1_start=(int)start_x;
	S1_end=(int)end_x;

	if (start_x-(int)start_x>0.5)q_S1-=(ped-h->GetBinContent(S1_start+1))*(start_x-(int)start_x-0.5);
	else q_S1-=(ped-h->GetBinContent(S1_start))*(start_x-(int)start_x-0.5);

	if (end_x-(int)end_x>0.5)q_S1+=(ped-h->GetBinContent(S1_end+1))*(end_x-(int)end_x-0.5);
	else q_S1+=(ped-h->GetBinContent(S1_end))*(end_x-(int)end_x-0.5);

	//lets_pause();
	for(int i = S1_start+1; i <= S1_end; i++)
	{
		q_S1+=ped-h->GetBinContent(i);
		
	}

	//cout << S1_start<<" " << S1_end<< " " << q_S1 << endl;	
	//lets_pause();
	return q_S1;
}

float get_S2_gaus(TH1F *h_plot, float &s2_width, float &tau_S2, float &amp_S2, float ped)
{
		TF1 *S2_gaus;	    

		S2_gaus = new TF1(Form("%s","S2_gaus_"),Form("%s","-gaus(0)+[3]"),300,1000);
		S2_gaus -> SetParameter(0,20);// charge
		S2_gaus -> SetParameter(1,400); // center
		S2_gaus -> SetParameter(2,40); // sigma
		S2_gaus -> SetParameter(3,ped); // offset
		S2_gaus -> SetParLimits(0,1,200);
		S2_gaus -> SetParLimits(1,300,950);
		S2_gaus -> SetParLimits(2,7,200);
		S2_gaus -> SetParLimits(3,ped-5,ped+5);

		TFitResultPtr res;
		res = h_plot->Fit(Form("%s","S2_gaus_"),"QRMES","goff",200,1000); // "QRMES"	 
		
		/*
		TCanvas *c1 = new TCanvas();
		h_plot->Draw(); 
		c1->Modified();
		c1->Update();
		lets_pause();   
		*/  

		TF1 *S2_gaus_2;
		S2_gaus_2 = new TF1(Form("%s","S2_gaus_2_"),Form("%s","-gaus(0)+[3]"),S2_gaus->GetParameter(1)-4*S2_gaus->GetParameter(2),S2_gaus->GetParameter(1)+4*S2_gaus->GetParameter(2));
		S2_gaus_2 -> SetParameter(0,S2_gaus->GetParameter(0));
		S2_gaus_2 -> SetParameter(1,S2_gaus->GetParameter(1));
		S2_gaus_2 -> SetParameter(2,S2_gaus->GetParameter(2));
		S2_gaus_2 -> SetParameter(3,S2_gaus->GetParameter(3));

		S2_gaus_2 -> SetParLimits(0,0.4*S2_gaus->GetParameter(0),S2_gaus->GetParameter(0)*2);
		S2_gaus_2 -> SetParLimits(1,S2_gaus->GetParameter(1)-S2_gaus->GetParameter(2),S2_gaus->GetParameter(1)+S2_gaus->GetParameter(2));
		S2_gaus_2 -> SetParLimits(2,0.4*S2_gaus->GetParameter(2),1.7*S2_gaus->GetParameter(2));
		S2_gaus_2 -> SetParLimits(3,ped-5,ped+5);

		TFitResultPtr res_2;

		res_2 = h_plot->Fit(Form("%s","S2_gaus_2_"),"QRMES","goff",S2_gaus->GetParameter(1)-4*S2_gaus->GetParameter(2),S2_gaus->GetParameter(1)+4*S2_gaus->GetParameter(2));

		/*
		h_plot->Draw(); 
		c1->Modified();
		c1->Update();
		lets_pause(); 
		*/   
		 
		 //q_S2_error=q_S2*TMath::Sqrt((S2_gaus_2->GetParError(0)/S2_gaus_2->GetParameter(0))*(S2_gaus_2->GetParError(0)/S2_gaus_2->GetParameter(0))+(S2_gaus_2->GetParError(2)/S2_gaus_2->GetParameter(2))*(S2_gaus_2->GetParError(2)/S2_gaus_2->GetParameter(2)));

		 
		 //q_S2_error*=(2.0/4096.0)*1e-6/50; // S2 in Coulombs
		 tau_S2=S2_gaus_2->GetParameter(1);//us / micros
		 //tau_S2_error=S2_gaus_2->GetParError(1); // in us
		 amp_S2=(ped-S2_gaus_2->Eval(S2_gaus_2->GetParameter(1)))*2/4096;
		 s2_width=S2_gaus_2->GetParameter(2);//us / micros;
		 
		return S2_gaus_2->GetParameter(0)*S2_gaus_2->GetParameter(2)*TMath::Sqrt(2*TMath::Pi());
}


float get_S2(TH1F *h, float &s2_width, float &tau_S2, float &amp_S2, int binpeak, float s1_width, float ped)
{

	double charge=0;

	//cout << "ped " << ped << " h " << h->GetBinContent(binpeak+(int)s1_width+1)<< endl;
	//amp_S2=ped-h->GetBinContent(binpeak);
	//cout << " amp s2 " << amp_S2<< " binpeak " << binpeak <<" width " << s1_width <<" rounded " << (int)s1_width <<  endl;
	amp_S2=ped-h->GetBinContent(binpeak+10*(int)s1_width);
	//cout << " amp s2 " << amp_S2<<  endl;
	for(int i = binpeak+10*(int)s1_width; i <  h->GetSize()-1; i++)
	{
		tau_S2+=TMath::Abs(ped-h->GetBinContent(i))*i;
		charge+=TMath::Abs(ped-h->GetBinContent(i));
		amp_S2=TMath::Max(amp_S2,ped-(float)h->GetBinContent(i));
			//cout << " amp s2 " << amp_S2<<  endl;
	}

	tau_S2=tau_S2/charge; //av. bin of S2


	
	 tau_S2=h->GetBinCenter((int)tau_S2);//us / micros
	amp_S2*=2.0/4096.0;//in Volts
	 s2_width=0.0;//us / micros;
	//cout << "tau " << tau_S2 << " amp " << amp_S2<< " charge "<< charge<< endl;
	//lets_pause();
	return charge;
}

