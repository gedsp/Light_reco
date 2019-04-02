#include "TCanvas.h"
#include "WaveformAnalysis.h"
#include "TH1.h"
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

#define max(x,y) x>y?x:y
#define min(x,y) x<y?x:y
#define hget(h,i) h->GetBinContent(i)
#define hcenter(h,i) h->GetXaxis()->GetBinCenter(i)
#define verase(v,i) v.erase(v.begin()+i)

//double S2_min_time=4.0;
//double S2_max_time=262144*4./1000.;

double WaveformAnalysis::baseline(const TH1* hist, double& rms, int binMin, int binMax)
{

  double base=0;
  rms=0;
  
  int n = binMax-binMin+1;
  for (int i = binMin; i<=binMax; i++) base+=hget(hist,i);
  base/=n;
  
  for (int i = binMin; i<=binMax; i++) rms+=pow(hget(hist,i)-base,2);
  rms/=(n-1);
  rms = sqrt(rms);
  return base;

}

int WaveformAnalysis::minbin(const TH1* hist)
{
  int minbin=0;
  double minval = DBL_MAX; 
  // neglects the underflow, the overflow, and the last bin

  for (int i = 1; i< hist->GetNbinsX(); i++)
  {
    if (hist->GetBinContent(i) < minval)
    {
      minval = hist->GetBinContent(i); 
      minbin = i; 
    }
  }
  
  return minbin;
} 

int WaveformAnalysis::maxbin(const TH1* hist)
{
  int maxbin=0;
  double maxval = -DBL_MAX; 
  // neglects the underflow, the overflow, and the last bin
  for (int i = 1; i< hist->GetNbinsX(); i++)
  {
    if (hist->GetBinContent(i) > maxval)
    {
      maxval = hist->GetBinContent(i); 
      maxbin = i; 
    }
  }
 
  return maxbin;
}

double WaveformAnalysis::integral_S1(const TH1* hist, double start, double end, double ped)
{
  start = max(hist->GetXaxis()->GetXmin(),start);
  end = min(hist->GetXaxis()->GetXmax(),end);
  int startBin = hist->GetXaxis()->FindBin(start);
  int endBin = hist->GetXaxis()->FindBin(end);

  double startCenter = hcenter(hist,startBin);
  double endCenter = hcenter(hist,endBin);
  
  double width = hist->GetXaxis()->GetBinWidth(1);

  double sum = 0;
  
  double startBinFrac = (start-startCenter+0.5*width)/width;
  double endBinFrac = (endCenter - 0.5*width-end)/width;
  
  sum += (ped-hget(hist,startBin)) * startBinFrac;
  sum += (ped-hget(hist,endBin)) * endBinFrac;
  for(int i = startBin+1; i<endBin; i++) sum += ped-hget(hist,i);
  
  return sum;
}


vector<int> WaveformAnalysis::peaks(const TH1* hist, double threshold, int minBin,int maxBin)
{

  minBin = minBin<0?1:minBin;
  maxBin = maxBin<0?hist->GetNbinsX():maxBin;
  vector<int> p;
  double v, lastV=hget(hist,minBin);
  bool isInPeak = false;
  for (int i = minBin; i<= maxBin; i++)
  {

    v = hget(hist,i);
    if (lastV<threshold&&v>=threshold){
      isInPeak=true;
      p.push_back(i);
    }else if (lastV>= threshold&&v<threshold) isInPeak=false;

    if (isInPeak){
      if (v>hget(hist,p[p.size()-1])) p[p.size()-1] = i;
    }

    lastV=v;
  }
  return p;

}

/*
vector<int> WaveformAnalysis::valleys(const TH1* hist, double threshold, int minBin,int maxBin)
{

  minBin = minBin<0?1:minBin;
  maxBin = maxBin<0?hist->GetNbinsX():maxBin;
  vector<int> p;
  double v, lastV=hget(hist,minBin);
  bool isInPeak = false;
  for (int i = minBin; i<= maxBin; i++)
  {

    v = hget(hist,i);
    if (lastV>threshold&&v<=threshold){
      isInPeak=true;
      p.push_back(i);
    }else if (lastV<= threshold&&v>threshold) isInPeak=false;
    
    if (isInPeak){
      if (v<hget(hist,p[p.size()-1])) p[p.size()-1] = i;
    }
    
    lastV=v;
  }
  return p;
}
*/

vector<int> WaveformAnalysis::valleys(const TH1* hist, int minBin,int maxBin)
{
  minBin = minBin<0?1:minBin;
  maxBin = maxBin<0?hist->GetNbinsX():maxBin;
  vector<int> p;
  double v, lastV=hget(hist,minBin);
  bool isInPeak = false;
  int min=0;
  double vmin=0;
  for (int i = minBin; i<= maxBin; i++)
  {
    v = hget(hist,i);
	double bw = hist->GetBinWidth(i);
	double dvdt = (v-lastV)/bw;
    if (dvdt<-60. && !isInPeak) { 
		isInPeak=true; // sharp negative derivative found
		min=i;
		vmin=hget(hist,min);
		int halfwidth=0;
		for (int j=i+1; j<=i+6; j++) if (hget(hist,j)<vmin) min=j;
		// now check if wf recovers up to 60% of peak height in subsequent 2 bins
		for (int j=min+1; j<=min+3; j++) 
		{
			double ph=lastV-vmin;
			//cout << "min at " << hcenter(hist,min) << endl;
			//cout << "    j = " << j-min << ": " << (hget(hist,j)-hget(hist,min))/ph << endl;
			if ((hget(hist,j)-vmin)>0.6*ph)
			{
				isInPeak=false;
				p.push_back(min);
				i=j;
				v=hget(hist,i);
				break;
			}
			isInPeak=false;
		}
	} 
	
    if (isInPeak){
	  if (vmin<hget(hist,p[p.size()-1])) p[p.size()-1] = min;
    }
    
    lastV=v;
  }
  return p;

}


double WaveformAnalysis::calc_S1_charge_m2(const TH1* hist, double ped, int binpeak, int &endbin)
{
	endbin=0;
	
	int rn=6; // 'rebin' factor

	for (int i=binpeak; i<hist->GetNbinsX()+1; i++) 
	{
		//if ((ped-hget(h,i))<0.02*halfamp) { S1_end2=i; break; }
		if (hget(hist,i)>ped) 
		{
			double avg=0;
			for (int j=i-rn; j<=i+rn; j++) avg+=hget(hist,j);
			avg/=(2.*rn+1);
			if (avg>ped) { endbin=i; break; }
		}
	}
	
	double tmin=hcenter(hist,binpeak);
	
	if (endbin) return integral_S1(hist,tmin-0.05,hcenter(hist,endbin),ped);
	else return 0.;
}


int WaveformAnalysis::find_S1_binpeak(const TH1* hist, int minbin=0, int maxbin=0)
{
	int binpeak=0;
	double min=1.E9;
	
	if (!minbin) minbin=1;
	if (!maxbin) maxbin=hist->GetNbinsX()+1;

	for (int i=minbin; i<=maxbin; i++)
	{
		if (hget(hist,i)<min)
		{
			binpeak=i;
			min=hget(hist,i);
		}
	}
	
	return binpeak;
}


double WaveformAnalysis::calc_S1_width(const TH1* hist, int binpeak, double ped)
{
	double width=0;
	
	double halfamp = 0.5*(ped-hget(hist,binpeak));
	
	// find bins corresponding to FWHM
	int S1_start=0;
	int S1_end=0;
	for (int i=binpeak; i<hist->GetNbinsX()+1; i++) 
	{
		if ((ped-hget(hist,i))<halfamp) { S1_end=i; break; }
	}

	for (int i=binpeak; i>0; i--) 
	{
		if ((ped-hget(hist,i))<halfamp) { S1_start=i; break; }
	}
	
	// interpolate to find fraction of x bin width where y = halfamp
	double y1=ped-hget(hist,S1_start);
	double y2=ped-hget(hist,S1_start+1);
	double deltaX = hist->GetBinWidth(S1_start);
	double m = y2-y1;
	
	double deltax_start = (y2-halfamp)/m;
		
	y1=ped-hget(hist,S1_end);
	y2=ped-hget(hist,S1_end-1);
	m=y2-y1;
	
	double deltax_end = (y2-halfamp)/m;
	
	width=(S1_end-1)-(S1_start+1)+deltax_start+deltax_end;
	
	return width;
}

/*
void WaveformAnalysis::calc_wvf_DC(const TH1* hist, double ped, double pedrms, double& DC_c, double& DC_d)
{
	DC_c=0;
	DC_d=0;
	
	int binpeak = hist->GetMaximumBin();
	double tpeak = hcenter(hist,binpeak);
	double binsOverPed=0;
	double binsUnderPed=0;
	double totbins=0;
	for (int i=1; i<binpeak; i++)
	{
		totbins++;
		double bc = hget(hist,i);
		if (bc<(ped-3.*pedrms)) binsUnderPed++;
		else if (bc>(ped+3.*pedrms)) binsOverPed++;
	}
	
	double binsum = binsOverPed+binsUnderPed;
	if (binsum>0)
	{
		DC_c = binsUnderPed/binsum;
		DC_d = binsOverPed/binsum;
	}
	
	cout << "tpeak =  " << tpeak << ", DCunder = " << DC_c << ", DCover = " <<  DC_d << endl;

	return;
}
*/

void WaveformAnalysis::correct_wvf_histo(const TH1* hist, TH1F*& hcorr, double ped, double RC_c, double RC_d)
{
	TH1F* hR = (TH1F*)hist->Clone("hR");
	
	for (int i=1; i<hR->GetNbinsX()+1; i++) hR->SetBinContent(i,hget(hR,i)-ped);
	
	hcorr->Reset();
	hcorr->Clear();
	hcorr->Add(hR);
	
	TH1F* hC   = (TH1F*)hR->Clone("hC");
	
	hC->SetBinContent(1,0.);
	hcorr->SetBinContent(1,0.);
	
	double dt = hcorr->GetBinWidth(1);
	//cout << "dt = " << dt << endl;
	for (int i=2; i<hcorr->GetNbinsX()+1; i++)
	{
		double vri=hget(hR,i);
		double vri1=hget(hR,i-1);
		double vpmti1=hget(hcorr,i-1);

		double vci1 = hget(hC,i-1);
		double vci=(vri<vri1)?vci1+dt/RC_c*(vpmti1-vci1):vci1+dt/RC_d*(vpmti1-vci1); 
		double vpmti=vci+vri;
		hC->SetBinContent(i,vci);
		hcorr->SetBinContent(i,vpmti);
	}
	
	// add pedestal back in for plotting, if desired
	for (int i=0; i<hcorr->GetNbinsX()+2; i++) hcorr->SetBinContent(i,hget(hcorr,i)+ped);
	
	delete hR;
	delete hC;
		
	return;
}

int WaveformAnalysis::find_S2_binpeak(const TH1* hist, double t_start, double t_end)
{
	int binpeak=0;
	double min=1.E9;
	
	int minbin,maxbin=0;
	
	/*
	if (!t_start) t_start=t_S1+S2_min_time;
	if (!t_end) t_end=S2_max_time;
	*/
		
	minbin=hist->GetXaxis()->FindBin(t_start);
	maxbin=hist->GetXaxis()->FindBin(t_end);
	
	for (int i=minbin; i<=maxbin; i++)
	{
		if (hget(hist,i)<min)
		{
			binpeak=i;
			min=hget(hist,i);
		}
	}
	
	return binpeak;
}


double WaveformAnalysis::find_S2_peak_coarse(const TH1* hist, double t_start, double t_end)
{
	// very coarse rebin to find largest S2 feature in wvf
	
	TH1D* htmp = (TH1D*)hist->Clone("htmp");
	TH1F* htmp_rebin = dynamic_cast<TH1F*>(htmp->Rebin(1024,"htmp_rebin"));
	
	float bw = htmp_rebin->GetBinWidth(1);
	
	int bp = find_S2_binpeak(htmp_rebin,t_start+bw,t_end);
	float t_bp = hcenter(htmp_rebin,bp);
	
	delete htmp_rebin;
	delete htmp;
	
	return (double)t_bp;
}

/*
int WaveformAnalysis::find_S2_binpeak(const TH1* hist, double t_S1)
{
	int binpeak=0;
	double min=1.E9;
	int binpeak2=0;
	double min2=1.E9;
	
	int rb=256;
	TH1D* h = (TH1D*)hist->Clone("htmp");
	TH1F* hist_rebin = dynamic_cast<TH1F*>(h->Rebin(rb,"htmp_rebin"));
	
	// first find coarse binpeak using rebinned wvf
	int minbin=hist_rebin->GetXaxis()->FindBin(t_S1+S2_min_time);
	int maxbin=hist_rebin->GetXaxis()->FindBin(S2_max_time);
	
	for (int i=minbin; i<=maxbin; i++)
	{
		if (hget(hist_rebin,i)<min)
		{
			binpeak=i;
			min=hget(hist_rebin,i);
		}
	}
	
	// now find fine binpeak using un-rebinned wvf
	int rebin = hist->GetNbinsX()/hist_rebin->GetNbinsX();
	int mid = hist->GetXaxis()->FindBin(hcenter(hist_rebin,binpeak));
	int minbin2=hist->GetXaxis()->FindBin(t_S1+S2_min_time);
	int maxbin2=hist->GetXaxis()->FindBin(S2_max_time);
	
	for (int i=mid-rebin/2; i<=mid+rebin/2; i++)
	{
		if (i<minbin2 || i>maxbin2) continue;
		
		if (hget(hist,i)<min2)
		{
			binpeak2=i;
			min2=hget(hist,i);
		}
	}
	
	delete hist_rebin;
	delete h;

	return binpeak2;
}
*/

/*
double WaveformAnalysis::integral_S2_corr2(const TH1* hist, double ped, double t_S1, double& Q1, double& Q2)
{
	int endBin = hist->GetMaximumBin();
	double tend_S2 = hcenter(hist,endBin);
	if (tend_S2<t_S1) { tend_S2=S2_max_time; endBin=hist->GetXaxis()->FindBin(S2_max_time);}
	
	double tstart_S2 = t_S1+S2_min_time;
	int startBin = hist->GetXaxis()->FindBin(tstart_S2);
	
	//cout << "tstart = " << tstart_S2 << ", tend = " << tend_S2 << ", " << startBin << ", " << endBin << endl;
	
	double sum=0;
	Q1=0;
	Q2=0;
	for(int i = startBin+1; i<endBin; i++)
	{
		if (ped>hget(hist,i)) sum += fabs(ped-hget(hist,i));
		if (ped>hget(hist,i)) Q1+=ped-hget(hist,i);
		if (ped<hget(hist,i)) Q2+=hget(hist,i)-ped;
	}

	return sum;
}
*/

int WaveformAnalysis::calc_S2_binavg(const TH1* hist, double start, double end, double ped)
{
	int binavg = 0;
    start = max(hist->GetXaxis()->GetXmin(),start);
    end = min(hist->GetXaxis()->GetXmax(),end);
    int startBin = hist->GetXaxis()->FindBin(start);
    int endBin = hist->GetXaxis()->FindBin(end);

    double startCenter = hcenter(hist,startBin);
    double endCenter = hcenter(hist,endBin);
  
    double width = hist->GetXaxis()->GetBinWidth(1);

    double sum,wsum = 0;
  
    double startBinFrac = (start-startCenter+0.5*width)/width;
    double endBinFrac = (endCenter - 0.5*width-end)/width;
	
	wsum+=(ped-hget(hist,startBin)) * startBinFrac * startBin;
	wsum+=(ped-hget(hist,endBin)) * endBinFrac * endBin;
	
	sum+=(ped-hget(hist,startBin)) * startBinFrac;
	sum+=(ped-hget(hist,endBin)) * endBinFrac;
	
  	for (int i = startBin+1; i<endBin; i++) 
	{
		if (ped>hget(hist,i)) 
		{ 
			wsum+=i*(ped-hget(hist,i)); 
			sum+=(ped-hget(hist,i)); 
		}
	}
	
	binavg = (int)(wsum/sum+0.5);
	
	return binavg;
}


double WaveformAnalysis::integral_S2(const TH1* hist, double start, double end, double ped, bool doWeight=false)
{
	// always set to false now that wvfs are corrected
	bool pb = false; // (k==2 || k==3) ? true:false;
	
  start = max(hist->GetXaxis()->GetXmin(),start);
  end = min(hist->GetXaxis()->GetXmax(),end);
  int startBin = hist->GetXaxis()->FindBin(start);
  int endBin = hist->GetXaxis()->FindBin(end);

  double startCenter = hcenter(hist,startBin);
  double endCenter = hcenter(hist,endBin);
  
  double width = hist->GetXaxis()->GetBinWidth(1);

  double sum = 0;
  
  double startBinFrac = (start-startCenter+0.5*width)/width;
  double endBinFrac = (endCenter - 0.5*width-end)/width;
  
  if (doWeight)
  {
	  startBinFrac*=(double)startBin;
	  endBinFrac*=(double)endBin;
  }
  
  if (pb)
  {
	  sum += TMath::Abs(ped-hget(hist,startBin)) * startBinFrac;
	  sum += TMath::Abs(ped-hget(hist,endBin)) * endBinFrac;
	  for(int i = startBin+1; i<endBin; i++) 
	  {
		  double w = doWeight?(double)i:1.;
		  sum += w*TMath::Abs(ped-hget(hist,i));
	  }
  }
  else
  {
	  sum += (ped-hget(hist,startBin)) * startBinFrac;
	  sum += (ped-hget(hist,endBin)) * endBinFrac;
	  for(int i = startBin+1; i<endBin; i++) 
	  {
		  double w = doWeight?(double)i:1.;
		  sum += w*(ped-hget(hist,i));
	  }	  
  }
  
  return sum;
}


/*
double WaveformAnalysis::integral_S2_2(const TH1* hist, double ped, double t_S1)
{
	int binpeak = hist->GetMaximumBin();
	double max = hist->GetMaximum();
	double tpeak = hcenter(hist,binpeak);
	
	double tstart_S2 = t_S1+S2_min_time;
	double tend_S2 = S2_max_time; 
	int startBin = hist->GetXaxis()->FindBin(tstart_S2);
	int endBin = hist->GetXaxis()->FindBin(tend_S2);

	double sum=0;
	
	TF1* fit3 = new TF1("fit3","[0]*[1]^x",0,1100); // "[0]*exp(223-x)/[1]",200,1100);
	fit3->SetParameter(0,ped/pow(ped/max,223./(223.-tpeak)));
	fit3->SetParameter(1,pow(ped/max,1./(223-tpeak)));
	
	for (int i = startBin+1; i<endBin; i++) 
	{
		sum+=fit3->Eval(hcenter(hist,i))-hget(hist,i);
	}
	
	return sum;
}
*/


double WaveformAnalysis::calc_S2_parameters(const TH1* hist, double ped, double tstart, double tend, int &binavg_S2)
{
	binavg_S2=0;
	//double d_binavg_S2=0;
	double charge=0;
	double weighted_charge=0.;
	
	/*
	int startbin = hist->GetXaxis()->FindBin(tstart_S2);
	int endbin = hist->GetXaxis()->FindBin(tend_S2);
	printf("startbin = %d, endbin = %d, tstart_S2 = %f, tend_S2 = %f\n",startbin,endbin,tstart_S2,tend_S2);
	*/

	charge = integral_S2(hist,tstart,tend,ped);
	
	binavg_S2 = calc_S2_binavg(hist,tstart,tend,ped);
	
	return charge;
}



double WaveformAnalysis::calc_S2_parameters_m2(const TH1* hist, const TH1* hist_rebin, double ped, double pedrms, double tS2min,
											   double tS2max, double t_binpeak, double &t_start, double &t_end, double &width)
{
	width=0;
	t_start=0;
	t_end=0;
	
	int S2_start,S2_end,S2_start2,S2_end2=0;	

	int minbin = hist_rebin->GetXaxis()->FindBin(tS2min);
	int maxbin = hist_rebin->GetXaxis()->FindBin(tS2max);
	
	double thresh = ped-0.5*pedrms;
	int startBin = hist_rebin->GetXaxis()->FindBin(t_binpeak); // hcenter(hist,binavg));
	
	for (int i=startBin; i>=minbin; i--)
	{
		if (hget(hist_rebin,i)>thresh)
		{
			if (hget(hist_rebin,i-1)>hget(hist_rebin,i)) continue; // wf still falling
			else { S2_start=i; break; }
		}
		if (i==minbin) S2_start=i;
	}
	
	for (int i=startBin; i<=maxbin; i++)
	{
		if (hget(hist_rebin,i)>thresh)
		{
			if (hget(hist_rebin,i+1)>hget(hist_rebin,i)) continue; // wf still falling
			else { S2_end=i; break; }
		}
		if (i==maxbin) S2_end=i;
	}
	
	if (S2_start) t_start = hcenter(hist_rebin,S2_start);
	if (S2_end) t_end = hcenter(hist_rebin,S2_end);
	
	//cout << "t_start = " << t_start << ", " << hget(hist_rebin,S2_start) << ", t_end = " << t_end << ", " << hget(hist_rebin,S2_end) << endl;
	
	/*
	// add bins until charge integral is constant to within X%
	double charge = integral_S2(hist_rebin,t_start,t_end,ped,pb);
	
	for (int i=S2_start-10; i>=minbin; i-=10)
	{
		if (pb) continue;
		double t_start_new = hcenter(hist_rebin,i);
		double charge_new = integral_S2(hist_rebin,t_start_new,t_end,ped,pb);
		if (fabs(charge_new-charge)/charge>0.01)
		{
			cout << "adding 10 bins to S2_start" << endl;
			S2_start=i;
			charge=charge_new;
			t_start=t_start_new;
		}
	}

	for (int i=S2_end+10; i<=maxbin; i+=10)
	{
		if (pb) continue;
		double t_end_new = hcenter(hist_rebin,i);
		double charge_new = integral_S2(hist_rebin,t_start,t_end_new,ped,pb);
		if (fabs(charge_new-charge)/charge>0.01)
		{
			cout << "adding 10 bins to S2_end" << endl;
			S2_end=i;
			charge=charge_new;
			t_end=t_end_new;
		}
	}
	*/
	
	// now fine tune

	S2_start2 = hist->GetXaxis()->FindBin(hcenter(hist_rebin,S2_start));
	S2_end2 = hist->GetXaxis()->FindBin(hcenter(hist_rebin,S2_end));	
	
	int minbin2 = hist->GetXaxis()->FindBin(tS2min);
	int maxbin2 = hist->GetXaxis()->FindBin(tS2max);
	
	for (int i=S2_start2; i>=minbin2; i--)
	{
		if (hget(hist,i)>thresh)
		{
			if (hget(hist,i-1)>hget(hist,i)) continue; // wf still falling
			else { S2_start2=i; break; }	
		}
		if (i==minbin2) S2_start2=i;
	}
	
	for (int i=S2_end2; i<=maxbin2; i++)
	{
		if (hget(hist,i)>thresh)
		{
			if (hget(hist,i+1)>hget(hist,i)) continue; // wf still falling
			else { S2_end2=i; break; }
		}
		if (i==maxbin2) S2_end2=i;
	}
	
	if (S2_start2) t_start = hcenter(hist,S2_start2);
	if (S2_end2) t_end = hcenter(hist,S2_end2);
	
	//cout << "t_start2 = " << t_start << ", " << hget(hist,S2_start2) << ", t_end2 = " << t_end << ", " << hget(hist,S2_end2) << endl;
	
	if (t_end<t_start) width=0.;
	else width=t_end-t_start;
	
	return integral_S2(hist,t_start,t_end,ped);
}
