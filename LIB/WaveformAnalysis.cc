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

double S2_max_time=900.;

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

double WaveformAnalysis::integral_S2_corrected(const TH1* hist, const TF1* fc, double start, double end, double ped, double q)
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
  
  	sum += ped-hget(hist,startBin) * startBinFrac;
	for(int i = startBin+1; i<endBin; i++)
	{
		double ped_diff = fc->Eval(q*1.E6*sum);
		double add = ped+ped_diff-hget(hist,i);
		sum += add;
	}
	double ped_diff = fc->Eval(q*1.E6*sum);
  	sum += (ped+ped_diff-hget(hist,endBin)) * endBinFrac;
  	  
    return sum;
}


double WaveformAnalysis::integral_S2(const TH1* hist, double start, double end, double ped=0, bool doWeight=false)
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
  
  if (doWeight)
  {
	  startBinFrac*=(double)startBin;
	  endBinFrac*=(double)endBin;
  }
  
  if (ped>0) 
  {
	  sum += TMath::Abs(ped-hget(hist,startBin)) * startBinFrac;
	  sum += TMath::Abs(ped-hget(hist,endBin)) * endBinFrac;
	  for(int i = startBin+1; i<endBin; i++) 
	  {
		  double w = doWeight?i:1.;
		  sum += w*TMath::Abs(ped-hget(hist,i));
	  }
  }
  else
  {
	  sum += hget(hist,startBin) * startBinFrac;
	  sum += hget(hist,endBin) * endBinFrac;
	  for(int i = startBin+1; i<endBin; i++)
	  {
		  double w = doWeight?i:1.;
		  sum += w*(hget(hist,i));
	  }
  }
  
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
	
	TH1F* h = (TH1F*)hist->Clone("htmp");
	
	if (!minbin && !maxbin) cout << "calc_S1_parameters::Warning: Using full range for S1 calculation" << endl;
	else h->GetXaxis()->SetRange(minbin,maxbin);
	
	binpeak=h->GetMinimumBin();
	
	delete h;
	
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


double WaveformAnalysis::calc_S2_corrected_charge(const TH1* hist, const TF1* fcorr, double ped, double q, double t_S1)
{
	double charge=0;
	
	TH1F* h = (TH1F*)hist->Clone("htmp");
	
	double tstart_S2 = t_S1+4.0;
	double tend_S2 = S2_max_time; 
	int startbin = h->GetXaxis()->FindBin(tstart_S2);
	int endbin = h->GetXaxis()->FindBin(tend_S2);
	
	h->GetXaxis()->SetRange(startbin,endbin-1);
	
	charge = integral_S2_corrected(h,fcorr,tstart_S2,tend_S2,ped,q);
	
	delete h;
	
	return charge;
}


double WaveformAnalysis::calc_S2_parameters(const TH1* hist, double ped, double t_S1, int &binpeak_S2, int &binavg_S2)
{
	binpeak_S2=0;
	binavg_S2=0;
	double d_binavg_S2=0;
	double charge=0;
	
	TH1F* h = (TH1F*)hist->Clone("htmp");
	
	double tstart_S2 = t_S1+4.0;
	double tend_S2 = S2_max_time; 
	int startbin = h->GetXaxis()->FindBin(tstart_S2);
	int endbin = h->GetXaxis()->FindBin(tend_S2);
	
	h->GetXaxis()->SetRange(startbin,endbin-1);
	
	binpeak_S2=h->GetMinimumBin();
	
	//printf("startbin = %d, endbin = %d, tstart_S2 = %f, tend_S2 = %f\n",startbin,endbin,tstart_S2,tend_S2);
	
	charge = integral_S2(h,tstart_S2,tend_S2,ped);
	d_binavg_S2 = integral_S2(h,tstart_S2,tend_S2,ped,true);
    d_binavg_S2/=charge;
	binavg_S2=(int)d_binavg_S2; 
	
	delete h;
	
	return charge;
}

double WaveformAnalysis::calc_S2_parameters_m2(const TH1* hist, double t_S1, double &t_start, double &t_end, double &width)
{
	width=0;
	t_start=0;
	t_end=0;
	
	double pedrms=0;
	double ped = WaveformAnalysis::baseline(hist,pedrms,hist->GetXaxis()->FindBin(218),hist->GetXaxis()->FindBin(228)-1); 
	
	int bin_S1 = hist->GetXaxis()->FindBin(t_S1+4.0);
	int bin_end = hist->GetXaxis()->FindBin(S2_max_time);

	int S2_start=0;
	int S2_end=0;
	
	double thresh = 10.*pedrms;
		
	for (int i=bin_S1; i<bin_end; i++) 
	{
		if ((ped-hget(hist,i))>thresh)
		{
			if (hget(hist,i+1)>hget(hist,i)) continue; // wf falling
			else { S2_start=i; break; }
		}
	}

	for (int i=bin_end; i>bin_S1; i--) 
	{
		if (fabs(ped-hget(hist,i))>thresh) { S2_end=i; break; }
	}
	
	if (S2_start) t_start=hcenter(hist,S2_start);
	if (S2_end) t_end=hcenter(hist,S2_end);
	
	if (t_end<t_start) width=0.;
	else width=t_end-t_start;
	
	return integral_S2(hist,t_start,t_end,ped);
}
