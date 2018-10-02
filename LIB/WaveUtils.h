#ifndef _WAVEUTILS_
#define _WAVEUTILS_


#include <iostream>


#include <vector>
#include <math.h>
#include <algorithm>

#include <TMath.h>
#include <TVector.h>
#include <TArray.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

using std::cout;
using std::endl;


class WaveUtils: public TObject {

 public:
  WaveUtils();
  ~WaveUtils() {}
  
  int   FindPeaks(const TH1F&,     int, int, int, float=0, float=0);
  int   FindPeaks(const TProfile&, int, int, int);
  float FindPedestal(const TH1F&, int, int);
  float GetPedestalRMS();

  void Init();
  void Reset();
  
 public:

  struct peak_t {
    int imax;
    int ifirst;
    int ilast;
    float amp;
    float integ;
    bool isolated;
    bool edge;
    int quality;
    TH1F hPeak;
  };

  std::vector<peak_t> peak_list;
  
  int fNumOfPeaks;
  int fNumOfIsoPeaks;
  int fPeakMaxWidth;
  int fPeakMinWidth;
  int fPeakSig;
  

  bool fPedestalOK;
  double fPedestal;
  double fPedestalRMS;
  
  bool fSaveSinglePeakHist;
  
  TH1F fHistPeaks;
  TH1F fHistWave;
  
};


inline WaveUtils::WaveUtils()
{

  Init();
  Reset();

}


inline void WaveUtils::Init()
{


  fHistPeaks.SetNameTitle("hPeaks", "Waveform; time [ns]; amplitude [ADC counts]");
  fHistWave.SetNameTitle("hWave", "Waveform; time [ns]; amplitude [ADC counts]");
  
  Reset();

  fPedestalOK = false;
  fPeakMaxWidth = 8; // 100
  fPeakMinWidth = 1; // 3
  fPeakSig = 20; // 2
  
  fSaveSinglePeakHist = true;
}


inline void WaveUtils::Reset()
{

  peak_list.clear();

  fHistWave.Reset();
  fHistPeaks.Reset();
  
  fNumOfPeaks = 0;
  fNumOfIsoPeaks=0;
  fPedestal   = 0;
  fPedestalRMS = 0;

}


inline float WaveUtils::GetPedestalRMS()
{ return fPedestalRMS;
}
inline float WaveUtils::FindPedestal(const TH1F& hwave, int start, int stop)
{

  fPedestal    = 0;
  fPedestalRMS = 0;
  
  start = TMath::Max(1,start);
  stop  = TMath::Min(stop, hwave.GetNbinsX() );

  if(start>stop) return fPedestal;
  
  for(int i=start; i<=stop; ++i) fPedestal += hwave.GetBinContent(i);

  float norm = 1./(stop-start+1.);
  
  fPedestal *= norm;
  
  for(int i=start; i<=stop; ++i) fPedestalRMS += pow(hwave.GetBinContent(i)-fPedestal,2);
  
  fPedestalRMS = TMath::Sqrt(fPedestalRMS/(stop-start));

  if( fPedestal != 0. && TMath::Abs(fPedestalRMS/fPedestal) > 2e-3) {
    fPedestalOK = false;
  }else{
    fPedestalOK = true;
  }
  
  return fPedestal;


}//FindPedestal



inline int WaveUtils::FindPeaks(const TProfile& hp, int sign, int start, int stop)
{

  //temporary
  
  TH1F h(hp.GetName(),hp.GetTitle()
	 ,hp.GetNbinsX()
	 ,hp.GetXaxis()->GetXmin()
	 ,hp.GetXaxis()->GetXmax()
	 );

  for(int i=1;i<=hp.GetNbinsX();++i)
    h.SetBinContent(i,hp.GetBinContent(i));
  
  return FindPeaks(h,sign,start,stop);
  
  
}//FindPeaks



inline int WaveUtils::FindPeaks(const TH1F& hwave, int sign, int start, int stop, float ped=0, float pedrms=0)
{

  Reset();
  
  if (ped==0 || pedrms==0) FindPedestal(hwave,1,20);
  else { fPedestal=ped; fPedestalRMS=pedrms; fPedestalOK=true;}
  //FindPedestal(hwave,1,start); 
  
  if( !fPedestalOK ) {
	 cout<<"Unable to find Pedestal "<<endl;
    return 0;
  }

  
  int Nb = hwave.GetNbinsX();
    
  
  if ( start<1 )  start = 1;
  
  if ( stop<=0 )  stop = Nb;  
  else stop = TMath::Min(stop,Nb);

  if( start >= stop) return 0;
  
  cout << "start = " << start << ", stop = " << stop << endl;


  fHistPeaks.SetBins(hwave.GetNbinsX(), hwave.GetXaxis()->GetXmin(), hwave.GetXaxis()->GetXmax() );
  fHistWave.SetBins(hwave.GetNbinsX(), hwave.GetXaxis()->GetXmin(), hwave.GetXaxis()->GetXmax() );
  
  TH1F hw(hwave);
  
  TH1F hflag(hwave);
  hflag.Reset();
  
  if( TMath::Abs(sign)!= 1 ) sign = -1; 
     
  if(sign == -1) {
    
    for(int i=1; i<=Nb; ++i) {
      hw.SetBinContent( i, fPedestal+sign*(hw.GetBinContent(i)-fPedestal) );
    }
    
  }

  for(int i=1; i<=Nb; ++i) 
    fHistWave.SetBinContent( i, (hw.GetBinContent(i)-fPedestal) );

  
			     
  int NSIGMA = fPeakSig;
	 
  int minWidth = fPeakMinWidth;
  int maxWidth = fPeakMaxWidth;

  float minAmp = fPedestal + fPedestalRMS*NSIGMA;


  float tol = 0.05;
  
  int maxNPeaks = Nb;

  float pamp=0, camp=0, dmin=0;
  float lamp=0, uamp=0;


  int nTry = 0;
  
  fNumOfPeaks = 0;
  fNumOfIsoPeaks = 0;
  
  while ( fNumOfPeaks < maxNPeaks && nTry<Nb ) {

    nTry++;
    
    int imax = hw.GetMaximumBin();

    if( imax < start || imax > stop ) {
 
		//cout<<" outside range:  "<<" "<<nTry<<" "
		//		  <<fPedestal<<" "<<imax<<" "<<hw.GetBinCenter(imax)<<" "<<start<<" "<<stop
		//			  	  <<" "<<hw.GetBinContent(imax)<<" "<<minAmp<<endl;
      hw.SetBinContent(imax,fPedestal);
      continue;      
    }
    
    float amp = hw.GetBinContent(imax);

    
    if( amp <= minAmp ) {
		cout<<" min amp reached "<<endl;
      break;
    }


    float posmax = hw.GetXaxis()->GetBinCenter(imax);
    
    
    int ulim = TMath::Min(imax+1,stop);
    int llim = TMath::Max(start,imax-1);


    cout<<" try "<<fNumOfPeaks<<"; "<<imax<<" "<<posmax<<" "<<amp<<" "<<minAmp
			<<"; "<<llim<<" "<<ulim<<"; "<<fPedestal<<endl;

    
    uamp = amp;
    while(    ulim<=stop
	      && uamp>=minAmp
	      && hw.GetBinContent(ulim)<uamp
	   ) {
      uamp = hw.GetBinContent(ulim);
      ulim++;
    }
    ulim--;


    lamp = amp;
    while(    llim>=start
	      && lamp>=minAmp
	      && hw.GetBinContent(llim)<lamp
	   ) {

      lamp = hw.GetBinContent(llim);
      llim--;
    }
    llim++;
    

    if(nTry==2) {

		cout<<" --> "<<amp<<"; "<<imax<<" "<<llim<<" "<<ulim<<endl;
		cout<<" "<<hw.GetBinCenter(imax)<<" "<<hw.GetBinCenter(llim)<<" "<<hw.GetBinCenter(ulim)<<endl;
    }
    
    float amin = TMath::Min(lamp,uamp);


    bool isIsolated = TMath::Abs( (lamp+uamp)*0.5 - fPedestal ) <= fPedestalRMS*2;
    
    int pwd = ulim - llim + 1;

    
    if( pwd < minWidth || pwd > maxWidth ) {

		cout<<" narrow/wide peak skip "<<pwd<<" "<<llim<<" "<<ulim<<endl;

      for(int i=llim; i<=ulim; ++i) {
	hw.SetBinContent(i,amin);
	hflag.SetBinContent(i,-1);
      }
      
      continue;
    }
	
    fNumOfPeaks++;

    if(isIsolated) fNumOfIsoPeaks++;

    cout<<" peak: "<<fNumOfPeaks<<" "<<imax<<" "<<posmax<<" "<<amp<<" "<<pwd<<endl;

    
    peak_t pk;
    pk.imax = imax;
    pk.ifirst = llim;
    pk.ilast = ulim;
    pk.amp = amp - fPedestal;
    pk.isolated = isIsolated;
    pk.edge = false;
    pk.quality = 0;

    if(fSaveSinglePeakHist) {
      pk.hPeak.SetNameTitle(Form("hPeak_%i",fNumOfPeaks),";time sample [ns];Amplitude [ADC counts]");
      pk.hPeak.SetBins(hw.GetNbinsX(),hw.GetXaxis()->GetXmin(),hw.GetXaxis()->GetXmax());
    }
    
    for(int i=llim; i<=ulim; ++i)
      if(hflag.GetBinContent(i)==-1) pk.quality++;
    
    for(int i=llim; i<=ulim; ++i) {
      if(i==start || i==stop) pk.edge = true; break;
    }

    float sum = 0;
    
    for(int i=llim; i<=ulim; ++i) {

      //if(isIsolated && pk.quality==0)
      fHistPeaks.SetBinContent(i, hw.GetBinContent(i) - fPedestal) ;
      
      if(fSaveSinglePeakHist) 
	pk.hPeak.SetBinContent(i, hw.GetBinContent(i) - fPedestal) ;
      
      sum += hw.GetBinContent(i) - fPedestal;

      //keep as last
      hw.SetBinContent(i,amin);
      hflag.SetBinContent(i,-1);
    }

    pk.integ = sum;

    peak_list.push_back(pk);
    
  }//while

  cout<<" "<<endl;
  cout<<" Pedestal: mean "<<fPedestal<<" rms "<<fPedestalRMS<<endl;
  cout<<" Min Amp: "<<minAmp<<endl;
  
  cout<<" Num of Peaks: "<<fNumOfPeaks<<" isolated "<<fNumOfIsoPeaks<<endl;

  return fNumOfPeaks;

}//FindPeaks


 

#endif
