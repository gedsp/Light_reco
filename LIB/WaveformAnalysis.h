#ifndef WAVEFORM_ANALYSIS
#define WAVEFORM_ANALYSIS

#include <vector>
#include "TMath.h"
class TH1;

namespace WaveformAnalysis
{

	/**Gets the baseline and baseline rms from the desired range of bins.
	   \param hist the histogram
	   \param rms the rms of the bin values (returned to user)
	   \param binMin the first bin to use
	   \param binMax the last bin to use
	   \return the mean value of the bins
	*/
	double baseline(const TH1* hist, double& rms, int binMin=1, int binMax=100);

	/**Gets the minimum/maximum bin the histogram neglecting the last bin,
	   the underflow, and the overflow.
	   \param hist the histogram
	*/
	int minbin(const TH1* hist);
	int maxbin(const TH1* hist);
	
	
	/**Takes the sum of all bins between the start and end times, subtracting pedestal
	  \param hist the histogram
	  \param start the start time
	  \param end the end time
	  \param ped the pedestal value to be subtracted 
	  \return the histogram integral between these points (not taking the bin width into account)
	*/
	double integral_S1(const TH1* hist, double start, double end, double ped);

	/**Takes the absolute valued sum of all bins between the start and end times, subtracting pedestal
	  \param hist the histogram
	  \param start the start time
	  \param end the end time
	  \param ped the pedestal value to be subtracted (if > 0)
	  \return the histogram integral between these points (not taking the bin width into account)
	*/
	double integral_S2(const TH1* hist, double start, double end, double ped=0, bool doWeight=false);

	/**Find a list of peaks between the given bin limits.  A peak here is 
	   defined as the highest point in each pulse passing over the threshold.
	  \param hist the histogram
	  \param threshold the threshold for finding pulses
	  \param minBin the first bin to use
	  \param maxBin the last bin to use
	  \return a list of bins containing the maximum value of each pulse
	*/
	std::vector<int> peaks(const TH1* hist, double threshold, int minBin=-1,int maxBin=-1);

	/**Find a list of valleys between the given bin limits.  A valley here is 
	   defined as the lowest point in each pulse passing below the threshold.
	  \param hist the histogram
	  \param threshold the threshold for finding pulses
	  \param minBin the first bin to use
	  \param maxBin the last bin to use
	  \return a list of bins containing the minimum value of each pulse
	*/
	//std::vector<int> valleys(const TH1* hist, double threshold, int minBin=-1,int maxBin=-1);
	
	std::vector<int> valleys(const TH1* hist, int minBin=-1,int maxBin=-1);
	
	
	// Returns (charge) integral of S1 (over pedestal) in the desired time range (in units of ADC*binwidth)
	inline double calc_S1_charge(const TH1* hist, double ped, double tstart=-0.05, double tend=0.95) 
		{ return integral_S1(hist,tstart,tend,ped); }

	// Returns (charge) integral of S1 (over pedestal) up to bin where wf reaches pedestal
	double calc_S1_charge_m2(const TH1* hist, double ped, int binpeak, int &endbin);

	// Looks for S1 binpeak between minbin and maxbin
	// If minbin and maxbin arguments are not provided, looks for S1 in entire bin range
	int find_S1_binpeak(const TH1* hist, int minbin=0, int maxbin=0);

	// Calculates FWHM (in # of bins)
	double calc_S1_width(const TH1* hist, int binpeak, double ped);
	
	
	double calc_S2_parameters(const TH1* hist, double ped, double t_S1, int &binpeak_S2, int &binavg_S2);
	
	double calc_S2_parameters_m2(const TH1* hist, double t_S1, double &t_start, double &t_end, double &width);
	
	// Returns absolute value (charge) integral of S2 over pedestal in the desired time range (in units of ADC*binwidth)
	inline double calc_S2_charge_m2(const TH1* hist, double ped, double tstart, double tend) 
		{ return integral_S2(hist,tstart,tend,ped); }
	
	TH1F* correct_wvf_histo(const TH1* hist, double ped, double RC_c, double RC_d);

}

#endif
