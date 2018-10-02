#include "../LIB/pause.h"
#include "../LIB/gains.h"

#include "dpdMaker.h"
#include "../config_reco.h"


void get_voltages(int myRun, int mySubrun, double volts[N_PMT], bool &lemON)
{
	TFile *fdb = new TFile(db_file.c_str());
	TTree *t = (TTree*)fdb->Get("ntuple");
	Int_t run, subrun;
	Float_t pmt[N_PMT]={0};
	Float_t lemu[N_LEM]={0};
	Float_t lemd[N_LEM]={0};
	Float_t sum_lem = 0.;
	
	t->SetBranchAddress("run",&run);
	t->SetBranchAddress("subrun",&subrun);
	t->SetBranchAddress("pmt",&pmt); 
	t->SetBranchAddress("lemu",&lemu);
	t->SetBranchAddress("lemd",&lemd);
	
	Int_t nentries = (Int_t)t->GetEntries();
	Int_t i=0;
	bool found=false;

	while(i<nentries && !found)
	{
		t->GetEntry(i);
		if ((int)run==myRun && (int)subrun==mySubrun)
		{
			found=true;
			cout << "Entry found for run " << myRun << " and subrun " << mySubrun << endl;
			cout << "Sum of all LEM voltages = " << sum_lem << endl;
			if (sum_lem>0.) lemON=true;
			for (uint k=0; k<N_PMT; k++) volts[k]=pmt[k];
			printf("volts = {");
			for (int k=0; k<N_PMT; k++) printf("%0.2f,",volts[k]);
			printf("\b}\n");
			for (uint k=0; k<N_LEM; k++) sum_lem+=lemu[k]+lemd[k];
		}
		i++;
	}
	
	if (!found)
	{
		printf("Run %d and subrun %d not found in db file %s. Exiting\n",myRun,mySubrun,db_file.c_str());
		gSystem->Exit(0);
	}
}


void dpdMaker(int run, int subrun)
{
	double voltages[N_PMT]={0};
	double gains[N_PMT]={0};
	bool lemsON=false;
	get_voltages(run,subrun,voltages,lemsON);
	get_gains(gains,voltages,lemsON); // get gains of the pmts from the voltages.

  
	TString outfile = Form("%s/dpd-%d-%d.root",dpd_dir.c_str(),run,subrun);
	cout << "DPD outfile = " << outfile.Data() << endl;
	
	//analyze the run and store the resulting ntuple in the file specified by the last argument
	make_dpd(run,subrun,gains,outfile.Data()); 

}
