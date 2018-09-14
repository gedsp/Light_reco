#include "../LIB/mylib.h"
#include "../LIB/gains.h"

#include "dpdMaker.h"
#include "../config_reco.h"


void get_voltages(int myRun, int mySubrun, float volts[N_PMT])
{
	TFile *fdb = new TFile(db_file.c_str());
	TTree *t = (TTree*)fdb->Get("ntuple");
	Int_t run, subrun;
	Float_t pmt[N_PMT]={0};
	
	t->SetBranchAddress("run",&run);
	t->SetBranchAddress("subrun",&subrun);
	t->SetBranchAddress("pmt",&pmt); 
	
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
			for (uint k=0; k<N_PMT; k++) volts[k]=pmt[k];
			printf("volts = {");
			for (int k=0; k<N_PMT; k++) printf("%0.2f,",volts[k]);
			printf("\b}\n");
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
	float voltages[N_PMT]={0};
	float gains[N_PMT]={0};
	get_voltages(run,subrun,voltages);
	get_gains(gains,voltages); // get gains of the pmts from the voltages.

  
	TString outfile = Form("%s/dpd-%d-%d.root",dpd_dir.c_str(),run,subrun);
	cout << "DPD outfile = " << outfile.Data() << endl;
	
	//analyze the run and store the resulting ntuple in the file specified by the last argument
	make_dpd(run,subrun,gains,outfile.Data()); 



}
