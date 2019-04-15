#include "../LIB/pause.h"
#include "../LIB/gains.h"

#include "dpdMaker.h"
#include "../config_reco.h"

double voltages[N_PMT]={0};
double gains[N_PMT]={0};
double drift_field=0;
bool lemsON=false;
int trig_conf=-2;

bool get_voltages(const string dbfile, int myRun, int mySubrun=-1)
{
	TFile *fdb = new TFile(dbfile.c_str());
	TTree *t = (TTree*)fdb->Get("ntuple");
	Int_t run, subrun;
	Float_t drift_field_corr=0.;
	Float_t pmt[N_PMT]={0};
	Float_t lemu[N_LEM]={0};
	Float_t lemd[N_LEM]={0};
	Float_t sum_lemu,sum_lemd = 0.;
	
	t->SetBranchAddress("run",&run);
	if (mySubrun>=0) t->SetBranchAddress("subrun",&subrun);
	t->SetBranchAddress("trig_conf",&trig_conf);
	t->SetBranchAddress("drift_field_corr",&drift_field_corr);
	t->SetBranchAddress("pmt",&pmt); 
	t->SetBranchAddress("lemu",&lemu);
	t->SetBranchAddress("lemd",&lemd);
	
	Int_t nentries = (Int_t)t->GetEntries();
	Int_t i=0;
	bool found=false;

	while(i<nentries && !found)
	{
		t->GetEntry(i);
		if ((int)run==myRun && (mySubrun==-1 || ((int)subrun==mySubrun)))
		{
			found=true;
			if (mySubrun==-1) cout << "Entry found for light run " << myRun << endl;
			else cout << "Entry found for matched (charge) run " << myRun << " and subrun " << mySubrun << endl;
			
			drift_field = drift_field_corr; 
			
			cout << "Drift field = " << drift_field << " kV/cm" << endl;
			for (uint k=0; k<N_LEM; k++) { sum_lemu+=lemu[k]; sum_lemd+=lemd[k]; }
			cout << "Sum of all LEMu voltages = " << sum_lemu << endl;
			cout << "Sum of all LEMd voltages = " << sum_lemd << endl;
			
			// Always use July 2017 calibration file 
			// Set lemsON=true to use October 2017 "stressed" calibration
			if (sum_lemd>1500. && sum_lemu>25.) lemsON=false; 
			for (uint k=0; k<N_PMT; k++) voltages[k]=pmt[k];
			printf("voltages = {");
			for (int k=0; k<N_PMT; k++) printf("%0.2f,",voltages[k]);
			printf("\b}\n");
		}
		i++;
	}
	
	return found;
}


void dpdMaker(int light_run, int subrun, bool doLight)
{
	cout << "dpdMaker:: making DPD for light run " << light_run << endl;
	
	//TString infilename = Form("%s/output0000%0.4d.root",light_data_dir.c_str(),light_run);
	TString infilename = Form("%s/output0000%0.4d_reprocessed.root",light_data_dir.c_str(),light_run);
	
	cout << "Light data input file: " << infilename.Data() << endl;
	
	TFile ifile(infilename);
	if (ifile.IsZombie())
	{
	       cout << "Error opening input file " << infilename.Data() << endl;
	       gSystem->Exit(-1);
	}
	
 	bool sc = get_voltages(db_light_file,light_run);
	
	if (!sc)
	{
		printf("Light run %d not found in db file %s. Exiting\n",light_run,db_light_file.c_str());
		gSystem->Exit(-1);
	}
	
	get_gains(gains,voltages,lemsON); // get gains of the pmts from the voltages.
	
	// double check that all events belong to the same run since this info is in a separate tree
	TChain * r = new TChain("run_info");
	r->Add(infilename);
	Int_t RunNum;
	r->SetBranchAddress("RunNum",&RunNum);
	Int_t nentries = (Int_t)r->GetEntries();
	int i=0;
	while(i<nentries)
	{
		r->GetEntry(i);
		if (RunNum != light_run) 
		{
			cout << "ERROR: dpdMaker (light): event " << i << " has RunNum = " << RunNum << endl;
			gSystem->Exit(-1);
		}
		i++;
	}
	delete r;
	
	TChain * t = new TChain("midas_data_new");
	t->Add(infilename);
	
	gSystem->Exec(Form("mkdir -p %s/light",dpd_dir.c_str()));
	TString outfile;
	if (subrun>=0) outfile = Form("%s/light/dpd-light-%d-%0.3d.root",dpd_dir.c_str(),light_run,subrun);
	else outfile = outfile = Form("%s/light/dpd-light-%d.root",dpd_dir.c_str(),light_run);
	
	cout << "DPD outfile: " << outfile.Data() << endl;
	
	//analyze the run and store the resulting ntuple in the file specified by the last argument
	make_dpd(t,light_run,trig_conf,drift_field,gains,outfile.Data(),subrun); 
}

void dpdMaker(int light_run)
{	
	dpdMaker(light_run,-1,true);
}

void dpdMaker(int charge_run, int subrun)
{	
	cout << "dpdMaker:: making DPD for matched (charge) run " << charge_run << endl;
		
	TString infilename = Form("%s/%d-%d-RecoFast-Parser.root",matched_data_dir.c_str(),charge_run,subrun);
	
	cout << "Matched data input file: " << infilename.Data() << endl;
	
	TFile ifile(infilename);
	if (ifile.IsZombie())
	{
	       cout << "Error opening input file " << infilename.Data() << endl;
	       gSystem->Exit(-1);
	}
	
	bool sc = get_voltages(db_charge_file,charge_run,subrun);
	
	if (!sc)
	{
		printf("Matched (charge) run %d and subrun %d not found in db file %s. Exiting\n",charge_run,subrun,db_charge_file.c_str());
		gSystem->Exit(0);
	}
	
	get_gains(gains,voltages,lemsON); // get gains of the pmts from the voltages.
	
	TChain * t = new TChain("analysistree/anatree");
	int nf = t->Add(infilename);
	if (nf==0) 
	{
		cout << "ERROR adding file " << infilename << " to TChain" << endl;
		gSystem->Exit(-1);
	}

	gSystem->Exec(Form("mkdir -p %s/matched",dpd_dir.c_str()));	
	TString outfile = Form("%s/matched/dpd-matched-%d-%d.root",dpd_dir.c_str(),charge_run,subrun);
	cout << "DPD outfile: " << outfile.Data() << endl;
	
	//analyze the run and store the resulting ntuple in the file specified by the last argument
	make_dpd(t,charge_run,trig_conf,drift_field,gains,outfile.Data()); 

}
