#include "../config_reco.h"

void IsZombie(int light_run)
{	
	TString infilename = Form("%s/light/dpd-light-%d.root",dpd_dir.c_str(),light_run);
	
	TFile ifile(infilename);
	if (ifile.IsZombie()) cout << "Run " << light_run << " is a zombie" << endl;
	else
	{
		TChain *t = new TChain("dpd");
		t->Add(infilename);
		cout << "Run " << light_run << ": " << t->GetEntries() << " events" << endl;
		ifile.Close();
	}
}
