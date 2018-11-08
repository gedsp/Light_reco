#include "TSystem.h"

#include "../config_reco.h"

void IsZombie(const string outfilename="zombies.txt")
{	
	
	TString dir = dpd_dir+"/light/";
	void* dirp = gSystem->OpenDirectory(dir);
	
	const char* entry;
	const char* filename[1500];
	
	const char* ext = ".root";
	
	Int_t n=0;
	TString str;
	
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
      str = entry;
      if(str.EndsWith(ext)) 
        filename[n++] = entry;
    }
	
	ofstream myfile;
	myfile.open(outfilename);

    for (Int_t i = 0; i < n; i++)
	{
		TFile ifile(gSystem->ConcatFileName(dir,filename[i]));
		if (ifile.IsZombie()) myfile << filename[i] << endl;
	}
	myfile.close();
}
