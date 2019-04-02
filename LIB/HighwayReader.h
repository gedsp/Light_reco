#include <iostream>
#include "Riostream.h"
#include <sstream>
#include <string>

#include "../config_reco.h"

std::vector<std::vector<float>> HighwayReader(int runrun, int subrunrun);

std::vector<std::vector<float>> HighwayReader(int runrun, int subrunrun)
{
	TString infile = Form("%s/%d/%d-%d-Highway.txt",highway_data_dir.c_str(),runrun,runrun,subrunrun);
	
	cout << "HighwayReader()::infile = " << infile << endl;
	
	//#	Run	Subrun	Event	Track	ChargesSmallBoxView0	ChargesLargeBoxView0	ChargesSmallBoxView1	ChargesLargeBoxView1
	
	int run,subrun,evt,track;
	float csbv0[KMAXNSEG]={0};
	float clbv0[KMAXNSEG]={0};
	float csbv1[KMAXNSEG]={0};
	float clbv1[KMAXNSEG]={0};
	float cbr0[KMAXNSEG]={0};
	float cbr1[KMAXNSEG]={0};
	
	ifstream in;
	std::string tmp;
	
    Int_t nlines = 0;
    in.open(infile);
  
  	char line[50];
    for (int rrrr=0; rrrr<9; rrrr++) {
       in >> line;
       //cout<<rrrr<<"\t"<<line<<endl;
     }
	 
	 //int nevents=0;
	 int ntracks[KMAXNEVENTS]={0};
	 
 	// initialize and allocate result vector
  	std::vector<std::vector<float>> vv_maxCBR;
	vv_maxCBR.resize(KMAXNEVENTS);
	 
    while (in.eof()==false) {
		
		in >> run >> subrun >> evt >> track >>
			csbv0[0] >> csbv0[1] >> csbv0[2] >> csbv0[3] >> csbv0[4] >> csbv0[5] >> csbv0[6] >> csbv0[7] >> csbv0[8] >> csbv0[9] >> 
				clbv0[0] >> clbv0[1] >> clbv0[2] >> clbv0[3] >> clbv0[4] >> clbv0[5] >> clbv0[6] >> clbv0[7] >> clbv0[8] >> clbv0[9] >>
					csbv1[0] >> csbv1[1] >> csbv1[2] >> csbv1[3] >> csbv1[4] >> csbv1[5] >> csbv1[6] >> csbv1[7] >> csbv1[8] >> csbv1[9] >> 
						clbv1[0] >> clbv1[1] >> clbv1[2] >> clbv1[3] >> clbv1[4] >> clbv1[5] >> clbv1[6] >> clbv1[7] >> clbv1[8] >> clbv1[9];
		
  	  //if (!in.good()) break;
	  
	  // calculations
		float max_CBR=ERRVAL;
			
		for (int i=0; i<KMAXNSEG; i++)
		{
			if (clbv0[i]==ERRVAL || csbv0[i]==ERRVAL || csbv0[i]==0) continue;
			float cbr0 = (clbv0[i]-csbv0[i])/csbv0[i];
			if (cbr0>max_CBR) max_CBR=cbr0; 
		}
		for (int i=0; i<KMAXNSEG; i++)
		{
			if (clbv1[i]==ERRVAL || csbv1[i]==ERRVAL || csbv1[i]==0) continue;
			float cbr1 = (clbv1[i]-csbv1[i])/csbv1[i];
			if (cbr1>max_CBR) max_CBR=cbr1;
		}
		
		// now fill result vector
		ntracks[evt]=track+1;
		vv_maxCBR.at(evt).resize(ntracks[evt]);
		vv_maxCBR.at(evt).at(ntracks[evt]-1)=max_CBR;
		
		//if (nlines<10 || nlines>2720) printf("line %d: %d %d %d %d %f %f %f %f\n",nlines,run,subrun,evt,track,csbv0[0],clbv0[0],csbv1[0],clbv1[0]);
		
		//if (nlines<50) printf("%d\t%d\t%d\t%d\t%f\n",run,subrun,evt,track,max_CBR);
		nlines++;
    }
	
	cout << "HighWayReader()::total nlines = " << nlines << endl;
	in.close();
	
	/*
	for (int ev=0; ev<nevents; ev++)
	{
		vv_maxCBR.at(ev).resize(ntracks[ev]);
		for (int t=0; t<ntracks[ev]; t++)
		{
			if (arr_maxCBR[ev][t]==0 || arr_maxCBR[ev][t]==ERRVAL) vv_maxCBR.at(ev).at(t)=ERRVAL;
			else vv_maxCBR.at(ev).at(t)=arr_maxCBR[ev][t];
		}
	}
	*/

	cout << "HighwayReader()::done" << endl;

	return vv_maxCBR;
}
