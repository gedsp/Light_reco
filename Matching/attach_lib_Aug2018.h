#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <iostream>
#include <cstdlib>

const int MAXLRUNS = 100;

void lets_pause();

void attach_charge_run(int runnum_C, int Lruns, int run_number_L[100],int Lstarttime[100], int Lendtime[100]); // read all the subruns and call one by one the attach function
void attach_light_charge(int, int , string ); // attach one subrun of charge all contained in a light run 
int assing_light_run(string charge_run,int runnum_C,int LRUNS, int light_runs[MAXLRUNS],int Lstarttime[MAXLRUNS],int Lendtime[MAXLRUNS]);
void timing_charge(string charge_run, int run_number, int &starttime, int &endtime,string option); // gives the start and endtime of a charge file
void timing_light(int run_number, int &starttime, int &endtime,string option); // gives the start and end time of a light file



void attach_charge_run(int runnum_C, int Lruns, int run_number_L[100],int Lstarttime[100], int Lendtime[100])
{
        cout << "Analyzing charge run: " << runnum_C << " " << Form("ls %s%i%s",indirC.c_str(),runnum_C,"/>files_.txt")<<endl;
        std::system(Form("ls %s%i/>files_%i.txt",indirC.c_str(),runnum_C,runnum_C));

	ifstream ifile(Form("files_%i.txt",runnum_C));
	string charge_run;
	while (getline(ifile,charge_run))
	{
		 cout << endl<< endl;
		 cout <<"----Processing run: "<< charge_run << endl;

		int assignedrunL = assing_light_run(charge_run,runnum_C,Lruns,run_number_L,Lstarttime, Lendtime);
		cout <<"\tassingedrunL: "<< assignedrunL << endl;

		if (assignedrunL!=-1) {attach_light_charge(assignedrunL,runnum_C,charge_run);std::system(Form("mv %s%s %s%s",tempfolder.c_str(),charge_run.c_str(),outfolder.c_str(),charge_run.c_str()));}
		else cout << "\tsubrun not processed..." << endl;

		
	}

}

void attach_light_charge(int run_number_L, int runnum_C, string ifile){

	bool debug=true;
	debug=false;

	std::system(Form("%s%s%i%s%s%s%s","cp -v ",indirC.c_str(),runnum_C,"/",ifile.c_str()," ",tempfolder.c_str()));



	TFile *fin=TFile::Open(TString::Format("%s%s%08d%s",indirL.c_str(),"output",run_number_L,".root"),"READ"); //PMTrun file
	
	if (!fin || fin->IsZombie()) {
	   cout <<"The file could not be opened!";
	   delete fin;
	} 
	TTree *t = (TTree*)fin->Get("midas_data");
	

	//actual branches at the ntuple:

	short _adc_value[5][300000]; //short or int!! check!
	int _PCTimeTag[3];
	int _nsamples;
	int _time_stamp;
    	int _nevent;
    	int _nchannels;
  	int _crt_daq_match;
  	int _crt_reco;
   	int _time_sample;

	TBranch * _b_time_stamp;  	t->SetBranchAddress("TimeStamp"   , &_time_stamp   , &_b_time_stamp  );
	//TBranch * _b_time_stamp;  t->SetBranchAddress("TimeStamp"   , &_time_stamp   , &_b_time_stamp  );
	TBranch * _b_event;	  	t->SetBranchAddress("event"       , &_nevent       , &_b_event       );
	TBranch * _b_nchannels;   	t->SetBranchAddress("nchannels"   , &_nchannels    , &_b_nchannels   );
	TBranch * _b_time_sample; 	t->SetBranchAddress("TimeSample"  , &_time_sample  , &_b_time_sample );
	TBranch * _b_nsamples;	  	t->SetBranchAddress("nsamples"    , &_nsamples     , &_b_nsamples    );	
	TBranch * _b_adc_value_0; 	t->SetBranchAddress("adc_value_0" , _adc_value[0]  , &_b_adc_value_0 );
	TBranch * _b_adc_value_1; 	t->SetBranchAddress("adc_value_1" , _adc_value[1]  , &_b_adc_value_1 );
	TBranch * _b_adc_value_2; 	t->SetBranchAddress("adc_value_2" , _adc_value[2]  , &_b_adc_value_2 );
	TBranch * _b_adc_value_3; 	t->SetBranchAddress("adc_value_3" , _adc_value[3]  , &_b_adc_value_3 );
	TBranch * _b_adc_value_4; 	t->SetBranchAddress("adc_value_4" , _adc_value[4]  , &_b_adc_value_4 );
	TBranch * _b_PCTimeTag;   	t->SetBranchAddress("PCTimeTag"   , _PCTimeTag     , &_b_PCTimeTag   );
	TBranch * _b_crt_daq_match;	t->SetBranchAddress("crt_daq_match"	, &_crt_daq_match	, &_b_crt_daq_match   );
	TBranch * _b_crt_reco;	t->SetBranchAddress("crt_reco"	, &_crt_reco	, &_b_crt_reco   );
	TBranch * _b_crt_adc      ; int    _crt_adc[4][32]      ; t->SetBranchAddress("crt_adc"       , _crt_adc          , &_b_crt_adc       );

	TBranch * _b_crt_track_pos0      ; float    _crt_track_pos0[3]     ; t->SetBranchAddress("crt_track_pos0"       , _crt_track_pos0          , &_b_crt_track_pos0       );
	TBranch * _b_crt_track_pos1      ; float    _crt_track_pos1[3]     ; t->SetBranchAddress("crt_track_pos1"       , _crt_track_pos1          , &_b_crt_track_pos1       );
	TBranch * _b_crt_ToF     ; float    _crt_ToF     ; t->SetBranchAddress("crt_ToF"       , &_crt_ToF         , &_b_crt_ToF       );

TBranch * _b_crt_track_hits_id   ; int    _crt_track_hits_id[4]   ; t->SetBranchAddress("crt_track_hits_id"   , _crt_track_hits_id       , &_b_crt_track_hits_id       );
TBranch * _b_crt_track_param     ; float  _crt_track_param[2]     ; t->SetBranchAddress("crt_track_param"     , _crt_track_param         , &_b_crt_track_param       );
//TBranch * _b_crt_track_pos0      ; float  _crt_track_pos0[3]      ; t->SetBranchAddress("crt_track_pos0"      , _crt_track_pos0          , &_b_crt_track_pos0       );
//TBranch * _b_crt_track_pos1      ; float  _crt_track_pos1[3]      ; t->SetBranchAddress("crt_track_pos1"      , _crt_track_pos1          , &_b_crt_track_pos1       );
TBranch * _b_crt_isFV            ; int    _crt_isFV               ; t->SetBranchAddress("crt_isFV"            , &_crt_isFV               , &_b_crt_isFV       );
TBranch * _b_crt_track_lenFV     ; float    _crt_track_lenFV      ; t->SetBranchAddress("crt_track_lenFV"     , &_crt_track_lenFV        , &_b_crt_track_lenFV       );
TBranch * _b_crt_point_in        ; float    _crt_point_in[3]      ; t->SetBranchAddress("crt_point_in"        , _crt_point_in            , &_b_crt_point_in       );
TBranch * _b_crt_point_out       ; float    _crt_point_out[3]     ; t->SetBranchAddress("crt_point_out"       , _crt_point_out           , &_b_crt_point_out       );
TBranch * _b_crt_pmt_dist        ; float    _crt_pmt_dist[5]      ; t->SetBranchAddress("crt_pmt_dist"        , _crt_pmt_dist            , &_b_crt_pmt_dist       );



	TFile *fin2=TFile::Open(Form("%s%s",tempfolder.c_str(),ifile.c_str()),"UPDATE"); //PMTrun file
	TTree *t2 = (TTree*)fin2->Get("analysistree/anatree");

	int _event;
    	int _eventtime_seconds;
    	int _eventtime_nanoseconds;
    	//int _no_hits;

	TBranch * _b2_event;	  		t2->SetBranchAddress("EventNumberInRun"       		, &_event       		, &_b2_event       );
	TBranch * _b2_eventtime_seconds;	t2->SetBranchAddress("EventTimeSeconds"		, &_eventtime_seconds		, &_b2_eventtime_seconds       );
	TBranch * _b2_eventtime_nanoseconds;	t2->SetBranchAddress("EventTimeNanoseconds"	, &_eventtime_nanoseconds       , &_b2_eventtime_nanoseconds       );
	//TBranch * _b2_no_hits;		t2->SetBranchAddress("NumberOfHits"		, &_no_hits			, &_b2_no_hits       );

	//newbranches
	const Int_t MaxSamples = 300000;
	short _adc_value_2[5][MaxSamples];
	int _PCTimeTag_2[3];
	int _nsamples_2;
	int _time_stamp_2;
    	int _nevent_2;
    	int _nchannels_2;
  	int _crt_daq_match_2;
  	int _crt_reco_2;
    	int _time_sample_2;
	int _run_2;
        float _deltatime;

	TBranch * _b2_deltatime =       t2->Branch("deltatime"		, &_deltatime	   	, "deltatime/F"       			);
	TBranch * _b2_run = 		t2->Branch("light_run"		, &_run_2	   	, "light_run/I"       			);
	TBranch * _b2_light_event = 	t2->Branch("light_event"	, &_nevent_2		, "light_event/I"     			);
	TBranch * _b2_nchannels = 	t2->Branch("light_nchannels"	, &_nchannels_2		, "light_nchannels/I"   		);
	TBranch * _b2_time_sample = 	t2->Branch("light_TimeSample" 	, &_time_sample_2  	, "light_time_sample/I" 		);
	TBranch * _b2_nsamples = 	t2->Branch("light_nsamples"   	, &_nsamples_2     	, "light_nsamples/I"    		);	
	TBranch * _b2_adc_value_0 = 	t2->Branch("light_adc_value_0"	, _adc_value_2[0]  	, "light_adc_value_0[light_nsamples]/S" 	);
	TBranch * _b2_adc_value_1 = 	t2->Branch("light_adc_value_1"	, _adc_value_2[1]  	, "light_adc_value_1[light_nsamples]/S" 	);
	TBranch * _b2_adc_value_2 = 	t2->Branch("light_adc_value_2"	, _adc_value_2[2]  	, "light_adc_value_2[light_nsamples]/S" 	);
	TBranch * _b2_adc_value_3 =  	t2->Branch("light_adc_value_3"	, _adc_value_2[3]  	, "light_adc_value_3[light_nsamples]/S" 	);
	TBranch * _b2_adc_value_4 =  	t2->Branch("light_adc_value_4"	, _adc_value_2[4]  	, "light_adc_value_4[light_nsamples]/S" 	);
	TBranch * _b2_PCTimeTag = 	t2->Branch("light_PCTimeTag"  	, _PCTimeTag_2     	, "light_PCTimeTag[3]/I"   		);
	TBranch * _b2_crt_daq_match = 	t2->Branch("light_crt_daq_match", &_crt_daq_match_2	, "light_crt_daq_match/I"		);
	TBranch * _b2_crt_reco = 	t2->Branch("light_crt_reco"	, &_crt_reco_2		, "light_crt_reco/I"			);
	int    _crt_adc2[4][32];	TBranch * _b2_crt_adc = 	t2->Branch("crt_adc" 	, _crt_adc2, 		"light_crt_adc[4][32]/I"     );
	float    _crt_track_pos0_2[3];	TBranch * _b2_crt_track_pos0 = 	t2->Branch("crt_track_pos0" , _crt_track_pos0_2, "crt_track_pos0[3]/F"     );
	float    _crt_track_pos1_2[3];	TBranch * _b2_crt_track_pos1 = 	t2->Branch("crt_track_pos1" , _crt_track_pos1_2, "crt_track_pos1[3]/F"     );
	float    _crt_ToF_2;	 	TBranch * _b2_crt_ToF = 	t2->Branch("crt_ToF" 	, &_crt_ToF_2, 		"crt_ToF/F"     );

//	t2->Branch("crt_ToF"        , _crt_ToF      );
//	t2->Branch("crt_adc"        , _crt_adc      ); 
//	t2->Branch("crt_track_pos0"      , _crt_track_pos0    );
//	t2->Branch("crt_track_pos1"      , _crt_track_pos1    );
        TBranch * _b2_track_hits_id = t2->Branch("crt_track_hits_id"   , _crt_track_hits_id );
	TBranch * _b2_crt_track_param = t2->Branch("crt_track_param"     , _crt_track_param   );
	TBranch * _b2_crt_isFV = t2->Branch("crt_isFV"            , &_crt_isFV         );
	TBranch * _b2_crt_track_lenFV = t2->Branch("crt_track_lenFV"     , &_crt_track_lenFV  );
	TBranch * _b2_crt_point_in = t2->Branch("crt_point_in"        , _crt_point_in      );
	TBranch * _b2_crt_point_out = t2->Branch("crt_point_out"       , _crt_point_out     );
	TBranch * _b2_crt_pmt_dist = t2->Branch("crt_pmt_dist"        , _crt_pmt_dist      );


	int numentries=t->GetEntriesFast();

	t2->GetEntry(t2->GetEntries()-1);
	t->GetEntry(t->GetEntries()-1);
	cout << "Last time for both files: t_charge: "  <<_eventtime_seconds <<" t_light: "<<_PCTimeTag[1]<< endl;
	if(_PCTimeTag[1] < _eventtime_seconds-100){cout << "This light file is not gonna work, we need an newer one. Please cancel." <<_eventtime_seconds <<" "<<_PCTimeTag[1]<<   endl; lets_pause();};

	int var=0;
	t2->GetEntry(0); //light
	t->GetEntry(0);  //charge
	cout << "Starting time for both files: t_charge: "  <<_eventtime_seconds <<" t_light: "<<_PCTimeTag[1]<< endl;
	if(_PCTimeTag[1] > _eventtime_seconds+100){cout << "This var file is not gonna work, we need an older one. Please cancel. "  <<_eventtime_seconds <<" "<<_PCTimeTag[1]<<   endl; lets_pause();};

	double time_charge=0.0;
	double time_light=0.0;
	double old_time_light=0.0;
	for(int ev=0; ev < t2->GetEntries() ; ++ev) {

		t2->GetEntry(ev);
		time_charge=(double)_eventtime_seconds-35.0+(double)_eventtime_nanoseconds*1e-9;

		if(debug) cout << "we are working now with charge event: \t" << ev  <<" \ttime_stamp: \t" << Form("%9.3f",time_charge) << endl;

		t->GetEntry(var);
		time_light=(double)_PCTimeTag[0]+(double)_PCTimeTag[2]*1e-3;
		if(debug) cout << "\tLet's find the match in light: \t" << _nevent  << " \ttime_stamp: \t" << Form("%9.3f",time_light) << " " << time_light-time_charge << endl;
		if (time_charge==time_light) {}
		else
		{
			do{
			
				if(time_light-time_charge<-3000) {var+=1000;}
				if(time_light-time_charge<-500) {var+=100;}
				var++;
				t->GetEntry(var);
				time_light=(double)_PCTimeTag[0]+(double)_PCTimeTag[2]*1e-3;
				if(debug) cout << "\tLet's find the match in light: \t" << _nevent  << " \ttime_stamp: \t" << Form("%9.3f",time_light) << " " << time_light-time_charge << endl;
				//if(debug) lets_pause();

				old_time_light=time_light;

				//if (ev>85) { cout << "WARNING! EVENTS TOO FAR AWAY" << endl; lets_pause();cout << "\tLet's find the match in light: \t" << var  << " \ttime_stamp: \t" << Form("%9.3f",time_light) << " " << time_light-time_charge << endl;}
					
			}while(time_light-time_charge<-0.002 && var<t->GetEntries());
			

			t2->GetEntry(var);
		}
		if(debug) cout << "This is the match found: ev "<< _nevent  <<"\t time_charge: " << Form("%9.3f",time_charge) << "\t time_light: "<< Form("%9.3f",time_light)<<"\tdt:\t"<<Form("%1.2f",time_light-time_charge)<<  endl;

		//if (time_light-time_charge>2) continue;

		//if (time_light-time_charge>0.02) { cout << "WARNING! EVENTS TOO FAR AWAY" << endl; cout << "This is the math found for ev "<< ev  <<"\t time_charge: " << Form("%9.3f",time_charge) <<"\t light event: "<< var << "\t time_light: "<< Form("%9.3f",time_light)<<  endl;lets_pause();}
					
		
		//lets_pause();
                _deltatime=time_light-time_charge;
		_run_2=run_number_L;
		_nevent_2=_nevent;
		_nchannels_2=_nchannels;
		_time_sample_2=_time_sample;
		_nsamples_2=_nsamples;

		if(debug) cout << "lets assing arrays. "<<  endl;
		//lets_pause();
		//memcpy(_adc_value_2[0], _adc_value[0], sizeof(_adc_value[0]));
		//memcpy(_adc_value_2[1], _adc_value[1], sizeof(_adc_value[1]));
		//memcpy(_adc_value_2[2], _adc_value[2], sizeof(_adc_value[2]));
		//memcpy(_adc_value_2[3], _adc_value[3], sizeof(_adc_value[3]));
		//memcpy(_adc_value_2[4], _adc_value[4], sizeof(_adc_value[4]));
		std::copy(std::begin(_adc_value[0]), std::end(_adc_value[0]), _adc_value_2[0]);
		std::copy(std::begin(_adc_value[1]), std::end(_adc_value[1]), _adc_value_2[1]);
		std::copy(std::begin(_adc_value[2]), std::end(_adc_value[2]), _adc_value_2[2]);
		std::copy(std::begin(_adc_value[3]), std::end(_adc_value[3]), _adc_value_2[3]);
		std::copy(std::begin(_adc_value[4]), std::end(_adc_value[4]), _adc_value_2[4]);

		std::copy(std::begin(_crt_track_pos0), std::end(_crt_track_pos0), _crt_track_pos0_2);
		std::copy(std::begin(_crt_track_pos1), std::end(_crt_track_pos1), _crt_track_pos1_2);
		_crt_ToF_2=_crt_ToF;

		std::copy(std::begin(_crt_adc[0]), std::end(_crt_adc[0]), _crt_adc2[0]);
		std::copy(std::begin(_crt_adc[1]), std::end(_crt_adc[1]), _crt_adc2[1]);
		std::copy(std::begin(_crt_adc[2]), std::end(_crt_adc[2]), _crt_adc2[2]);
		std::copy(std::begin(_crt_adc[3]), std::end(_crt_adc[3]), _crt_adc2[3]);

//		memcpy(_adc_value_2[4], _adc_value[4], sizeof(_adc_value[4]));
		memcpy(_PCTimeTag_2, _PCTimeTag, sizeof(_PCTimeTag_2));
		
		_crt_daq_match_2=_crt_daq_match;
		_crt_reco_2=_crt_reco;

		if(debug) cout << "Variables assigned "<<  endl;
		//lets_pause();
                _b2_deltatime->Fill();
		_b2_run->Fill();
		_b2_light_event->Fill();
		_b2_nchannels->Fill();
		_b2_time_sample->Fill();
		_b2_nsamples->Fill();
		_b2_adc_value_0->Fill();
		_b2_adc_value_1->Fill();
		_b2_adc_value_2->Fill();
		_b2_adc_value_3->Fill();
		_b2_adc_value_4->Fill();
		_b2_crt_adc->Fill();
		_b2_PCTimeTag->Fill();
		_b2_crt_daq_match->Fill();
		_b2_crt_reco->Fill();
		_b2_crt_ToF->Fill();
		_b2_crt_track_pos0->Fill();
		_b2_crt_track_pos1->Fill();
        	_b2_track_hits_id->Fill();
		_b2_crt_track_param->Fill();
		_b2_crt_isFV->Fill();
		_b2_crt_track_lenFV->Fill();
		_b2_crt_point_in->Fill();
		_b2_crt_point_out->Fill();
		_b2_crt_pmt_dist->Fill();

		//if ( ev%1000==0 )cout << "Processing run "<<run_number<<". Evaluated " << ev << " events of " <<numentries<< " ("<< 100*(float)ev/(float)numentries<<"%)"<<  endl;


	} // End event loop

	cout << "TTree readed."<< endl;

	fin->Close();
	fin2->cd("analysistree/");
	t2->Write();
	fin2->Close();
	//fout->Close();

	cout << "End of macro."<<  endl;

		
}






int assing_light_run(string charge_run,int runnum_C,int LRUNS, int light_runs[MAXLRUNS],int Lstarttime[MAXLRUNS],int Lendtime[MAXLRUNS])
{


	int Cstarttime, Cendtime;
	timing_charge(charge_run,runnum_C,Cstarttime,Cendtime,"V");


	int j=0;
	bool found=false;
	do{

		//cout << "Let's compare ("<<i<<","<<j<<") C"<< charge_runs[i] << " with L" << light_runs[j] <<  endl;
		//cout << "   timing C("<<Cstarttime[i]<<","<<Cendtime[i]<<")  L("<<Lstarttime[j]<<","<<Lendtime[j]<<")"<<  endl;
		if (Lstarttime[j]<Cstarttime && Lendtime[j]>Cendtime) {cout << "CORRESPONDENCE FOUND! C"<< charge_run << " L" << light_runs[j] <<  endl; found=true; return light_runs[j];}
		j++;
	}while(j<LRUNS && !found);
	if (!found){ cout << "Not Light run found for this charge run C" << charge_run << endl; return -1;}

	return -1;

}



void timing_charge(string charge_run, int run_number, int &starttime, int &endtime,string option){
	cout << "--- CHARGE ana for RUN "  << run_number << " ----"<<  endl;	//
	//
	// Create the TChain to Display
	//
	TChain * t = new TChain("analysistree/anatree");

	//
	// Add root file to Display
	//
	t->Add( Form("%s%i%s%s",indirC.c_str(),run_number,"/",charge_run.c_str()) );
//cout << "t> "<<t << endl;
	//
	// Branches
	//
	int _event;
    	int _eventtime_seconds;
    	int _eventtime_nanoseconds;
    	int _no_hits;

	TBranch * _b_event;	  		t->SetBranchAddress("EventNumberInRun"       		, &_event       		, &_b_event       );
	TBranch * _b_eventtime_seconds;		t->SetBranchAddress("EventTimeSeconds"		, &_eventtime_seconds		, &_b_eventtime_seconds       );
	TBranch * _b_eventtime_nanoseconds;	t->SetBranchAddress("EventTimeNanoseconds"	, &_eventtime_nanoseconds       , &_b_eventtime_nanoseconds       );
	TBranch * _b_no_hits;			t->SetBranchAddress("NumberOfHits"		, &_no_hits			, &_b_no_hits       );

	double time_charge=0.0;
	double time_charge_old=0.0;

	int ev=0;
	t->GetEntry(ev);	
	cout << _event << "\t" << _eventtime_seconds << "\t" << _eventtime_nanoseconds << "\t" << _no_hits <<  endl;
	time_charge=(double)_eventtime_seconds-35.0+(double)_eventtime_nanoseconds*1e-9;
	time_charge_old=time_charge;
	starttime=_eventtime_seconds;

	t->GetEntry(t->GetEntries()-1);
	cout << _event << "\t" << _eventtime_seconds << "\t" << _eventtime_nanoseconds << "\t" << _no_hits <<  endl;
	endtime=_eventtime_seconds;

}


void timing_light(int run_number, int &starttime, int &endtime,string option){


	cout << "--- LIGHT ana for RUN "  << run_number << " ----"<<  endl;	//
	TChain * t = new TChain("midas_data");

	t->Add( TString::Format("%s%s%08d%s",indirL.c_str(),"output",run_number,".root") );

	short _adc_value[5][300000];
	int _PCTimeTag[3], _crt_SE_time[2];
	int _nsamples;
	int _time_stamp;
    	int _nevent;
    	int _nchannels;
    	int _time_sample;
	int _crt_daq_match;
	TBranch * _b_time_stamp;  	t->SetBranchAddress("TimeStamp"   , &_time_stamp   , &_b_time_stamp  );
	TBranch * _b_event;	  	t->SetBranchAddress("event"       , &_nevent       , &_b_event       );
	TBranch * _b_nchannels;   	t->SetBranchAddress("nchannels"   , &_nchannels    , &_b_nchannels   );
	TBranch * _b_time_sample; 	t->SetBranchAddress("TimeSample"  , &_time_sample  , &_b_time_sample );
	TBranch * _b_nsamples;	  	t->SetBranchAddress("nsamples"    , &_nsamples     , &_b_nsamples    );	
	TBranch * _b_adc_value_0; 	t->SetBranchAddress("adc_value_0" , _adc_value[0]  , &_b_adc_value_0 );
	TBranch * _b_adc_value_1; 	t->SetBranchAddress("adc_value_1" , _adc_value[1]  , &_b_adc_value_1 );
	TBranch * _b_adc_value_2; 	t->SetBranchAddress("adc_value_2" , _adc_value[2]  , &_b_adc_value_2 );
	TBranch * _b_adc_value_3; 	t->SetBranchAddress("adc_value_3" , _adc_value[3]  , &_b_adc_value_3 );
	TBranch * _b_adc_value_4; 	t->SetBranchAddress("adc_value_4" , _adc_value[4]  , &_b_adc_value_4 );
	TBranch * _b_PCTimeTag;   	t->SetBranchAddress("PCTimeTag"   , _PCTimeTag     , &_b_PCTimeTag   );
	TBranch * _b_crt_daq_match;  	t->SetBranchAddress("crt_daq_match"	, &_crt_daq_match	, &_b_crt_daq_match   );
	TBranch * _b_crt_SE_time; 	t->SetBranchAddress("crt_SE_time"	, _crt_SE_time		, &_b_crt_SE_time   );

	double time_light=0.0;
	double time_light_old=0.0;
	time_light=(double)_PCTimeTag[0]+0.001* _PCTimeTag[2];
	time_light_old=time_light;

	int ev=0;
	t->GetEntry(ev);	
	cout << run_number << " "<<_nevent << "\t" << _PCTimeTag[0] << "\t" << _PCTimeTag[2] << endl;
	starttime=_time_stamp;

	t->GetEntry(t->GetEntries()-1);
	cout << run_number << " " << _nevent << "\t" << _PCTimeTag[0] << "\t" << _PCTimeTag[2] << endl;
	endtime=_time_stamp;// lets_pause();
}

void lets_pause (){
	TTimer * timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
	timer->TurnOn();
	timer->Reset();
	std::cout << "q/Q to quit, other to continuee: ";
	char kkey;
	std::cin.get(kkey);
	if( kkey == 'q' || kkey == 'Q') 			gSystem->Exit(0);
	timer->TurnOff();
	delete timer;
}

