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
static const int MAXSAMPLES = 300000;

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
	//debug=false;

	std::system(Form("%s%s%i%s%s%s%s","cp -v ",indirC.c_str(),runnum_C,"/",ifile.c_str()," ",tempfolder.c_str()));



	TFile *fin=TFile::Open(TString::Format("%s%s%08d%s",indirL.c_str(),"output",run_number_L,"_reprocessed.root"),"READ"); //PMTrun file
	
	if (!fin || fin->IsZombie()) {
	   cout <<"The file could not be opened!";
	   delete fin;
	} 
	TTree *t = (TTree*)fin->Get("midas_data");
	

	//actual branches at the ntuple:

	Int_t    _TimeStamp;
	//Int_t    _TriggerTimeTag;
	Int_t    _PCTimeTag[3];
	Int_t    _nevent;
	Int_t    _nchannels;
	Int_t    _nsamples;
	Int_t    _TimeSample;
	//Int_t    _channel_nsat[8];
	//Float_t  _channel_pedestal[8];
	//Float_t  _channel_charge[8];
	Int_t    _crt_daq_match;
	Int_t    _crt_reco;
	Int_t    _crt_reco_sat;
	//Int_t    _crt_adc[4][32];
	Int_t    _crt_SE_time[2];
	//Int_t    _crt_TS_time[2];
	//Int_t    _crt_nHitAllTot;
	//Int_t    _crt_nHitTot;
	//Int_t    _crt_nHitSatTot;
	//Float_t  _crt_charge;
	Float_t  _crt_ToF;
	//Int_t    _crt_track_hits_id[4];
	
	Float_t  _crt_track_param[5];
	Float_t  _crt_track_param_p[5];
	Float_t  _crt_track_param_m[5];
	Float_t  _crt_track_door[3];
	Float_t  _crt_track_door_p[3];
	Float_t  _crt_track_door_m[3];
	Float_t  _crt_track_wall[3];
	Float_t  _crt_track_wall_p[3];
	Float_t  _crt_track_wall_m[3];
	Float_t  _crt_pmt_dist[5];
	Float_t  _crt_pmt_dist_p[5];
	Float_t  _crt_pmt_dist_m[5];
	Float_t  _crt_closest_coord[5][3];
	Int_t    _crt_isclosestpoint[5];
	Float_t  _crt_drift_len[5];
	Float_t  _crt_drift_len_p[5];
	Float_t  _crt_drift_len_m[5];
	Int_t    _crt_isFV;
	Float_t  _crt_point_door_fv[3];
	Float_t  _crt_point_door_fv_p[3];
	Float_t  _crt_point_door_fv_m[3];
	Float_t  _crt_point_wall_fv[3];
	Float_t  _crt_point_wall_fv_p[3];
	Float_t  _crt_point_wall_fv_m[3];
	Int_t    _crt_isFC;
	Float_t  _crt_point_door_fc[3];
	Float_t  _crt_point_door_fc_p[3];
	Float_t  _crt_point_door_fc_m[3];
	Float_t  _crt_point_wall_fc[3];
	Float_t  _crt_point_wall_fc_p[3];
	Float_t  _crt_point_wall_fc_m[3];

	Short_t  _adc_value[5][300000];

	TBranch        *b_TimeStamp; t->SetBranchAddress("TimeStamp", &_TimeStamp, &b_TimeStamp);   
	//TBranch        *b_TriggerTimeTag; t->SetBranchAddress("TriggerTimeTag", &TriggerTimeTag, &b_TriggerTimeTag);   
	TBranch        *b_PCTimeTag; t->SetBranchAddress("PCTimeTag", _PCTimeTag, &b_PCTimeTag);   
	TBranch        *b_nevent; t->SetBranchAddress("event", &_nevent, &b_nevent); 
	TBranch        *b_nchannels; t->SetBranchAddress("nchannels", &_nchannels, &b_nchannels);    
	TBranch        *b_nsamples;  t->SetBranchAddress("nsamples", &_nsamples, &b_nsamples);  
	TBranch        *b_TimeSample; t->SetBranchAddress("TimeSample", &_TimeSample, &b_TimeSample);   
	//TBranch        *b_channel_nsat; t->SetBranchAddress("channel_nsat", channel_nsat, &b_channel_nsat);
	//TBranch        *b_channel_pedestal; t->SetBranchAddress("channel_pedestal", channel_pedestal, &b_channel_pedestal);   
	//TBranch        *b_channel_charge; t->SetBranchAddress("channel_charge", channel_charge, &b_channel_charge);  
	TBranch        *b_crt_daq_match; t->SetBranchAddress("crt_daq_match", &_crt_daq_match, &b_crt_daq_match);   
	TBranch        *b_crt_reco; t->SetBranchAddress("crt_reco", &_crt_reco, &b_crt_reco); 
	TBranch        *b_crt_reco_sat; t->SetBranchAddress("crt_reco_sat", &_crt_reco_sat, &b_crt_reco_sat);   
	//TBranch        *b_crt_adc;  t->SetBranchAddress("crt_adc", crt_adc, &b_crt_adc);  
	TBranch        *b_crt_SE_time; t->SetBranchAddress("crt_SE_time", _crt_SE_time, &b_crt_SE_time);   
	//TBranch        *b_crt_TS_time; t->SetBranchAddress("crt_TS_time", crt_TS_time, &b_crt_TS_time);    
	//TBranch        *b_crt_nHitAllTot; t->SetBranchAddress("crt_nHitAllTot", &_crt_nHitAllTot, &b_crt_nHitAllTot);  
	//TBranch        *b_crt_nHitTot; t->SetBranchAddress("crt_nHitTot", &_crt_nHitTot, &b_crt_nHitTot);    
	//TBranch        *b_crt_nHitTot; t->SetBranchAddress("crt_nHitSatTot", &_crt_nHitSatTot, &b_crt_nHitTot);   
	//TBranch        *b_crt_charge; t->SetBranchAddress("crt_charge", &_crt_charge, &b_crt_charge);    
	TBranch        *b_crt_ToF; t->SetBranchAddress("crt_ToF", &_crt_ToF, &b_crt_ToF);  
	//TBranch        *b_crt_track_hits_id; t->SetBranchAddress("crt_track_hits_id", crt_track_hits_id, &b_crt_track_hits_id);

	TBranch        *b_crt_track_param; t->SetBranchAddress("crt_track_param", _crt_track_param, &b_crt_track_param);   
	TBranch        *b_crt_track_param_p; t->SetBranchAddress("crt_track_param_p", _crt_track_param_p, &b_crt_track_param_p);   
	TBranch        *b_crt_track_param_m;  t->SetBranchAddress("crt_track_param_m", _crt_track_param_m, &b_crt_track_param_m);
	TBranch        *b_crt_track_door; t->SetBranchAddress("crt_track_door", _crt_track_door, &b_crt_track_door);    
	TBranch        *b_crt_track_door_p; t->SetBranchAddress("crt_track_door_p", _crt_track_door_p, &b_crt_track_door_p);   
	TBranch        *b_crt_track_door_m;   t->SetBranchAddress("crt_track_door_m", _crt_track_door_m, &b_crt_track_door_m);
	TBranch        *b_crt_track_wall;    t->SetBranchAddress("crt_track_wall", _crt_track_wall, &b_crt_track_wall);
	TBranch        *b_crt_track_wall_p;  t->SetBranchAddress("crt_track_wall_p", _crt_track_wall_p, &b_crt_track_wall_p);  
	TBranch        *b_crt_track_wall_m;  t->SetBranchAddress("crt_track_wall_m", _crt_track_wall_m, &b_crt_track_wall_m);  
	TBranch        *b_crt_pmt_dist;   t->SetBranchAddress("crt_pmt_dist", _crt_pmt_dist, &b_crt_pmt_dist); 
	TBranch        *b_crt_pmt_dist_p;   t->SetBranchAddress("crt_pmt_dist_p", _crt_pmt_dist_p, &b_crt_pmt_dist_p); 
	TBranch        *b_crt_pmt_dist_m;  t->SetBranchAddress("crt_pmt_dist_m", _crt_pmt_dist_m, &b_crt_pmt_dist_m);   
	TBranch        *b_crt_cloosest_coord; t->SetBranchAddress("crt_closest_coord", _crt_closest_coord, &b_crt_cloosest_coord);   
	TBranch        *b_isclosestpoint; t->SetBranchAddress("crt_isclosestpoint", _crt_isclosestpoint, &b_isclosestpoint);   
	TBranch        *b_crt_drift_len;   t->SetBranchAddress("crt_drift_len", _crt_drift_len, &b_crt_drift_len); 
	TBranch        *b_crt_drift_len_p; t->SetBranchAddress("crt_drift_len_p", _crt_drift_len_p, &b_crt_drift_len_p);   
	TBranch        *b_crt_drift_len_m;   t->SetBranchAddress("crt_drift_len_m", _crt_drift_len_m, &b_crt_drift_len_m);
	TBranch        *b_crt_isFV; t->SetBranchAddress("crt_isFV", &_crt_isFV, &b_crt_isFV);

	TBranch        *b_crt_point_door_fv;  t->SetBranchAddress("crt_point_door_fv", _crt_point_door_fv, &b_crt_point_door_fv);  
	TBranch        *b_crt_point_door_fv_p;   t->SetBranchAddress("crt_point_door_fv_p", _crt_point_door_fv_p, &b_crt_point_door_fv_p);
	TBranch        *b_crt_point_door_fv_m;  t->SetBranchAddress("crt_point_door_fv_m", _crt_point_door_fv_m, &b_crt_point_door_fv_m);  
	TBranch        *b_crt_point_wall_fv;  t->SetBranchAddress("crt_point_wall_fv", _crt_point_wall_fv, &b_crt_point_wall_fv);  
	TBranch        *b_crt_point_wall_fv_p;   t->SetBranchAddress("crt_point_wall_fv_p", _crt_point_wall_fv_p, &b_crt_point_wall_fv_p); 
	TBranch        *b_crt_point_wall_fv_m;    t->SetBranchAddress("crt_point_wall_fv_m", _crt_point_wall_fv_m, &b_crt_point_wall_fv_m);
	TBranch        *b_isFC;    t->SetBranchAddress("crt_isFC", &_crt_isFC, &b_isFC);
	TBranch        *b_crt_point_door_fc;  t->SetBranchAddress("crt_point_door_fc", _crt_point_door_fc, &b_crt_point_door_fc);  
	TBranch        *b_crt_point_door_fc_p;  t->SetBranchAddress("crt_point_door_fc_p", _crt_point_door_fc_p, &b_crt_point_door_fc_p);  
	TBranch        *b_crt_point_door_fc_m;  t->SetBranchAddress("crt_point_door_fc_m", _crt_point_door_fc_m, &b_crt_point_door_fc_m);  
	TBranch        *b_crt_point_wall_fc;    t->SetBranchAddress("crt_point_wall_fc", _crt_point_wall_fc, &b_crt_point_wall_fc);
	TBranch        *b_crt_point_wall_fc_p;    t->SetBranchAddress("crt_point_wall_fc_p", _crt_point_wall_fc_p, &b_crt_point_wall_fc_p);
	TBranch        *b_crt_point_wall_fc_m;   t->SetBranchAddress("crt_point_wall_fc_m", _crt_point_wall_fc_m, &b_crt_point_wall_fc_m); 
	TBranch        *b_adc_value_0;   t->SetBranchAddress("adc_value_0", _adc_value[0], &b_adc_value_0); 
	TBranch        *b_adc_value_1; t->SetBranchAddress("adc_value_1", _adc_value[1], &b_adc_value_1);   
	TBranch        *b_adc_value_2;  t->SetBranchAddress("adc_value_2", _adc_value[2], &b_adc_value_2);  
	TBranch        *b_adc_value_3; t->SetBranchAddress("adc_value_3", _adc_value[3], &b_adc_value_3);   
	TBranch        *b_adc_value_4; t->SetBranchAddress("adc_value_4", _adc_value[4], &b_adc_value_4);  


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

	// for safety
	int _run_2;
        float _deltatime;

	TBranch * _b2_deltatime =       t2->Branch("deltatime"		, &_deltatime	   	, "deltatime/F"       			);
	TBranch * _b2_run = 		t2->Branch("light_run"		, &_run_2	   	, "light_run/I"       			);
	TBranch * _b2_light_event = 	t2->Branch("light_event"	, &_nevent		, "light_event/I"     			);
	

	TBranch * _b2_TimeStamp = t2->Branch("TimeStamp", &_TimeStamp, "TimeStamp/I");
	TBranch * _b2_PCTimeTag = t2->Branch("PCTimeTag", _PCTimeTag, "PCTimeTag[3]/I");
	TBranch * _b2_nchannels = t2->Branch("nchannels", &_nchannels, "nchannels/I");
	TBranch * _b2_nsamples = t2->Branch("nsamples", &_nsamples, "nsamples/I");
	TBranch * _b2_TimeSample = t2->Branch("TimeSample", &_TimeSample, "TimeSample/I");
	TBranch * _b2_crt_daq_match = t2->Branch("crt_daq_match", &_crt_daq_match, "crt_daq_match/I");
	TBranch * _b2_crt_reco = t2->Branch("crt_reco", &_crt_reco, "crt_reco/I");
	TBranch * _b2_crt_reco_sat = t2->Branch("crt_reco_sat", &_crt_reco_sat, "crt_reco_sat/I");
	TBranch * _b2_crt_SE_time = t2->Branch("crt_SE_time", _crt_SE_time, "crt_SE_time[2]/I");
	TBranch * _b2_crt_ToF = t2->Branch("crt_ToF", &_crt_ToF, "crt_ToF/F");

	TBranch * _b2_crt_track_param = t2->Branch("crt_track_param", _crt_track_param, "crt_track_param[5]/F");
	TBranch * _b2_crt_track_param_p = t2->Branch("crt_track_param_p", _crt_track_param_p, "crt_track_param_p[5]/F");
	TBranch * _b2_crt_track_param_m = t2->Branch("crt_track_param_m", _crt_track_param_m, "crt_track_param_m[5]/F");
	TBranch * _b2_crt_track_door = t2->Branch("crt_track_door", _crt_track_door, "crt_track_door[3]/F");
	TBranch * _b2_crt_track_door_p = t2->Branch("crt_track_door_p", _crt_track_door_p, "crt_track_door_p[3]/F");
	TBranch * _b2_crt_track_door_m = t2->Branch("crt_track_door_m", _crt_track_door_m, "crt_track_door_m[3]/F");
	TBranch * _b2_crt_track_wall = t2->Branch("crt_track_wall", _crt_track_wall, "crt_track_wall[3]/F");
	TBranch * _b2_crt_track_wall_p = t2->Branch("crt_track_wall_p", _crt_track_wall_p, "crt_track_wall_p[3]/F");
	TBranch * _b2_crt_track_wall_m = t2->Branch("crt_track_wall_m", _crt_track_wall_m, "crt_track_wall_m[3]/F");
	TBranch * _b2_crt_pmt_dist = t2->Branch("crt_pmt_dist", _crt_pmt_dist, "crt_pmt_dist[5]/F");
	TBranch * _b2_crt_pmt_dist_p = t2->Branch("crt_pmt_dist_p", _crt_pmt_dist_p, "crt_pmt_dist_p[5]/F");
	TBranch * _b2_crt_pmt_dist_m = t2->Branch("crt_pmt_dist_m", _crt_pmt_dist_m, "crt_pmt_dist_m[5]/F");
	TBranch * _b2_crt_closest_coord = t2->Branch("crt_closest_coord", _crt_closest_coord, "crt_cloosest_coord[5][3]/F");
	TBranch * _b2_crt_isclosestpoint = t2->Branch("crt_isclosestpoint", _crt_isclosestpoint, "crt_isclosestpoint[5]/I");
	TBranch * _b2_crt_drift_len = t2->Branch("crt_drift_len", _crt_drift_len, "crt_drift_len[5]/F");
	TBranch * _b2_crt_drift_len_p = t2->Branch("crt_drift_len_p", _crt_drift_len_p, "crt_drift_len_p[5]/F");
	TBranch * _b2_crt_drift_len_m = t2->Branch("crt_drift_len_m", _crt_drift_len_m, "crt_drift_len_m[5]/F");
	TBranch * _b2_crt_isFV = t2->Branch("crt_isFV", &_crt_isFV, "crt_isFV/I");
	TBranch * _b2_crt_point_door_fv = t2->Branch("crt_point_door_fv", _crt_point_door_fv, "crt_point_door_fv[3]/F");
	TBranch * _b2_crt_point_door_fv_p = t2->Branch("crt_point_door_fv_p", _crt_point_door_fv_p, "crt_point_door_fv_p[3]/F");
	TBranch * _b2_crt_point_door_fv_m = t2->Branch("crt_point_door_fv_m", _crt_point_door_fv_m, "crt_point_door_fv_m[3]/F");
	TBranch * _b2_crt_point_wall_fv = t2->Branch("crt_point_wall_fv", _crt_point_wall_fv, "crt_point_wall_fv[3]/F");
	TBranch * _b2_crt_point_wall_fv_p = t2->Branch("crt_point_wall_fv_p", _crt_point_wall_fv_p, "crt_point_wall_fv_p[3]/F");
	TBranch * _b2_crt_point_wall_fv_m = t2->Branch("crt_point_wall_fv_m", _crt_point_wall_fv_m, "crt_point_wall_fv_m[3]/F");
	TBranch * _b2_crt_isFC = t2->Branch("crt_isFC", &_crt_isFC, "isFC/I");
	TBranch * _b2_crt_point_door_fc = t2->Branch("crt_point_door_fc", _crt_point_door_fc, "crt_point_door_fc[3]/F");
	TBranch * _b2_crt_point_door_fc_p = t2->Branch("crt_point_door_fc_p", _crt_point_door_fc_p, "crt_point_door_fc_p[3]/F");
	TBranch * _b2_crt_point_door_fc_m = t2->Branch("crt_point_door_fc_m", _crt_point_door_fc_m, "crt_point_door_fc_m[3]/F");
	TBranch * _b2_crt_point_wall_fc = t2->Branch("crt_point_wall_fc", _crt_point_wall_fc, "crt_point_wall_fc[3]/F");
	TBranch * _b2_crt_point_wall_fc_p = t2->Branch("crt_point_wall_fc_p", _crt_point_wall_fc_p, "crt_point_wall_fc_p[3]/F");
	TBranch * _b2_crt_point_wall_fc_m = t2->Branch("crt_point_wall_fc_m", _crt_point_wall_fc_m, "crt_point_wall_fc_m[3]/F");
	TBranch * _b2_adc_value_0 = t2->Branch("adc_value_0", _adc_value[0], "adc_value_0[nsamples]/S");
	TBranch * _b2_adc_value_1 = t2->Branch("adc_value_1", _adc_value[1], "adc_value_1[nsamples]/S");
	TBranch * _b2_adc_value_2 = t2->Branch("adc_value_2", _adc_value[2], "adc_value_2[nsamples]/S");
	TBranch * _b2_adc_value_3 = t2->Branch("adc_value_3", _adc_value[3], "adc_value_3[nsamples]/S");
	TBranch * _b2_adc_value_4 = t2->Branch("adc_value_4", _adc_value[4], "adc_value_4[nsamples]/S");



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


		if(debug) cout << "Variables assigned "<<  endl;
		//cout << "_nevent = " << _nevent << endl;
		//lets_pause();

		_b2_run->Fill();
		_b2_deltatime->Fill();
		_b2_TimeStamp->Fill();      
		_b2_PCTimeTag->Fill();    
		_b2_light_event->Fill();    
		_b2_nchannels->Fill();    
		_b2_nsamples->Fill();    
		_b2_TimeSample->Fill();       
		_b2_crt_daq_match->Fill();    
		_b2_crt_reco->Fill();    
		_b2_crt_reco_sat->Fill();    
		_b2_crt_SE_time->Fill();      
		_b2_crt_ToF->Fill();    
		_b2_crt_track_param->Fill();    
		_b2_crt_track_param_p->Fill();    
		_b2_crt_track_param_m->Fill();    
		_b2_crt_track_door->Fill();    
		_b2_crt_track_door_p->Fill();    
		_b2_crt_track_door_m->Fill();    
		_b2_crt_track_wall->Fill();    
		_b2_crt_track_wall_p->Fill();    
		_b2_crt_track_wall_m->Fill();    
		_b2_crt_pmt_dist->Fill();    
		_b2_crt_pmt_dist_p->Fill();    
		_b2_crt_pmt_dist_m->Fill();    
		_b2_crt_closest_coord->Fill();    
		_b2_crt_isclosestpoint->Fill();    
		_b2_crt_drift_len->Fill();    
		_b2_crt_drift_len_p->Fill();    
		_b2_crt_drift_len_m->Fill();    
		_b2_crt_isFV->Fill();    
		_b2_crt_point_door_fv->Fill();    
		_b2_crt_point_door_fv_p->Fill();    
		_b2_crt_point_door_fv_m->Fill();    
		_b2_crt_point_wall_fv->Fill();    
		_b2_crt_point_wall_fv_p->Fill();    
		_b2_crt_point_wall_fv_m->Fill();    
		_b2_crt_isFC->Fill();    
		_b2_crt_point_door_fc->Fill();    
		_b2_crt_point_door_fc_p->Fill();    
		_b2_crt_point_door_fc_m->Fill();    
		_b2_crt_point_wall_fc->Fill();    
		_b2_crt_point_wall_fc_p->Fill();    
		_b2_crt_point_wall_fc_m->Fill();    
		_b2_adc_value_0->Fill();    
		_b2_adc_value_1->Fill();    
		_b2_adc_value_2->Fill();    
		_b2_adc_value_3->Fill();    
		_b2_adc_value_4->Fill();   

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

	t->Add( TString::Format("%s%s%08d%s",indirL.c_str(),"output",run_number,"_reprocessed.root") );

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

