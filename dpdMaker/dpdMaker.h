#include "../LIB/WaveformAnalysis.cc"
#include "../config_reco.h"
//#include "THDet.C"

bool debug=false;
bool display_tracks=false;
bool display_waveforms=false;

double ADC_to_volts = 2./4096.;
double charge_e = 1.602E-19;

static const int NMAXPEAKS=50;
static const int KMAXNHITS=100000;
static const int KMAXNTRACKS=1000;

void make_dpd(TChain* t2, int runNum, double gains[N_PMT], string outfilename);

void make_dpd(TChain* t2, int runNum, double gains[N_PMT], string outfilename)
{
	bool doCharge = (TString(outfilename).Contains("matched"))?true:false;
	
	cout << "doCharge = " << doCharge << ", runNum = " << runNum << endl;
	
	TFile ofile(outfilename.c_str(),"RECREATE");
	if (ofile.IsZombie())
	{
	       cout << "Error opening output file " << outfilename.c_str() << endl;
	       exit(-1);
	}

	TTree *dpd= new TTree("dpd","result tree");
	
	
	//#include "../LIB/branches_2018Feb05.h"
	
	// Input charge variables
	int _runcharge;
	int _subruncharge;
	int _event;
    int _eventtime_seconds;
    int _eventtime_nanoseconds;
	
	int _no_hits;
	float _Hit_ChargeIntegral[KMAXNHITS];
	short _hit_trkid[KMAXNHITS];
	short _NumberOfTracks_pmtrack;
	short _nclusters;

	float _Track_PitchInViews_pmtrack[KMAXNTRACKS][2];
	short _Track_NumberOfHitsPerView_pmtrack[KMAXNTRACKS][2];

	float _Track_TPC_pmtrack[KMAXNHITS];
	float _Track_HitX_pmtrack[KMAXNHITS];
	float _Track_HitY_pmtrack[KMAXNHITS];
	float _Track_HitZ_pmtrack[KMAXNHITS];
	float _Track_StartX_pmtrack[KMAXNTRACKS];
	float _Track_StartY_pmtrack[KMAXNTRACKS];
	float _Track_StartZ_pmtrack[KMAXNTRACKS];
	
	float _Track_Hit_ChargeSummedADC_pmtrack[KMAXNHITS]; 
	float _Track_Hit_ChargeIntegral_pmtrack[KMAXNHITS];
	float _Track_Hit_dx_LocalTrackDirection_pmtrack[KMAXNHITS];
	float _Track_Hit_dx_3DPosition_pmtrack[KMAXNHITS];

	float _Track_EndX_pmtrack[KMAXNTRACKS];
	float _Track_EndY_pmtrack[KMAXNTRACKS];
	float _Track_EndZ_pmtrack[KMAXNTRACKS];
	float _Track_StartDirection_Theta_pmtrack[KMAXNTRACKS];
	float _Track_Length_pmtrack[KMAXNTRACKS];
	float _Track_StartDirection_Phi_pmtrack[KMAXNTRACKS];
	float _Hit_PeakTime[KMAXNHITS];
	short _Hit_Channel[KMAXNHITS];
	short _Hit_View[KMAXNHITS];
	short _Hit_TPC[KMAXNHITS];
	
	// Input light Variables
	int _PCTimeTag[3];
	int _nsamples;
	int _time_stamp;
    int _nevent;
    int _nchannels;
  	int _crt_daq_match;
  	int _crt_reco;
    int _time_sample;
	int _runlight;
	short _adc_value[5][300000];
	int	_crt_adc[4][32];
	float	_crt_track_pos0[3];
	float	_crt_track_pos1[3];
	float	_crt_ToF; 
	int		_crt_isFV;
	float	_crt_track_lenFV;
	float	_crt_point_in[3];
	float _crt_point_out[3];
	float _crt_pmt_dist[3];
	

	if (doCharge)
	{
		 // Input charge + light branches
		
		TBranch * _b2_runcharge;	  	t2->SetBranchAddress("Run"       		, &_runcharge       		, &_b2_runcharge       );
		TBranch * _b2_subruncharge;	  	t2->SetBranchAddress("Subrun"       		, &_subruncharge       		, &_b2_subruncharge       );
		TBranch * _b2_event;	  		t2->SetBranchAddress("EventNumberInRun"       		, &_event       		, &_b2_event       );
		TBranch * _b2_eventtime_seconds;	t2->SetBranchAddress("EventTimeSeconds"		, &_eventtime_seconds		, &_b2_eventtime_seconds       );
		TBranch * _b2_eventtime_nanoseconds;	t2->SetBranchAddress("EventTimeNanoseconds"	, &_eventtime_nanoseconds       , &_b2_eventtime_nanoseconds       );

	
		TBranch * _b2_no_hits;			t2->SetBranchAddress("NumberOfHits"		, &_no_hits			, &_b2_no_hits       );
		TBranch * _b2_Hit_ChargeIntegral;	t2->SetBranchAddress("Hit_ChargeIntegral"	, _Hit_ChargeIntegral		, &_b2_Hit_ChargeIntegral      );
		TBranch * _b2_hit_trkid;		t2->SetBranchAddress("Hit_TrackID"		, _hit_trkid			, &_b2_hit_trkid      );
		TBranch * _b2_NumberOfTracks_pmtrack;	t2->SetBranchAddress("NumberOfTracks_pmtrack"	, &_NumberOfTracks_pmtrack	, &_b2_NumberOfTracks_pmtrack      );
		TBranch * _b2_nclusters;		t2->SetBranchAddress("NumberOfClusters"		, &_nclusters			, &_b2_nclusters       );
		TBranch * _b2_Track_PitchInViews_pmtrack;t2->SetBranchAddress("Track_PitchInViews_pmtrack"		, _Track_PitchInViews_pmtrack			, &_b2_Track_PitchInViews_pmtrack       );

		TBranch * _b2_Hit_TPC;			t2->SetBranchAddress("Hit_TPC"		, _Hit_TPC		, &_b2_Hit_TPC		);
		TBranch * _b2_Hit_View;			t2->SetBranchAddress("Hit_View"		, _Hit_View		, &_b2_Hit_View		);
		TBranch * _b2_Hit_Channel;		t2->SetBranchAddress("Hit_Channel"	, _Hit_Channel		, &_b2_Hit_Channel	);
		TBranch * _b2_Hit_PeakTime;		t2->SetBranchAddress("Hit_PeakTime"	, _Hit_PeakTime		, &_b2_Hit_PeakTime	);

		TBranch * _b2_Track_NumberOfHitsPerView_pmtrack;t2->SetBranchAddress("Track_NumberOfHitsPerView_pmtrack", _Track_NumberOfHitsPerView_pmtrack, &_b2_Track_NumberOfHitsPerView_pmtrack);
	
		TBranch * _b2_Track_Hit_ChargeSummedADC_pmtrack ;t2->SetBranchAddress("Track_Hit_ChargeSummedADC_pmtrack"		, _Track_Hit_ChargeSummedADC_pmtrack			, &_b2_Track_Hit_ChargeSummedADC_pmtrack       );
		TBranch * _b2_Track_Hit_ChargeIntegral_pmtrack ;t2->SetBranchAddress("Track_Hit_ChargeIntegral_pmtrack"		, _Track_Hit_ChargeIntegral_pmtrack			, &_b2_Track_Hit_ChargeIntegral_pmtrack       );


		TBranch * _b2_Track_Hit_dx_LocalTrackDirection_pmtrack ; t2->SetBranchAddress("Track_Hit_dx_LocalTrackDirection_pmtrack"		, _Track_Hit_dx_LocalTrackDirection_pmtrack			, &_b2_Track_Hit_ChargeIntegral_pmtrack       );
		TBranch * _b2_Track_Hit_dx_3DPosition_pmtrack ;t2->SetBranchAddress("Track_Hit_dx_3DPosition_pmtrack"		, _Track_Hit_dx_3DPosition_pmtrack			, &_b2_Track_Hit_dx_3DPosition_pmtrack       );


		TBranch * _b2_Track_TPC_pmtrack		;t2->SetBranchAddress("Track_Hit_TPC_pmtrack"		, _Track_TPC_pmtrack			, &_b2_Track_TPC_pmtrack       );
		TBranch * _b2_Track_HitX_pmtrack	;t2->SetBranchAddress("Track_Hit_X_pmtrack"		, _Track_HitX_pmtrack			, &_b2_Track_HitX_pmtrack       );
		TBranch * _b2_Track_HitY_pmtrack	;t2->SetBranchAddress("Track_Hit_Y_pmtrack"		, _Track_HitY_pmtrack			, &_b2_Track_HitY_pmtrack       );
		TBranch * _b2_Track_HitZ_pmtrack	;t2->SetBranchAddress("Track_Hit_Z_pmtrack"		, _Track_HitZ_pmtrack			, &_b2_Track_HitZ_pmtrack       );

		TBranch * _b2_Track_StartX_pmtrack	;t2->SetBranchAddress("Track_StartPoint_X_pmtrack"		, _Track_StartX_pmtrack			, &_b2_Track_StartX_pmtrack       );
		TBranch * _b2_Track_StartY_pmtrack	;t2->SetBranchAddress("Track_StartPoint_Y_pmtrack"		, _Track_StartY_pmtrack			, &_b2_Track_StartY_pmtrack       );
		TBranch * _b2_Track_StartZ_pmtrack	;t2->SetBranchAddress("Track_StartPoint_Z_pmtrack"		, _Track_StartZ_pmtrack			, &_b2_Track_StartZ_pmtrack       );

		TBranch * _b2_Track_EndX_pmtrack	;t2->SetBranchAddress("Track_EndPoint_X_pmtrack"		, _Track_EndX_pmtrack			, &_b2_Track_EndX_pmtrack       );
		TBranch * _b2_Track_EndY_pmtrack	;t2->SetBranchAddress("Track_EndPoint_Y_pmtrack"		, _Track_EndY_pmtrack			, &_b2_Track_EndY_pmtrack       );
		TBranch * _b2_Track_EndZ_pmtrack	;t2->SetBranchAddress("Track_EndPoint_Z_pmtrack"		, _Track_EndZ_pmtrack			, &_b2_Track_EndZ_pmtrack       );


		TBranch * _b2_Track_StartDirection_Theta_pmtrack	;t2->SetBranchAddress("Track_StartDirection_Theta_pmtrack"		, _Track_StartDirection_Theta_pmtrack			, &_b2_Track_StartDirection_Theta_pmtrack      );
		TBranch * _b2_Track_Length_pmtrack	;t2->SetBranchAddress("Track_Length_pmtrack"		, _Track_Length_pmtrack			, &_b2_Track_Length_pmtrack       );
		TBranch * _b2_Track_StartDirection_Phi_pmtrack	;	t2->SetBranchAddress("Track_StartDirection_Phi_pmtrack"		, _Track_StartDirection_Phi_pmtrack			, &_b2_Track_StartDirection_Phi_pmtrack       );



		TBranch * _b2_runlight;		t2->SetBranchAddress("light_run"		, &_runlight			, &_b2_runlight       	);
	
		TBranch * _b2_light_event;	t2->SetBranchAddress("light_event"		, &_nevent    		, &_b2_light_event	);
		TBranch * _b2_time_sample;	t2->SetBranchAddress("light_TimeSample"		, &_time_sample		, &_b2_time_sample   	);	
		TBranch * _b2_nchannels;	t2->SetBranchAddress("light_nchannels"		, &_nchannels		, &_b2_nchannels   	);
		TBranch * _b2_PCTimeTag;	t2->SetBranchAddress("light_PCTimeTag"  	, _PCTimeTag  		, &_b2_PCTimeTag	);

		TBranch * _b2_nsamples;		t2->SetBranchAddress("light_nsamples"   	, &_nsamples  		,  &_b2_nsamples	);	
	
		TBranch * _b2_adc_value_0; 	t2->SetBranchAddress("light_adc_value_0"	, _adc_value[0]  	, &_b2_adc_value_0 	);
		TBranch * _b2_adc_value_1; 	t2->SetBranchAddress("light_adc_value_1"	, _adc_value[1]  	, &_b2_adc_value_1 	);
		TBranch * _b2_adc_value_2; 	t2->SetBranchAddress("light_adc_value_2"	, _adc_value[2]  	, &_b2_adc_value_2 	);
		TBranch * _b2_adc_value_3;  	t2->SetBranchAddress("light_adc_value_3"	, _adc_value[3]  	, &_b2_adc_value_3 	);
		TBranch * _b2_adc_value_4;  	t2->SetBranchAddress("light_adc_value_4"	, _adc_value[4]  	, &_b2_adc_value_4 	);

		TBranch * _b2_crt_daq_match;	t2->SetBranchAddress("light_crt_daq_match", &_crt_daq_match	, &_b2_crt_daq_match		);
		TBranch * _b2_crt_reco;		t2->SetBranchAddress("light_crt_reco"	, &_crt_reco		, &_b2_crt_reco		);
	
		TBranch * _b2_crt_adc      	;  t2->SetBranchAddress("crt_adc"       , _crt_adc          , &_b2_crt_adc       );
		TBranch * _b2_crt_track_pos0    ;  t2->SetBranchAddress("crt_track_pos0"       , _crt_track_pos0          , &_b2_crt_track_pos0       );
		TBranch * _b2_crt_track_pos1    ;  t2->SetBranchAddress("crt_track_pos1"       , _crt_track_pos1          , &_b2_crt_track_pos1       );
		TBranch * _b2_crt_ToF      	; t2->SetBranchAddress("crt_ToF"       , &_crt_ToF          , &_b2_crt_ToF       );

		TBranch * _b2_crt_isFV      	;  t2->SetBranchAddress("crt_isFV"        , &_crt_isFV         , &_b2_crt_isFV        	);
		TBranch * _b2_crt_track_lenFV   ;  t2->SetBranchAddress("crt_track_lenFV" , &_crt_track_lenFV  , &_b2_crt_track_lenFV	);
		TBranch * _b2_crt_point_in      ;  t2->SetBranchAddress("crt_point_in"    , _crt_point_in      , &_b2_crt_point_in       );
		TBranch * _b2_crt_point_out     ;  t2->SetBranchAddress("crt_point_out"    , _crt_point_out      , &_b2_crt_point_out    );
		TBranch * _b2_crt_pmt_dist     ;  t2->SetBranchAddress("crt_pmt_dist"    , _crt_pmt_dist      , &_b2_crt_pmt_dist       );
	}
	
	else 
	{ 
		// Input light-only branches
		
		TBranch * _b2_light_event;	t2->SetBranchAddress("event"		, &_nevent    		, &_b2_light_event	);
		TBranch * _b2_time_sample;	t2->SetBranchAddress("TimeSample"		, &_time_sample		, &_b2_time_sample   	);	
		TBranch * _b2_nchannels;	t2->SetBranchAddress("nchannels"		, &_nchannels		, &_b2_nchannels   	);
		TBranch * _b2_PCTimeTag;	t2->SetBranchAddress("PCTimeTag"  	, _PCTimeTag  		, &_b2_PCTimeTag	);
	
		TBranch * _b2_nsamples;		t2->SetBranchAddress("nsamples"   	, &_nsamples  		,  &_b2_nsamples	);	
	
		TBranch * _b2_adc_value_0; 	t2->SetBranchAddress("adc_value_0"	, _adc_value[0]  	, &_b2_adc_value_0 	);
		TBranch * _b2_adc_value_1; 	t2->SetBranchAddress("adc_value_1"	, _adc_value[1]  	, &_b2_adc_value_1 	);
		TBranch * _b2_adc_value_2; 	t2->SetBranchAddress("adc_value_2"	, _adc_value[2]  	, &_b2_adc_value_2 	);
		TBranch * _b2_adc_value_3;  	t2->SetBranchAddress("adc_value_3"	, _adc_value[3]  	, &_b2_adc_value_3 	);
		TBranch * _b2_adc_value_4;  	t2->SetBranchAddress("adc_value_4"	, _adc_value[4]  	, &_b2_adc_value_4 	);

		TBranch * _b2_crt_daq_match;	t2->SetBranchAddress("crt_daq_match", &_crt_daq_match	, &_b2_crt_daq_match		);
		TBranch * _b2_crt_reco;		t2->SetBranchAddress("crt_reco"	, &_crt_reco		, &_b2_crt_reco		);
	
		TBranch * _b2_crt_adc      	;  t2->SetBranchAddress("crt_adc"       , _crt_adc          , &_b2_crt_adc       );
		TBranch * _b2_crt_track_pos0    ;  t2->SetBranchAddress("crt_track_pos0"       , _crt_track_pos0          , &_b2_crt_track_pos0       );
		TBranch * _b2_crt_track_pos1    ;  t2->SetBranchAddress("crt_track_pos1"       , _crt_track_pos1          , &_b2_crt_track_pos1       );
		TBranch * _b2_crt_ToF      	; t2->SetBranchAddress("crt_ToF"       , &_crt_ToF          , &_b2_crt_ToF       );

		TBranch * _b2_crt_isFV      	;  t2->SetBranchAddress("crt_isFV"        , &_crt_isFV         , &_b2_crt_isFV        	);
		TBranch * _b2_crt_track_lenFV   ;  t2->SetBranchAddress("crt_track_lenFV" , &_crt_track_lenFV  , &_b2_crt_track_lenFV	);
		TBranch * _b2_crt_point_in      ;  t2->SetBranchAddress("crt_point_in"    , _crt_point_in      , &_b2_crt_point_in       );
		TBranch * _b2_crt_point_out     ;  t2->SetBranchAddress("crt_point_out"    , _crt_point_out      , &_b2_crt_point_out    );
		TBranch * _b2_crt_pmt_dist     ;  t2->SetBranchAddress("crt_pmt_dist"    , _crt_pmt_dist      , &_b2_crt_pmt_dist       );
	}
	
	// DPD variables

	int _run_charge_out;
	int _subrun_charge_out;
	int _nev_charge_out;
	
	int _run_light_out;
	int _nev_light_out;

	double _time_charge;
	double _time_light;
	
	double _pmt_charge[N_PMT+1];
	double _pmt_npe[N_PMT+1];
	double _pmt_ped[N_PMT];
	double _pmt_pedrms[N_PMT];
	double _pmt_ped2[N_PMT];
	double _pmt_pedrms2[N_PMT];
	double _pmt_ped_end[N_PMT];
	double _pmt_ped_end_corr[N_PMT];
	double _pmt_pedrms_end[N_PMT];
	double _pmt_pedrms_end_corr[N_PMT];
	int	   _pmt_wvf_ADC_sat[N_PMT];
	double _pmt_wvf_end_fit_c0[N_PMT];
	double _pmt_wvf_end_fit_c1[N_PMT];
	double _pmt_wvf_end_fit_c2[N_PMT];
	double _pmt_wvf_end_fit_chi2[N_PMT];
	int    _pmt_wvf_end_fit_ndof[N_PMT];
	double _pmt_wvf_properties[N_PMT][2];
	int    _pmt_npeaks[N_PMT]={0};
	double _pmt_peaks_tau[N_PMT][NMAXPEAKS]={0};

	double _pmt_S1_charge_1us[N_PMT+1];
	double _pmt_S1_charge_4us[N_PMT+1];
	double _pmt_S1_charge_80ns[N_PMT+1];
	double _pmt_S1_charge_m2[N_PMT+1];
	double _pmt_S1_charge_4us_corr[N_PMT+1];
	double _pmt_S1_npe_1us[N_PMT+1];
	double _pmt_S1_npe_4us[N_PMT+1];
	double _pmt_S1_npe_80ns[N_PMT+1];
	double _pmt_S1_npe_m2[N_PMT+1];
	double _pmt_S1_npe_4us_corr[N_PMT+1];
	double _pmt_S1_width[N_PMT];
	double _pmt_S1_amp[N_PMT];
	double _pmt_S1_tau[N_PMT];
	double _pmt_S1_tau_end[N_PMT];

	double _pmt_S2_charge[N_PMT+1];
	double _pmt_S2_charge_m2[N_PMT+1];
	double _pmt_S2_charge_corr[N_PMT+1];
	double _pmt_S2_charge_corr_p[N_PMT+1];
	double _pmt_S2_charge_corr_m[N_PMT+1];
	double _pmt_S2_npe[N_PMT+1];
	double _pmt_S2_npe_m2[N_PMT+1];
	double _pmt_S2_npe_corr[N_PMT+1];
	double _pmt_S2_npe_corr_p[N_PMT+1];
	double _pmt_S2_npe_corr_m[N_PMT+1];
	double _pmt_S2_width[N_PMT];
	double _pmt_S2_amp[N_PMT];
	double _pmt_S2_tau[N_PMT];
	double _pmt_S2_tau_start[N_PMT];
	double _pmt_S2_tau_end[N_PMT];
	double _pmt_S2_tau_avg[N_PMT];

	double _tpc_totcharge;
	double _tpc_totrecocharge;
	double _tpc_track_charge[1000]={0};

	double _crt_mYX;
	double _crt_mZX;
	double _crt_nYX;
	double _crt_nZX;
	double _crt_ToF_n;
	double _crt_track_lenFV_n;
	double _crt_pmt_dist_n[3];
	int    _crt_isFV_n;

	short _ntracks;
	short _n_sel_tracks;
	double _tpc_track_mYX[KMAXNTRACKS]={0};
	double _tpc_track_mZX[KMAXNTRACKS]={0};
	double _tpc_track_nYX[KMAXNTRACKS]={0};
	double _tpc_track_nZX[KMAXNTRACKS]={0};
	double _tpc_track_startX[KMAXNTRACKS];
	double _tpc_track_endX[KMAXNTRACKS];
	
	double _tpc_track_start_theta[KMAXNTRACKS]={0};
	double _tpc_track_start_phi[KMAXNTRACKS]={0};
	
	int _tpc_track_fitresult_yx[KMAXNTRACKS]={-1};
	int _tpc_track_fitresult_zx[KMAXNTRACKS]={-1};
	int _crt_matchreco;
	double _tpc_track_length_n[KMAXNTRACKS]={0};

	double _tpc_drift_time_at_pmt_pos[N_PMT]={0};
	
	//DPD branches
	
	if (doCharge)
	{
		TBranch * _bn_run_charge	= dpd->Branch("run_charge"      , &_run_charge_out         , "run_charge/I"      );
		TBranch * _bn_ev_charge		= dpd->Branch("ev_charge"      , &_nev_charge_out         , "ev_charge/I"      );
		TBranch * _bn_subrun_charge = dpd->Branch("subrun_charge"      , &_subrun_charge_out         , "subrun_charge/I"      );
		TBranch * _bn_time_charge 	= dpd->Branch("time_charge", &_time_charge  , "time_charge/D");	
	}
	
	TBranch * _bn_run_light		= dpd->Branch("run_light"      , &_run_light_out         , "run_light/I"      );
	TBranch * _bn_ev_light		= dpd->Branch("ev_light"      , &_nev_light_out         , "ev_light/I"      );
	TBranch * _bn_nsamples 		= dpd->Branch("nsamples", &_nsamples , "nsamples/I");
	TBranch * _bn_time_light 	= dpd->Branch("time_light", &_time_light  , "time_light/D");
	
	if (doCharge)
	{
		TBranch * _bn_ntracks			= dpd->Branch("ntracks"   	,&_ntracks      	, "ntracks/S"   );
		TBranch * _bn_n_sel_tracks		= dpd->Branch("n_sel_tracks"   	,&_n_sel_tracks      	, "n_sel_tracks/S"   );
		TBranch * _bn_tpc_totcharge 	= dpd->Branch("tpc_totcharge"   , &_tpc_totcharge       , "tpc_totcharge/D"    );	
		TBranch * _bn_tpc_totrecocharge = dpd->Branch("tpc_totrecocharge"   , &_tpc_totrecocharge       , "tpc_totrecocharge/D"    );	
	}
	
	TBranch * _bn_crt_matchreco		= dpd->Branch("crt_matchreco"      , &_crt_matchreco         , "crt_matchreco/I"      );
	TBranch * _bn_crt_mYX			= dpd->Branch("crt_mYX"   	,&_crt_mYX      	, "crt_mYX/D"   );
	TBranch * _bn_crt_mZX			= dpd->Branch("crt_mZX"   	,&_crt_mZX      	, "crt_mZX/D"   );
	TBranch * _bn_crt_nYX			= dpd->Branch("crt_nYX"   	,&_crt_nYX      	, "crt_nYX/D"   );
	TBranch * _bn_crt_nZX			= dpd->Branch("crt_nZX"   	,&_crt_nZX      	, "crt_nZX/D"   );
	TBranch * _bn_crt_ToF			= dpd->Branch("crt_ToF"   	,&_crt_ToF_n      	, "crt_ToF/D"   );
	TBranch * _bn_crt_track_lenFV	= dpd->Branch("crt_track_lenFV"   	,&_crt_track_lenFV_n      	, "crt_track_lenFV/D"   );
	TBranch * _bn_crt_pmt_dist_n	= dpd->Branch("crt_pmt_dist"   	,&_crt_pmt_dist_n      	, "crt_pmt_dist[3]/D"   );
	TBranch * _bn_crt_isFV			= dpd->Branch("crt_isFV"   	,&_crt_isFV_n      	, "crt_isFV/I"   );
	
	TBranch * _bn_pmt_wvf_properties = dpd->Branch("pmt_wvf_properties"    , _pmt_wvf_properties       , "pmt_wvf_properties[5][2]/D"    );

	TBranch * _bn_pmt_charge 		= dpd->Branch("pmt_charge"    , _pmt_charge       , "pmt_charge[6]/D"    );	
	TBranch * _bn_pmt_npe			= dpd->Branch("pmt_npe"   , _pmt_npe      , "pmt_npe[6]/D"   );
	TBranch * _bn_pmt_ped       	= dpd->Branch("pmt_ped"   ,  _pmt_ped      , "pmt_ped[5]/D"   );
	TBranch * _bn_pmt_pedrms    	= dpd->Branch("pmt_pedrms", _pmt_pedrms   , "pmt_pedrms[5]/D"   );
	TBranch * _bn_pmt_ped2       	= dpd->Branch("pmt_ped2"   ,  _pmt_ped2      , "pmt_ped2[5]/D"   );
	TBranch * _bn_pmt_pedrms2    	= dpd->Branch("pmt_pedrms2", _pmt_pedrms2   , "pmt_pedrms2[5]/D"   );
	TBranch * _bn_pmt_ped_end   	= dpd->Branch("pmt_ped_end"   ,  _pmt_ped_end      , "pmt_ped_end[5]/D"   );
	TBranch * _bn_pmt_ped_end_corr   	= dpd->Branch("pmt_ped_end_corr"   ,  _pmt_ped_end_corr      , "pmt_ped_end_corr[5]/D"   );
	TBranch * _bn_pmt_pedrms_end   	= dpd->Branch("pmt_pedrms_end"   ,  _pmt_pedrms_end      , "pmt_pedrms_end[5]/D"   );
	TBranch * _bn_pmt_pedrms_end_corr   	= dpd->Branch("pmt_pedrms_end_corr"   ,  _pmt_pedrms_end_corr      , "pmt_pedrms_end_corr[5]/D"   );
	TBranch * _bn_pmt_wvf_ADC_sat  	= dpd->Branch("pmt_wvf_ADC_sat"   ,  _pmt_wvf_ADC_sat      , "pmt_wvf_ADC_sat[5]/I"   );
	TBranch * _bn_pmt_wvf_end_fit_c0   	= dpd->Branch("pmt_wvf_end_fit_c0"   ,  _pmt_wvf_end_fit_c0      , "pmt_wvf_end_fit_c0[5]/D"   );
	TBranch * _bn_pmt_wvf_end_fit_c1   	= dpd->Branch("pmt_wvf_end_fit_c1"   ,  _pmt_wvf_end_fit_c1      , "pmt_wvf_end_fit_c1[5]/D"   );
	TBranch * _bn_pmt_wvf_end_fit_c2   	= dpd->Branch("pmt_wvf_end_fit_c2"   ,  _pmt_wvf_end_fit_c2      , "pmt_wvf_end_fit_c2[5]/D"   );
	TBranch * _bn_pmt_wvf_end_fit_chi2  = dpd->Branch("pmt_wvf_end_fit_chi2"   ,  _pmt_wvf_end_fit_chi2      , "pmt_wvf_end_fit_chi2[5]/D"   );
	TBranch * _bn_pmt_wvf_end_fit_ndof  = dpd->Branch("pmt_wvf_end_fit_ndof"   ,  _pmt_wvf_end_fit_ndof      , "pmt_wvf_end_fit_ndof[5]/I"   );
	TBranch * _bn_pmt_npeaks	    = dpd->Branch("pmt_npeaks"   ,  _pmt_npeaks    , "pmt_npeaks[5]/I"   );
	TBranch * _bn_pmt_peaks_tau 	= dpd->Branch("pmt_peaks_tau"   , _pmt_peaks_tau   ,  Form("pmt_peaks_tau[5][%d]/D",NMAXPEAKS) );
	
	TBranch * _bn_pmt_S1_charge_1us	= dpd->Branch("pmt_S1_charge_1us"   ,_pmt_S1_charge_1us      , "pmt_S1_charge_1us[6]/D"   );
	TBranch * _bn_pmt_S1_charge_4us	= dpd->Branch("pmt_S1_charge_4us"   ,_pmt_S1_charge_4us      , "pmt_S1_charge_4us[6]/D"   );
	TBranch * _bn_pmt_S1_charge_80ns= dpd->Branch("pmt_S1_charge_80ns"   ,_pmt_S1_charge_80ns      , "pmt_S1_charge_80ns[6]/D"   );
	TBranch * _bn_pmt_S1_charge_m2	= dpd->Branch("pmt_S1_charge_m2"   ,_pmt_S1_charge_m2      , "pmt_S1_charge_m2[6]/D"   );
	TBranch * _bn_pmt_S1_charge_4us_corr	= dpd->Branch("pmt_S1_charge_4us_corr"   ,_pmt_S1_charge_4us_corr      , "pmt_S1_charge_4us_corr[6]/D"   );
	TBranch * _bn_pmt_S1_npe_1us	= dpd->Branch("pmt_S1_npe_1us"   , _pmt_S1_npe_1us      , "pmt_S1_npe_1us[6]/D"   );
	TBranch * _bn_pmt_S1_npe_4us	= dpd->Branch("pmt_S1_npe_4us"   , _pmt_S1_npe_4us      , "pmt_S1_npe_4us[6]/D"   );
	TBranch * _bn_pmt_S1_npe_80ns	= dpd->Branch("pmt_S1_npe_80ns"   , _pmt_S1_npe_80ns      , "pmt_S1_npe_80ns[6]/D"   );
	TBranch * _bn_pmt_S1_npe_m2		= dpd->Branch("pmt_S1_npe_m2"   , _pmt_S1_npe_m2      , "pmt_S1_npe_m2[6]/D"   );
	TBranch * _bn_pmt_S1_npe_4us_corr	= dpd->Branch("pmt_S1_npe_4us_corr"   , _pmt_S1_npe_4us_corr      , "pmt_S1_npe_4us_corr[6]/D"   );
	TBranch * _bn_pmt_S1_width		= dpd->Branch("pmt_S1_width"   , _pmt_S1_width      , "pmt_S1_width[5]/D"   );
	TBranch * _bn_pmt_S1_amp		= dpd->Branch("pmt_S1_amp"   , _pmt_S1_amp      , "pmt_S1_amp[5]/D"   );
	TBranch * _bn_pmt_S1_tau		= dpd->Branch("pmt_S1_tau"   , _pmt_S1_tau      , "pmt_S1_tau[5]/D"   );
	TBranch * _bn_pmt_S1_tau_end	= dpd->Branch("pmt_S1_tau_end"   , _pmt_S1_tau_end      , "pmt_S1_tau_end[5]/D"   );

	TBranch * _bn_pmt_S2_charge		= dpd->Branch("pmt_S2_charge"   ,_pmt_S2_charge      , "pmt_S2_charge[6]/D"   );
	TBranch * _bn_pmt_S2_charge_m2	= dpd->Branch("pmt_S2_charge_m2"   ,_pmt_S2_charge_m2      , "pmt_S2_charge_m2[6]/D"   );
	TBranch * _bn_pmt_S2_charge_corr = dpd->Branch("pmt_S2_charge_corr"   ,_pmt_S2_charge_corr      , "pmt_S2_charge_corr[6]/D"   );
	TBranch * _bn_pmt_S2_charge_corr_p = dpd->Branch("pmt_S2_charge_corr_p"   ,_pmt_S2_charge_corr_p      , "pmt_S2_charge_corr_p[6]/D"   );
	TBranch * _bn_pmt_S2_charge_corr_m = dpd->Branch("pmt_S2_charge_corr_m"   ,_pmt_S2_charge_corr_m      , "pmt_S2_charge_corr_m[6]/D"   );
	TBranch * _bn_pmt_S2_npe		= dpd->Branch("pmt_S2_npe"   , _pmt_S2_npe      , "pmt_S2_npe[6]/D"   );
	TBranch * _bn_pmt_S2_npe_m2		= dpd->Branch("pmt_S2_npe_m2"   , _pmt_S2_npe_m2      , "pmt_S2_npe_m2[6]/D"   );
	TBranch * _bn_pmt_S2_npe_corr	= dpd->Branch("pmt_S2_npe_corr"   , _pmt_S2_npe_corr      , "pmt_S2_npe_corr[6]/D"   );
	TBranch * _bn_pmt_S2_npe_corr_p	= dpd->Branch("pmt_S2_npe_corr_p"   , _pmt_S2_npe_corr_p      , "pmt_S2_npe_corr_p[6]/D"   );
	TBranch * _bn_pmt_S2_npe_corr_m	= dpd->Branch("pmt_S2_npe_corr_m"   , _pmt_S2_npe_corr_m      , "pmt_S2_npe_corr_m[6]/D"   );
	TBranch * _bn_pmt_S2_width		= dpd->Branch("pmt_S2_width"   , _pmt_S2_width      , "pmt_S2_width[5]/D"   );
	TBranch * _bn_pmt_S2_amp		= dpd->Branch("pmt_S2_amp"   , _pmt_S2_amp      , "pmt_S2_amp[5]/D"   );
	TBranch * _bn_pmt_S2_tau		= dpd->Branch("pmt_S2_tau"   , _pmt_S2_tau      , "pmt_S2_tau[5]/D"   );
	TBranch * _bn_pmt_S2_tau_avg 	= dpd->Branch("pmt_S2_tau_avg"   , _pmt_S2_tau_avg      , "pmt_S2_tau_avg[5]/D"   );
	TBranch * _bn_pmt_S2_tau_start	= dpd->Branch("pmt_S2_tau_start"   , _pmt_S2_tau_start      , "pmt_S2_tau_start[5]/D"   );
	TBranch * _bn_pmt_S2_tau_end	= dpd->Branch("pmt_S2_tau_end"   , _pmt_S2_tau_end      , "pmt_S2_tau_end[5]/D"   );
	
	if (doCharge)
	{
		TBranch * _bn_tpc_track_charge	= dpd->Branch("tpc_track_charge",&_tpc_track_charge     , "tpc_track_charge[ntracks]/D"   );

		TBranch * _bn_tpc_track_mYX		= dpd->Branch("tpc_track_mYX"   	,_tpc_track_mYX      	, "tpc_track_mYX[ntracks]/D"   );
		TBranch * _bn_tpc_track_mZX		= dpd->Branch("tpc_track_mZX"   	,_tpc_track_mZX      	, "tpc_track_mZX[ntracks]/D"   );
		TBranch * _bn_tpc_track_nYX		= dpd->Branch("tpc_track_nYX"   	,_tpc_track_nYX      	, "tpc_track_nYX[ntracks]/D"   );
		TBranch * _bn_tpc_track_nZX		= dpd->Branch("tpc_track_nZX"   	,_tpc_track_nZX      	, "tpc_track_nZX[ntracks]/D"   );

		TBranch * _bn_tpc_track_startX	= dpd->Branch("tpc_track_startX"   	,_tpc_track_startX      	, "tpc_track_startX[ntracks]/D"   );
		TBranch * _bn_tpc_track_endX	= dpd->Branch("tpc_track_endX"   	,_tpc_track_endX      	, "tpc_track_endX[ntracks]/D"   );

		TBranch * _bn_tpc_track_start_theta		= dpd->Branch("tpc_track_start_theta"   	,_tpc_track_start_theta     	, "tpc_track_start_theta[ntracks]/D"   );
		TBranch * _bn_tpc_track_start_phi		= dpd->Branch("tpc_track_start_phi"   	,_tpc_track_start_phi      	, "tpc_track_start_phi[ntracks]/D"   );

		TBranch * _bn_tpc_track_fitresult_yx	= dpd->Branch("tpc_track_fitresult_yx"   	,_tpc_track_fitresult_yx      	, "tpc_track_fitresult_yx[ntracks]/I"   );
		TBranch * _bn_tpc_track_fitresult_zx	= dpd->Branch("tpc_track_fitresult_zx"   	,_tpc_track_fitresult_zx      	, "tpc_track_fitresult_zx[ntracks]/I"   );

		TBranch * _bn_tpc_track_length_n		= dpd->Branch("tpc_track_length"		, _tpc_track_length_n		, "tpc_track_length[ntracks]/D"       );
		TBranch * _bn_tpc_drift_time_at_pmt_pos	= dpd->Branch("tpc_drift_time_at_pmt_pos"		, _tpc_drift_time_at_pmt_pos		, "tpc_drift_time_at_pmt_pos[5]/D"       );
	}

	if (debug) cout << "Number of events: \t" << t2->GetEntries()  << endl;

	for(int ev=0; ev < t2->GetEntries() ; ++ev)
	{	
		cout << "dpdMaker: Event = " << ev << endl;
		
		if(ev%100==0) cout << ev << " over " <<t2->GetEntries()<<" events!" << endl;

 		t2->GetEntry(ev);
		
		if (!doCharge) _runlight=runNum;
		
		if (debug)
		{
			if (doCharge) printf("Charge (run,subrun,event) = %d, %d, %d\n",_runcharge,_subruncharge,_event); 
			printf("Light (run,event) = %d, %d\n",_runlight,_nevent);
		}
		
		// calculate delta t
		double time_charge=0.0;
		double time_light=0.0;
		
		time_light=(double)_PCTimeTag[0]+(double)_PCTimeTag[2]*1e-3;
		
		if (doCharge)
		{
			time_charge=(double)_eventtime_seconds-35.0+(double)_eventtime_nanoseconds*1e-9;
			
			if (debug)
			{
				cout << "\tTime: C_" << _eventtime_seconds << " L_" << _PCTimeTag[0] << endl;
				cout << "\tTime C_ns_" << _eventtime_nanoseconds << " L_ms_" << _PCTimeTag[2] << endl;
				cout << "\tdt: " <<time_charge-time_light << endl<< endl;
			}
			
			//
			// Charge stuff
			//
	
			TH2F * view[2];
			view[0] = new TH2F("view0","View 0;Channel;Time",320,0,320,1667,0,1667);
			view[1] = new TH2F("view1","View 1;Channel;Time",960,0,960,1667,0,1667);

			double pmt_pos[N_PMT]={0};

			pmt_pos[0]=150.0-92.246; // pmt position in Y coordinate (larsoft coordinate systems) in cm
			pmt_pos[1]=150.0-46.123;
			pmt_pos[2]=150.0-0;
			pmt_pos[3]=150.0+46.123;
			pmt_pos[4]=150.0+92.246;


			double S2_width_tolerance_channel=5;//cm

			double tpc_drift_time_at_pmt_pos[N_PMT]={0};
			int tpc_drift_time_at_pmt_pos_counter[N_PMT]={0};
			//if (debug) cout << _no_hits << endl;

			double totq=0;
			double totrecoq=0;
			double trackcharge[1000]={0};
			int number_reco_hits=0;
	
			TGraph *trackYX[KMAXNTRACKS], *trackZX[KMAXNTRACKS];
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) trackYX[j]= new TGraph();
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) trackZX[j]= new TGraph();
	
			//Dor (int j=0; j<_no_hits; j++)	if (debug) cout << "Ev "<<ev << " hit: " << j <<"/"<<_no_hits<< "  _charge:" <<  _Hit_ChargeIntegral[j] <<  endl;

			if (debug) cout << "Charge Event "<<_event <<" - n.hits "<<_no_hits<< " - n.tracks " << _NumberOfTracks_pmtrack << endl;
			for (int j=0; j<_no_hits; j++)
			{
		
				//if (debug) cout << "\t hit "<< j << " over " << _no_hits << endl;
				if(_Hit_ChargeIntegral[j]==_Hit_ChargeIntegral[j])
				{
					//if (debug) cout << "entro "  << endl;
					totq+=_Hit_ChargeIntegral[j];
					if(_hit_trkid[j]!=-1)
					{
						//if (debug) cout << "Hit belongs to the track " << j << " "  << _hit_trkid[j]   << endl;
						totrecoq+=_Hit_ChargeIntegral[j];
						number_reco_hits++;
						trackcharge[_hit_trkid[j]]+=_Hit_ChargeIntegral[j];
					}
					//if (debug) cout << "\t\t hitview "<< _Hit_View[j]  << endl;
					if( _Hit_View[j]==0) view[_Hit_View[j]]->Fill(_Hit_Channel[j],0.4*(1667-_Hit_PeakTime[j]),_Hit_ChargeIntegral[j]);
					else view[_Hit_View[j]]->Fill(_Hit_Channel[j]-320,0.4*(1667-_Hit_PeakTime[j]),_Hit_ChargeIntegral[j]);
				
				}
			}


			if (debug) cout << "Charge Event "<<_event << " Reco proportion (" <<Form("%.0d%s",100*number_reco_hits/_no_hits,"%") <<" hits) - "<< Form("%.2f%s",100*totrecoq/totq,"%") << "IADC"  <<  endl;
			//if (debug) cout << "Track " << 0 << " hits in view 0 "<< _Track_NumberOfHitsPerView_pmtrack[0][0] <<  endl;


			if (debug) for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "Track " << j << " hits in view 0 "<< _Track_NumberOfHitsPerView_pmtrack[j][0] <<  endl;
			if (debug) for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "Track " << j << " hits in view 1 "<< _Track_NumberOfHitsPerView_pmtrack[j][1] <<  endl;
			if (debug) for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "Total Charge per track " << j << ": "<< trackcharge[j] <<  endl;

	
			int counter=0;
			for (int j=0; j<_NumberOfTracks_pmtrack; j++)
			{
				for (int k=0; k<_Track_NumberOfHitsPerView_pmtrack[j][0]; k++)
				{
					trackYX[j]->SetPoint(trackYX[j]->GetN(),_Track_HitY_pmtrack[counter],_Track_HitX_pmtrack[counter]);
					trackZX[j]->SetPoint(trackZX[j]->GetN(),_Track_HitZ_pmtrack[counter],_Track_HitX_pmtrack[counter]);
					counter++;
				}
				for (int l=0; l<_Track_NumberOfHitsPerView_pmtrack[j][1]; l++)
				{
					trackYX[j]->SetPoint(trackYX[j]->GetN(),_Track_HitY_pmtrack[counter],_Track_HitX_pmtrack[counter]);
					trackZX[j]->SetPoint(trackZX[j]->GetN(),_Track_HitZ_pmtrack[counter],_Track_HitX_pmtrack[counter]);
					counter++;
				}
			}

			//cout << "\t... sorting 3d tracks by their charge." << endl;

			short sorted_tracks[KMAXNTRACKS]; // tracks sorted by total charge in the track

			std::size_t n(0);
			std::generate(std::begin(sorted_tracks), std::end(sorted_tracks), [&]{ return n++; });

			std::sort(  std::begin(sorted_tracks), std::end(sorted_tracks), [&](int i1, int i2) { return trackcharge[i1] > trackcharge[i2]; } );

			//Dor (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "\t\tSorted track " << j << " " << sorted_tracks[j] << " " << trackcharge[sorted_tracks[j]] <<"IAC"<< endl;
		
			//cout << "\t... 3d tracks sorted." << endl;
		
			// count number of tracks passing some basic selection criteria
			short NumberOfSelectedTracks=0;
			for (int j=0; j<_NumberOfTracks_pmtrack; j++)
			{
				// track is greater than 100 cm and reconstructed end-points 
				// must be closer than 2 cm to anode and cathode
					if (_Track_Length_pmtrack[j]>100. && 
						((fabs(_Track_StartX_pmtrack[j]-50.)<2. && fabs(_Track_EndX_pmtrack[j]+50.)<2.) || 
						 (fabs(_Track_StartX_pmtrack[j]+50.)<2. && fabs(_Track_EndX_pmtrack[j]-50.)<2.)))
								NumberOfSelectedTracks++;
			} 

			double t_yx_m[KMAXNTRACKS], t_zx_m[KMAXNTRACKS], t_yx_n[KMAXNTRACKS], t_zx_n[KMAXNTRACKS];
			int fail_yx[KMAXNTRACKS], fail_zx[KMAXNTRACKS];
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) trackYX[j]->LeastSquareLinearFit(trackYX[j]->GetN()-1,t_yx_n[j],t_yx_m[j], fail_yx[j],-50,50);
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) trackZX[j]->LeastSquareLinearFit(trackZX[j]->GetN()-1,t_zx_n[j],t_zx_m[j], fail_zx[j],0,300);
			if (debug) for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "TrackYX_m,n: "<< t_yx_m[j] << "  " << t_yx_n[j]<<endl;
			if (debug) for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "TrackZX_m,n: "<< t_zx_m[j] << "  " << t_zx_n[j]<<endl;
			if (debug) for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "Fit Results: "<< fail_yx[j] << "  " << fail_zx[j]<<endl;
			if (debug) for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "Number of points in track " << j << ": "<< trackYX[j]->GetN() << ", " << trackZX[j]->GetN() << " totq = " <<trackcharge[j]<<endl;
	
		 //LeastSquareLinearFit(Int_t n, Double_t& constant, Double_t& slope, Int_t& ifail, Double_t xmin = 0, Double_t xmax = 0)
			if (display_tracks) {
				for (int j=0; j<_NumberOfTracks_pmtrack; j++)
				{
					cout << "Track " << j << endl;
					TCanvas* c5 = new TCanvas("c5","c5",1200,600);
					c5->Divide(2,1);
					c5->cd(1);
					trackYX[j]->Draw("AP");
					TLine lin1 = TLine(-50.,t_yx_n[j]-50*t_yx_m[j],50.,t_yx_n[j]+50.*t_yx_m[j]);
					lin1.SetLineColor(kRed);
					lin1.Draw("same");
			
					c5->cd(2);
					trackZX[j]->Draw("AP");
					TLine lin2 = TLine(0.,t_zx_n[j],300.,t_zx_n[j]+300.*t_zx_m[j]);
					lin2.SetLineColor(kRed);
					lin2.Draw("same");
			
					c5->Modified();
					c5->Update();
					lets_pause();
				
					delete c5;
				}
			}

			for (int j=0; j<_no_hits; j++)
			{
				if(_Hit_ChargeIntegral[j]==_Hit_ChargeIntegral[j])
				{
					for (int k=0; k<N_PMT; k++) if (_hit_trkid[j]==sorted_tracks[0]&&_Hit_View[j]==1 && (_Hit_Channel[j]-320)*300.0/960.0>pmt_pos[k]-S2_width_tolerance_channel && (_Hit_Channel[j]-320)*300.0/960.0<pmt_pos[k]+S2_width_tolerance_channel ) 
					{
						tpc_drift_time_at_pmt_pos[k]+=0.4*_Hit_PeakTime[j];
						tpc_drift_time_at_pmt_pos_counter[k]++;
						//cout << "channel " <<  _Hit_Channel[j]<< " in cm " << _Hit_Channel[j]*300.0/960.0<< " time " << 0.4*_Hit_PeakTime[j] << " assing to pmt " << k << endl;
					}
				}
			}

			for (int k=0; k<N_PMT; k++)
			{
				if (tpc_drift_time_at_pmt_pos_counter[k]!=0) tpc_drift_time_at_pmt_pos[k]/=tpc_drift_time_at_pmt_pos_counter[k];
				else tpc_drift_time_at_pmt_pos[k]=-1;
			}
		
			if (debug)
			{
				cout << "Time: charge: " << time_charge << ", light: " << time_light << ", delta: " << time_charge - time_light << endl;
				cout << "Charge (Run,SubRun,ev): (" <<_runcharge << ","<< _subruncharge<<","<<_event<<") /t totcharge: "<<  totq <<endl; 
				cout << "   Length: (";
				for (Int_t i=0;i<_NumberOfTracks_pmtrack;i++)  cout << _Track_Length_pmtrack[i] << ", ";
				cout << ")"<< endl;
				cout << "   Charge: (";
				for (Int_t i=0;i<_NumberOfTracks_pmtrack;i++)cout << trackcharge[i] << ", ";
				cout << ")"<< endl;
			}
		
			// prepare all charge DPD variables for filling later
		
			_run_charge_out=_runcharge;
			_subrun_charge_out=_subruncharge;
			_nev_charge_out=_event;
		
			_ntracks=_NumberOfTracks_pmtrack;
			_n_sel_tracks=NumberOfSelectedTracks;

			_time_charge = time_charge;
		
			_tpc_totcharge=totq;
			_tpc_totrecocharge=totrecoq;

			for (int j=0; j<_NumberOfTracks_pmtrack; j++) 
			{
				_tpc_track_charge[j] = trackcharge[sorted_tracks[j]];
				_tpc_track_fitresult_yx[j]=fail_yx[sorted_tracks[j]];
				_tpc_track_fitresult_zx[j]=fail_zx[sorted_tracks[j]];
				_tpc_track_mYX[j] = t_yx_m[sorted_tracks[j]];
				_tpc_track_mZX[j] = t_zx_m[sorted_tracks[j]];
				_tpc_track_nYX[j] = t_yx_n[sorted_tracks[j]];
				_tpc_track_nZX[j] = t_zx_n[sorted_tracks[j]];
				_tpc_track_startX[j] = _Track_StartX_pmtrack[sorted_tracks[j]];
				_tpc_track_endX[j] = _Track_EndX_pmtrack[sorted_tracks[j]];
				_tpc_track_start_theta[j] = _Track_StartDirection_Theta_pmtrack[sorted_tracks[j]];
				_tpc_track_start_phi[j] = _Track_StartDirection_Phi_pmtrack[sorted_tracks[j]];
				_tpc_track_length_n[j] = _Track_Length_pmtrack[sorted_tracks[j]];
			}
		
			for (int k=0; k<N_PMT; k++) _tpc_drift_time_at_pmt_pos[k]=tpc_drift_time_at_pmt_pos[k];
		
			if (view[0]) view[0]->Delete();
			if (view[1]) view[1]->Delete();
		
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) if (trackYX[j]) trackYX[j]->Delete();
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) if (trackZX[j]) trackZX[j]->Delete();
		}
		
		//
		// Light stuff
		//
		
 		TH1F *h[N_PMT];
		TH1F *h_plot[N_PMT];
		TH1F *h_corr[N_PMT];
		TH1F *h_corr_p[N_PMT];
		TH1F *h_corr_m[N_PMT];
		double totlight[N_PMT+1]={0};
		int pmt_wvf_ADC_sat[N_PMT]={0};
		
		if (debug) cout << "\tLight variables: Nsamples " << _nsamples  << "\t" << _time_sample<< endl;

		for (int k=0; k<N_PMT; k++) 
		{
			h[k] = new TH1F(Form("%s%i","h",k),Form("%s%i%s","Channel ",k," ;Time [#mus]; Voltage [ADC counts]"),_nsamples,0.0,0.001*_nsamples*_time_sample);

			for (int j=0; j<_nsamples; j++) 
			{
				h[k]->SetBinContent(j+1, _adc_value[k][j]); 
				if (_adc_value[k][j] <=0 || _adc_value[k][j]>=4095) pmt_wvf_ADC_sat[k]=1;
				totlight[k]+= _adc_value[k][j];
			}
			h[k]->Sumw2();
		}
		if (debug) cout << "\tpmt histograms created " << ev << endl;
		
		double ped[N_PMT]={0};
		double pedrms[N_PMT]={0};
		double ped2[N_PMT]={0};
		double pedrms2[N_PMT]={0};
		double ped_end[N_PMT]={0};
		double ped_end_corr[N_PMT]={0};
		double wvf_end_fit_c0[N_PMT]={0};
		double wvf_end_fit_c1[N_PMT]={0};
		double wvf_end_fit_c2[N_PMT]={0};
		double wvf_end_fit_chi2[N_PMT]={0};
		int    wvf_end_fit_ndof[N_PMT]={0};
		double pedrms_end[N_PMT]={0};
		double pedrms_end_corr[N_PMT]={0};
		
		double ped_start = _nsamples==1000?0.0:218.;
		double ped_stop  = _nsamples==1000?0.5:228.;
		double ped2_start = _nsamples==1000?0.0:228.5;
		double ped2_stop  = _nsamples==1000?0.5:229.;
		double ped_end_start = _nsamples==1000?3.5:910.;
		double ped_end_stop  = _nsamples==1000?4.0:920.;
		
		int rebinfactor =_nsamples==1000?4:64;
		double S1_mintime = _nsamples==1000?0.3:228;
		double S1_maxtime = _nsamples==1000?0.65:231;
		double S2_maxtime = _nsamples==1000?4.0:900.;	
		
		vector<int> pmt_valleys[N_PMT];
		double pmt_valleys_tau[N_PMT][NMAXPEAKS]={0};
		
		for (int k=0; k<N_PMT; k++)	
		{	
			if (debug) cout << "\t... calculating pedestals"<< endl;
			
			ped[k] = WaveformAnalysis::baseline(h[k],pedrms[k],h[k]->FindBin(ped_start),h[k]->FindBin(ped_stop)-1); 
			
			ped2[k] = WaveformAnalysis::baseline(h[k],pedrms2[k],h[k]->FindBin(ped2_start),h[k]->FindBin(ped2_stop)-1); 
			
			ped_end[k] = WaveformAnalysis::baseline(h[k],pedrms_end[k],h[k]->FindBin(ped_end_start),h[k]->FindBin(ped_end_stop)-1); 
			
			if (debug) cout << "\t... correcting wf histogram and calculating end pedestal" << endl;
			
			WaveformAnalysis::correct_wvf_histo(h[k],h_corr[k],ped[k]);	
			WaveformAnalysis::correct_wvf_histo(h[k],h_corr_p[k],ped[k],313.8,313.8);	
			WaveformAnalysis::correct_wvf_histo(h[k],h_corr_m[k],ped[k],256.5,258.6);	
			
			ped_end_corr[k] = WaveformAnalysis::baseline(h_corr[k],pedrms_end_corr[k],h[k]->FindBin(ped_end_start),h[k]->FindBin(ped_end_stop)-1); 
			
			if (debug) cout << "\t... rebinning wf histogram for plotting and fitting"<< endl;
			
			h_plot[k] = dynamic_cast<TH1F*>(h[k]->Rebin(rebinfactor,Form("%s%i","h",k)));
			h_plot[k]->Scale(1.0/rebinfactor);
			
			if (debug) cout << "\t... finding peaks"<< endl;
			
			/*
			double thresh = ped[k]-15.*pedrms[k];
			pmt_valleys[k] = WaveformAnalysis::valleys(h_plot[k],thresh,1,h_plot[k]->FindBin(S1_maxtime));
			for (int i=0; i<pmt_valleys[k].size(); i++) pmt_valleys_tau[k][i]=hcenter(h_plot[k],pmt_valleys[k].at(i));
			*/
			
			pmt_valleys[k] = WaveformAnalysis::valleys(h_plot[k],1,h_plot[k]->FindBin(S2_maxtime));
			for (int i=0; i<TMath::Min(NMAXPEAKS,(int)pmt_valleys[k].size()); i++) pmt_valleys_tau[k][i]=hcenter(h_plot[k],pmt_valleys[k].at(i));
			
			if (_nsamples>1000 && (k==2 || k==3)) // only for positive-base PMTs
			{
				if (debug) cout << "\t... fitting end of wf"<< endl;
				
				/*
				TF1* fit = new TF1("fit","pol1",800,1100);
				h_plot[k]->Fit("fit","QRMN","",900,1020);
				wvf_end_fit_n[k] = fit->GetParameter(0);
				wvf_end_fit_m[k] = fit->GetParameter(1);
				wvf_end_fit_chi2[k] = fit->GetChisquare();
				wvf_end_fit_ndof[k] = fit->GetNDF();
				
				delete fit;
				*/
				
				//gStyle->SetOptStat(0);
				
				double ped_diff = ped_end[k]-ped[k];
				
				TF1* fit = new TF1("fit","[0]+[1]*exp((900-x)/[2])",850,1100);
				fit->SetParameter(0,ped[k]);
				fit->SetParameter(1,ped_diff);
				fit->SetParameter(2,250.);
				fit->SetParLimits(0,0.9*ped[k],1.1*ped[k]);
				fit->SetParLimits(1,0.9*ped_diff,1.1*ped_diff);
				fit->SetParLimits(2,10.,3000);
				h_plot[k]->Fit("fit","QRMN","",900,1040);
				wvf_end_fit_c0[k] = fit->GetParameter(0);
				wvf_end_fit_c1[k] = fit->GetParameter(1);
				wvf_end_fit_c2[k] = fit->GetParameter(2);
				wvf_end_fit_chi2[k] = fit->GetChisquare();
				wvf_end_fit_ndof[k] = fit->GetNDF();
				
				/*
				if (pedrms_end[k]<2.0 && !pmt_wvf_ADC_sat[k] && wvf_end_fit_chi2[k]<4.0 && wvf_end_fit_c2[k]>20. && wvf_end_fit_c2[k]<1000)
				{
					printf("ped = %f, ped_diff = %f\n",ped[k],ped_end[k]-ped[k]);
					printf("c0 = %f, c1 = %f, c2 = %f, chi2 = %f, pedend_rms = %f\n",wvf_end_fit_c0[k],wvf_end_fit_c1[k],wvf_end_fit_c2[k],wvf_end_fit_chi2[k],pedrms_end[k]);
				
					gStyle->SetOptStat(0);
					TCanvas *c1 = new TCanvas("c1");
					
					//h_plot[k]->GetXaxis()->SetRange(h_plot[k]->FindBin(650),h_plot[k]->FindBin(1048));
					h_plot[k]->Draw("hist");
					fit->SetLineColor(kRed);
					fit->Draw("same");
					c1->Modified();
					c1->Update();
					lets_pause();
				}
				*/
				
				
				delete fit;
			
				
				/*
				TFitResultPtr res = h_plot[k]->Fit("pol1","QRMSN","",900,1020);
				wvf_end_fit_n[k] = res->Parameter(0);
				wvf_end_fit_m[k] = res->Parameter(1);
				wvf_end_fit_chi2[k] = res->Chi2();
				wvf_end_fit_ndof[k] = res->Ndf();
				*/
			}
			
			if (debug) 
			{ 
				printf("ped[%i] = %0.2f +/- %0.2f\n", k, ped[k], pedrms[k]); 
				printf("ped_end[%i] = %0.2f +/- %0.2f\n", k, ped_end[k], pedrms_end[k]); 
				printf("ch: %d: No of valleys: %lu\n",k,pmt_valleys[k].size());
				for (int i=0; i<TMath::Min(NMAXPEAKS,(int)pmt_valleys[k].size()); i++) printf("\t%d: %d, %f\n",i,pmt_valleys[k].at(i),pmt_valleys_tau[k][i]);
				
				/*
				gStyle->SetOptStat(0);
				TCanvas *c1 = new TCanvas("c1");
				h_plot[k]->Draw("hist");
				c1->Modified();
				c1->Update();
				lets_pause();
				*/
				

			}
			if (0 && pmt_valleys[k].size()!=1) 
			{
				cout << "WARNING: " << pmt_valleys[k].size() << " valleys found! " << endl;
				debug=true;
				display_waveforms=true;
			}

			totlight[k]=-totlight[k]+ped[k]*_nsamples;
		}

		for (int k=0; k<N_PMT; k++)	totlight[N_PMT]+=totlight[k];

		//cout << "\t... integrating S1 for all PMTs."<< endl;
		
		double q=ADC_to_volts*_time_sample*1.e-9/(50.); // from iADC to Coulombs
		
		int    binpeak_S1[N_PMT]={0};
		int    endbin_S1[N_PMT]={0};
		double q_S1_1us[N_PMT]={0.0};
		double q_S1_4us[N_PMT]={0.0};
		double q_S1_80ns[N_PMT]={0.0};
		double q_S1_m2[N_PMT]={0.0};
		double q_S1_4us_corr[N_PMT]={0.0};
		double amp_S1[N_PMT]={0.0};
		double width_S1[N_PMT]={0.0};
		double tau_S1[N_PMT]={0};
		double tau_S1_end[N_PMT]={0}; 
		
		for (int k=0; k<N_PMT; k++)
		{
			int S1_minbin=h[k]->FindBin(S1_mintime);
			int S1_maxbin=h[k]->FindBin(S1_maxtime);
			
			//q_S1[k]=q*get_S1full(h[k],width_S1[k], binpeak_S1[k], ped[k]); // s1 in coulombs
			//tau_S1[k]=h[k]->GetBinCenter(binpeak_S1[k]); // in microseconds
			//amp_S1[k]=(ped[k]-h[k]->GetBinContent(binpeak_S1[k]))*ADC_to_volts; // in volts
			//printf("Jose: q = %f, amp = %f, width=%f, tau=%f\n",1e6*q_S1[k],amp_S1[k],width_S1[k],tau_S1[k]);
			
			if (debug) cout << "\t... calculating S1 parameters for all PMTs."<< endl;
			
			binpeak_S1[k] = WaveformAnalysis::find_S1_binpeak(h[k],S1_minbin,S1_maxbin);
			tau_S1[k] 	  = hcenter(h[k],binpeak_S1[k]);
			amp_S1[k] 	  = ADC_to_volts*(ped[k]-hget(h[k],binpeak_S1[k]));
			width_S1[k]   = _time_sample*WaveformAnalysis::calc_S1_width(h[k],binpeak_S1[k],ped[k]);
			q_S1_1us[k]   = q*WaveformAnalysis::calc_S1_charge(h[k],ped[k],tau_S1[k]-0.05,tau_S1[k]+1.0);
			q_S1_4us[k]   = q*WaveformAnalysis::calc_S1_charge(h[k],ped[k],tau_S1[k]-0.05,tau_S1[k]+4.0);
			q_S1_80ns[k]  = q*WaveformAnalysis::calc_S1_charge(h[k],ped[k],tau_S1[k]-0.05,tau_S1[k]+0.08);
			q_S1_m2[k]	  = q*WaveformAnalysis::calc_S1_charge_m2(h[k],ped[k],binpeak_S1[k],endbin_S1[k]);
			q_S1_4us_corr[k] = q*WaveformAnalysis::calc_S1_charge(h_corr[k],0.,tau_S1[k]-0.05,tau_S1[k]+4.0); 
			tau_S1_end[k] = hcenter(h[k],endbin_S1[k]);
	
			//printf("S1: q_1us = %f, q_4us = %f, q_m2 = %f, width = %f, amp = %f, tau = %f, tau_end = %f\n", 
			//1e6*q_S1_1us[k],1e6*q_S1_4us[k],1e6*q_S1_m2[k],width_S1[k],amp_S1[k],tau_S1[k],tau_S1_end[k]);
		}

		// S2 calculations
		
		int    binpeak_S2[N_PMT]={0};
		int    binavg_S2[N_PMT]={0};
		double q_S2[N_PMT]={0.0};
		double q_S2_corr[N_PMT]={0.0};
		double q_S2_corr_p[N_PMT]={0.0};
		double q_S2_corr_m[N_PMT]={0.0};
		double q_S2_m2[N_PMT]={0.0};
		double amp_S2[N_PMT]={0.0};
		double width_S2[N_PMT]={0};
		double tau_S2[N_PMT]={0}; 
		double tau_S2_start[N_PMT]={0}; 
		double tau_S2_end[N_PMT]={0}; 
		double tau_S2_avg[N_PMT]={0};
		
		// only do S2 calculation if sufficient nsamples
		if (_nsamples>1000)
		{
			/*
			if (debug) cout << "\t... calculating S2 parameters for all PMTs."<< endl;
		
			for (int k=0; k<N_PMT; k++)	
			{
				q_S2[k]= q*get_S2(h[k],width_S2[k],tau_S2[k],amp_S2[k],binpeak_S1[k],width_S1[k],ped[k]);
				printf("Jose: q = %f, amp = %f, tau=%f\n",1e6*q_S2[k],amp_S2[k],tau_S2[k]);
			}
			
			if (debug) cout << "\t... calculating S2_gaus parameteres for all PMTs."<< endl;
		
			for (int k=0; k<N_PMT; k++)	q_S2_gaus[k]= (ADC_to_volts*1e-6/50.)*get_S2_gaus(h_plot[k],width_S2_gaus[k], tau_S2_gaus[k], amp_S2_gaus[k],ped[k]);
			*/
			
			if (debug) cout << "\t... calculating S2 parameteres for all PMTs."<< endl;
			
			for (int k=0; k<N_PMT; k++)	
			{
				q_S2[k]			= q*WaveformAnalysis::calc_S2_parameters(h[k],ped[k],tau_S1[k],binpeak_S2[k],binavg_S2[k]);
				tau_S2[k] 		= hcenter(h[k],binpeak_S2[k]);
				tau_S2_avg[k] 	= hcenter(h[k],binavg_S2[k]);
				amp_S2[k] 		= ADC_to_volts*(ped[k]-hget(h[k],binpeak_S2[k]));
				q_S2_m2[k]		= ((double)rebinfactor)*q*WaveformAnalysis::calc_S2_parameters_m2(h_plot[k],tau_S1[k],tau_S2_start[k],tau_S2_end[k],width_S2[k]);

				int blah,blah2;
				q_S2_corr[k]    = q*WaveformAnalysis::calc_S2_parameters(h_corr[k],0.,tau_S1[k],blah,blah2);
				q_S2_corr_p[k]  = q*WaveformAnalysis::calc_S2_parameters(h_corr_p[k],0.,tau_S1[k],blah,blah2);
				q_S2_corr_m[k]  = q*WaveformAnalysis::calc_S2_parameters(h_corr_m[k],0.,tau_S1[k],blah,blah2);
			}
		}
		
		if (debug) 
		{
			cout << "\t... done."<< endl << endl;
			
			cout << "Light: (Run,ev): ("<<_runlight <<","<<_nevent <<") /ttotlight: " << totlight[N_PMT] << " " << _nsamples << endl;
			cout << "S1 light 1us: " << q_S1_1us[0]+q_S1_1us[1]+q_S1_1us[2]+q_S1_1us[3]+q_S1_1us[4] << " S2 light: " << q_S2[0]+q_S2[1]+q_S2[2]+q_S2[3]+q_S2[4]<< endl;
		
			cout <<"Pos1-STAIRS_CRT = "<< _crt_track_pos0[0] << " "<<_crt_track_pos0[1] << " "<<_crt_track_pos0[2]  << endl;
			cout <<"Pos2-WALL_CRT   = "<< _crt_track_pos1[0] << " "<<_crt_track_pos1[1] << " "<<_crt_track_pos1[2]  << endl;
		}

		double panel_space=100+20;
		double det_dy = 7302.;
		double y_center = (det_dy + panel_space*2)/2.;
		y_center=1500;
    	double z_center = (1910.+1930.)/2 + 112.*8;
		double fc_zmin = 2129;
		double fc_zmax = 2129+950;
		z_center = (fc_zmin+fc_zmax)/2 ;
		double fv_dz   = 950 + (618.08-418.08);
		double fv_zmin = 2129 -(618.08-418.08);
		double fv_zmax = fv_zmin + fv_dz;
		z_center = (fv_zmin+fv_zmax)/2;

		double x1,y1,z1,x2,y2,z2;
		double crtYX_m=0;
		double crtYX_n=0;

		double crtZX_m=0;
		double crtZX_n=0;
		TF1 *crtYX_line;
		TF1 *crtZX_line;

		if (_crt_daq_match==1 && _crt_reco==1)
		{
			_crt_matchreco=1;

			x2=_crt_track_pos0[2]*0.1-z_center*0.1;
			y2=-_crt_track_pos0[0]*0.1;
			z2=_crt_track_pos0[1]*0.1+y_center*0.1;//150;
			x1=_crt_track_pos1[2]*0.1-z_center*0.1;
			y1=-_crt_track_pos1[0]*0.1;
			z1=_crt_track_pos1[1]*0.1+y_center*0.1;//150;

			if (debug) cout <<"Pos1-STAIRS_LS = "<< x1 << " "<< y1 << " "<< z1  << endl;
			if (debug) cout <<"Pos2-WALL_LS   = "<< x2 << " "<< y2 << " "<< z2  << endl;

			crtYX_m=(x2-x1)/(y2-y1);
			crtYX_n=x2-crtYX_m*y2;

			crtZX_m=(x2-x1)/(z2-z1);
			crtZX_n=x2-crtZX_m*z2;

			if (debug) cout <<"YXm = "<< crtYX_m << " YXn = "<< crtYX_n << endl;
			if (debug) cout <<"ZXm = "<< crtZX_m << " ZXn = "<< crtZX_n << endl;

			crtYX_line = new TF1("YX_crt_track",Form("%.4f%s%.4f",crtYX_m,"*x+",crtYX_n),-50,50);
			crtZX_line = new TF1("ZX_crt_track",Form("%.4f%s%.4f",crtZX_m,"*x+",crtZX_n),0,3e2);

			crtZX_line->SetLineColor(kRed);
			crtYX_line->SetLineColor(kRed);
			crtZX_line->SetLineStyle(8);
			crtYX_line->SetLineStyle(8);
			
			if (crtYX_line) crtYX_line->Delete();
			if (crtZX_line) crtZX_line->Delete();
		}
		else _crt_matchreco=0;
		
		// prepare DPD variables for filling
		
		_run_light_out=_runlight;
		_nev_light_out=_nevent;
		_time_light = time_light;
		
		for (int k=0; k<N_PMT; k++) 
		{
			_pmt_wvf_properties[k][0]=ped[k];//ADC Counts
			_pmt_wvf_properties[k][1]=pedrms[k]; //ADC Counts
		
			_pmt_ped[k]=ped[k];
			_pmt_pedrms[k]=pedrms[k];
			_pmt_ped2[k]=ped2[k];
			_pmt_pedrms2[k]=pedrms2[k];
			_pmt_ped_end[k]=ped_end[k];
			_pmt_ped_end_corr[k]=ped_end_corr[k];
			_pmt_pedrms_end[k]=pedrms_end[k];
			_pmt_pedrms_end_corr[k]=pedrms_end_corr[k];
			_pmt_wvf_ADC_sat[k]=pmt_wvf_ADC_sat[k];
			_pmt_wvf_end_fit_c0[k]=wvf_end_fit_c0[k];
			_pmt_wvf_end_fit_c1[k]=wvf_end_fit_c1[k];
			_pmt_wvf_end_fit_c2[k]=wvf_end_fit_c2[k];
			_pmt_wvf_end_fit_chi2[k]=wvf_end_fit_chi2[k];
			_pmt_wvf_end_fit_ndof[k]=wvf_end_fit_ndof[k];
			
			_pmt_npeaks[k] = TMath::Min(NMAXPEAKS,(int)pmt_valleys[k].size());
			//for (int i=0; i<TMath::Min(_pmt_npeaks[k],NMAXPEAKS); i++) _pmt_peaks_tau[k][i]=pmt_valleys_tau[k][i]; 
			for (int i=0; i<NMAXPEAKS; i++) _pmt_peaks_tau[k][i]=pmt_valleys_tau[k][i];
	
			_pmt_charge[k]=(double)q*totlight[k];
			_pmt_npe[k]=(double)q*totlight[k]/(charge_e)/gains[k];

			_pmt_S1_charge_1us[k]=q_S1_1us[k]; // in volts
			_pmt_S1_charge_4us[k]=q_S1_4us[k]; // in volts
			_pmt_S1_charge_80ns[k]=q_S1_80ns[k]; // in volts
			_pmt_S1_charge_m2[k]=q_S1_m2[k]; // in volts
			_pmt_S1_charge_4us_corr[k]=q_S1_4us_corr[k]; // in volts
			
			_pmt_S1_npe_1us[k]=q_S1_1us[k]/(charge_e)/gains[k];
			_pmt_S1_npe_4us[k]=q_S1_4us[k]/(charge_e)/gains[k];
			_pmt_S1_npe_80ns[k]=q_S1_80ns[k]/(charge_e)/gains[k];
			_pmt_S1_npe_m2[k]=q_S1_m2[k]/(charge_e)/gains[k];
			_pmt_S1_npe_4us_corr[k]=q_S1_4us_corr[k]/(charge_e)/gains[k];
			
			_pmt_S1_width[k]=width_S1[k]; // in ns
			_pmt_S1_amp[k]=amp_S1[k]; // in volts
			_pmt_S1_tau[k]=tau_S1[k]; // in us
			_pmt_S1_tau_end[k]=tau_S1_end[k]; // in us

			_pmt_S2_charge[k]=q_S2[k];
			_pmt_S2_charge_m2[k]=q_S2_m2[k];
			_pmt_S2_charge_corr[k]=q_S2_corr[k];
			_pmt_S2_charge_corr_p[k]=q_S2_corr_p[k];
			_pmt_S2_charge_corr_m[k]=q_S2_corr_m[k];
			
			_pmt_S2_npe[k]=q_S2[k]/(charge_e)/gains[k];
			_pmt_S2_npe_m2[k]=q_S2_m2[k]/(charge_e)/gains[k];
			_pmt_S2_npe_corr[k]=q_S2_corr[k]/(charge_e)/gains[k];
			_pmt_S2_npe_corr_p[k]=q_S2_corr_p[k]/(charge_e)/gains[k];
			_pmt_S2_npe_corr_m[k]=q_S2_corr_m[k]/(charge_e)/gains[k];
			
			_pmt_S2_width[k]=width_S2[k]; // in us
			_pmt_S2_amp[k]=amp_S2[k]; // in volts
			_pmt_S2_tau[k]=tau_S2[k]; // in us
			_pmt_S2_tau_start[k]=tau_S2_start[k]; // in us
			_pmt_S2_tau_end[k]=tau_S2_end[k]; // in us
			_pmt_S2_tau_avg[k]=tau_S2_avg[k]; // in us
		}
		
		_pmt_charge[N_PMT]=0;
		_pmt_npe[N_PMT]=0;
		
		_pmt_S1_charge_1us[N_PMT]=0;
		_pmt_S1_npe_1us[N_PMT]=0;
		
		_pmt_S1_charge_4us[N_PMT]=0;
		_pmt_S1_npe_4us[N_PMT]=0;
		
		_pmt_S1_charge_80ns[N_PMT]=0;
		_pmt_S1_npe_80ns[N_PMT]=0;
		
		_pmt_S1_charge_m2[N_PMT]=0;
		_pmt_S1_npe_m2[N_PMT]=0;
		
		_pmt_S1_charge_4us_corr[N_PMT]=0;
		_pmt_S1_npe_4us_corr[N_PMT]=0;
		
		_pmt_S2_charge[N_PMT]=0;
		_pmt_S2_npe[N_PMT]=0;
		
		_pmt_S2_charge_m2[N_PMT]=0;
		_pmt_S2_npe_m2[N_PMT]=0;
		
		_pmt_S2_charge_corr[N_PMT]=0;
		_pmt_S2_npe_corr[N_PMT]=0;
		
		_pmt_S2_charge_corr_p[N_PMT]=0;
		_pmt_S2_npe_corr_p[N_PMT]=0;
		
		_pmt_S2_charge_corr_m[N_PMT]=0;
		_pmt_S2_npe_corr_m[N_PMT]=0;
		
		
		for (int k=0; k<N_PMT; k++) 
		{
			_pmt_charge[N_PMT]+=_pmt_charge[k];
			_pmt_npe[N_PMT]+=_pmt_npe[k];
		
			_pmt_S1_charge_1us[N_PMT]+=_pmt_S1_charge_1us[k]; // in volts
			_pmt_S1_npe_1us[N_PMT]+=_pmt_S1_npe_1us[k];
			
			_pmt_S1_charge_4us[N_PMT]+=_pmt_S1_charge_4us[k]; // in volts
			_pmt_S1_npe_4us[N_PMT]+=_pmt_S1_npe_4us[k];
			
			_pmt_S1_charge_80ns[N_PMT]+=_pmt_S1_charge_80ns[k]; // in volts
			_pmt_S1_npe_80ns[N_PMT]+=_pmt_S1_npe_80ns[k];
			
			_pmt_S1_charge_m2[N_PMT]+=_pmt_S1_charge_m2[k]; // in volts
			_pmt_S1_npe_m2[N_PMT]+=_pmt_S1_npe_m2[k];
			
			_pmt_S1_charge_4us_corr[N_PMT]+=_pmt_S1_charge_4us_corr[k]; // in volts
			_pmt_S1_npe_4us_corr[N_PMT]+=_pmt_S1_npe_4us_corr[k];

			_pmt_S2_charge[N_PMT]+=_pmt_S2_charge[k]; // in volts
			_pmt_S2_npe[N_PMT]+=_pmt_S2_npe[k];
			
			_pmt_S2_charge_m2[N_PMT]+=_pmt_S2_charge_m2[k]; // in volts
			_pmt_S2_npe_m2[N_PMT]+=_pmt_S2_npe_m2[k];
			
			_pmt_S2_charge_corr[N_PMT]+=_pmt_S2_charge_corr[k]; // in volts
			_pmt_S2_npe_corr[N_PMT]+=_pmt_S2_npe_corr[k];
			
			_pmt_S2_charge_corr_p[N_PMT]+=_pmt_S2_charge_corr_p[k]; // in volts
			_pmt_S2_npe_corr_p[N_PMT]+=_pmt_S2_npe_corr_p[k];
			
			_pmt_S2_charge_corr_m[N_PMT]+=_pmt_S2_charge_corr_m[k]; // in volts
			_pmt_S2_npe_corr_m[N_PMT]+=_pmt_S2_npe_corr_m[k];
		}

		_crt_mYX = crtYX_m;
		_crt_mZX = crtZX_m;
		_crt_nYX = crtYX_n;
		_crt_nZX = crtZX_n;
		_crt_ToF_n = _crt_ToF;
		_crt_track_lenFV_n = _crt_track_lenFV;
		for (int k=0; k<3; k++) _crt_pmt_dist_n[k] = _crt_pmt_dist[k];
		_crt_isFV_n = _crt_isFV;
		
		dpd->Fill();
		
		
		
		// Wf displays	
		
		TCanvas *c1;
		//if (_pmt_S2_tau_end[3]<500. && fabs(_pmt_ped[3]-_pmt_ped_end[3])>5) display_waveforms=true;
					
		if (display_waveforms)
		{
			gStyle->SetOptStat(0);
		    gStyle->SetPadTickX(1);
		    gStyle->SetPadTickY(1);
			
			c1 = new TCanvas("c1","c1",1200,800);
			c1->Divide(3,2);
			for (int k=0; k<N_PMT; k++)
			{
				c1->cd(k+1);
				
				h_plot[k]->GetYaxis()->SetTitleOffset(1.25);
				h_plot[k]->GetYaxis()->SetTitleSize(0.05);
				h_plot[k]->GetXaxis()->SetTitleSize(0.05);
				
				h[k]->Draw("HIST");
			
				double x = 0.15;
				double y = 0.39;
				TLatex l;
				l.SetNDC();
				//l.DrawLatex(x,y+0.06,Form("Run %d, Event %d",_run_out,ev));
				l.DrawLatex(x,y,Form("ped = %0.2f +/- %0.2f",ped[k],pedrms[k]));
				l.DrawLatex(x,y-0.06,Form("ped_end = %0.2f +/- %0.2f",ped_end[k],pedrms_end[k]));
				l.DrawLatex(x,y-0.12,Form("S1: %0.2f nC, %0.2f V, %0.1f ns, %0.1f #mus",1e9*_pmt_S1_charge_1us[k],_pmt_S1_amp[k],_pmt_S1_width[k],_pmt_S1_tau[k]));
				if (_pmt_S2_amp[k]>0.) 
				{
					l.DrawLatex(x,y-0.18,Form("S2: %0.2f uC, %0.2f mV, %0.1f #mus",1e6*_pmt_S2_charge[k],_pmt_S2_amp[k]*1000.,_pmt_S2_tau[k]));
					l.DrawLatex(x,y-0.24,Form("S2_m2: %0.2f #mus, %0.1f #mus, %0.1f #mus",_pmt_S2_tau_start[k],_pmt_S2_tau_end[k],_pmt_S2_width[k]));
				}	
				//l.DrawLatex(x,y-0.30,Form("ped_diff = %0.2f",(ped_end[k]-ped[k])/pedrms[k]));
			}
			c1->Update();
			c1->Modified();
			//c1->Print(Form("working_v2/S1_amp_1/run%d_%d_ev%d.png",_run_out,_subrun_out,_event));	
			
			lets_pause();
		}
		
		//if (debug || display_waveforms) lets_pause();
		
		//if (c1) delete c1;
		
		//debug=false;
		//display_waveforms=false;

		for (int k=0; k<N_PMT; k++) if (h[k]) h[k]->Delete(); 
		for (int k=0; k<N_PMT; k++) if (h_plot[k]) h_plot[k]->Delete(); 
		for (int k=0; k<N_PMT; k++) if (h_corr[k]) h_corr[k]->Delete(); 
		for (int k=0; k<N_PMT; k++) if (h_corr_p[k]) h_corr_p[k]->Delete(); 
		for (int k=0; k<N_PMT; k++) if (h_corr_m[k]) h_corr_m[k]->Delete(); 

	} //event loop

	dpd->Write();
	ofile.Close();
}

