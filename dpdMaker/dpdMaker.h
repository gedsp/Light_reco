#include "../LIB/WaveformAnalysis.cc"
#include "../LIB/ChargeCalibration.cc"
#include "../config_reco.h"
//#include "THDet.C"
#include "../LIB/HighwayReader.h"
//#include "../LIB/myFunc.h"

bool debug=false;
bool display_tracks=false;
bool display_waveforms=false;
bool display_crt=false;

double ADC_to_volts = 2./4096.;
double charge_e = 1.602E-19;
double chargeCalib[2]={59.8053E15,66.6715E15}; // ADC*ticks / Coulomb (for each view)

static const float Z_CENTER = 1500.0; // 1577.75;

// DPD variables

int _run_charge_out;
int _subrun_charge_out;
int _nevent_charge_out;
int _trig_conf;

int _run_light_out;
int _nevent_light_out;

double _time_charge;
double _time_light;

double _pmt_charge[N_PMT+1];
double _pmt_npe[N_PMT+1];
double _pmt_ped[N_PMT];
double _pmt_pedrms[N_PMT];
double _pmt_ped_uncorr[N_PMT];
double _pmt_pedrms_uncorr[N_PMT];
double _pmt_ped2[N_PMT];
double _pmt_pedrms2[N_PMT];
double _pmt_ped_end[N_PMT];
double _pmt_ped_end_uncorr[N_PMT];
double _pmt_pedrms_end[N_PMT];
double _pmt_pedrms_end_uncorr[N_PMT];
int	   _pmt_wvf_ADC_sat[N_PMT];
double _pmt_wvf_end_fit_c0[N_PMT];
double _pmt_wvf_end_fit_c1[N_PMT];
double _pmt_wvf_end_fit_c2[N_PMT];
double _pmt_wvf_end_fit_chi2[N_PMT];
int    _pmt_wvf_end_fit_ndof[N_PMT];
double _pmt_wvf_corr_RC_eff[N_PMT];
double _pmt_wvf_max_uncorr[N_PMT];
double _pmt_wvf_max_corr[N_PMT];
double _pmt_wvf_properties[N_PMT][2];
int    _pmt_npeaks[N_PMT];
double _pmt_peaks_tau[N_PMT][KMAXNPEAKS];

double _pmt_S1_charge_1us[N_PMT+1];
double _pmt_S1_charge_4us[N_PMT+1];
double _pmt_S1_charge_80ns[N_PMT+1];
double _pmt_S1_charge_m2[N_PMT+1];
double _pmt_S1_npe_1us[N_PMT+1];
double _pmt_S1_npe_4us[N_PMT+1];
double _pmt_S1_npe_80ns[N_PMT+1];
double _pmt_S1_npe_m2[N_PMT+1];
double _pmt_S1_width[N_PMT];
double _pmt_S1_amp[N_PMT];
double _pmt_S1_tau[N_PMT];
double _pmt_S1_tau_end[N_PMT];

double _pmt_S2_charge[N_PMT+1];
double _pmt_S2_charge_m2[N_PMT+1];
double _pmt_S2_charge_corr_p[N_PMT+1];
double _pmt_S2_charge_corr_m[N_PMT+1];
double _pmt_S2_npe[N_PMT+1];
double _pmt_S2_npe_m2[N_PMT+1];
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
double _tpc_totcharge_anode;
double _tpc_totrecocharge_anode;
double _tpc_totcharge_anode_view[2];
double _tpc_totrecocharge_anode_view[2];
double _tpc_totcharge_LAr_view[2];
double _tpc_totrecocharge_LAr_view[2];

double _tpc_track_charge[KMAXNTRACKS];
double _tpc_track_charge_anode[KMAXNTRACKS];
double _tpc_track_charge_anode_view[KMAXNTRACKS][2];
double _tpc_track_charge_LAr_view[KMAXNTRACKS][2];
double _tpc_track_dE_view[KMAXNTRACKS][2];
double _tpc_track_avg_dEds_view[KMAXNTRACKS][2];

int   _crt_matchreco;
float _crt_fit_yx_m;
float _crt_fit_yx_n;
float _crt_fit_zx_m;
float _crt_fit_zx_n;

short _ntracks;
short _n_sel_tracks;

float _tpc_track_start_theta[KMAXNTRACKS];
float _tpc_track_start_phi[KMAXNTRACKS];
float _tpc_track_end_theta[KMAXNTRACKS];
float _tpc_track_end_phi[KMAXNTRACKS];

float _tpc_track_start_x[KMAXNTRACKS];
float _tpc_track_start_y[KMAXNTRACKS];
float _tpc_track_start_z[KMAXNTRACKS];
float _tpc_track_end_x[KMAXNTRACKS];
float _tpc_track_end_y[KMAXNTRACKS];
float _tpc_track_end_z[KMAXNTRACKS];

float _tpc_track_fit_yx_m[KMAXNTRACKS];
float _tpc_track_fit_yx_n[KMAXNTRACKS];
float _tpc_track_fit_yx_chi2[KMAXNTRACKS];
int	  _tpc_track_fit_yx_ndof[KMAXNTRACKS];

float _tpc_track_fit_zx_m[KMAXNTRACKS];
float _tpc_track_fit_zx_n[KMAXNTRACKS];
float _tpc_track_fit_zx_chi2[KMAXNTRACKS];
int	  _tpc_track_fit_zx_ndof[KMAXNTRACKS];

short _tpc_track_id[KMAXNTRACKS];
short _tpc_track_nhits_view[KMAXNTRACKS][2];
float _tpc_track_momentum[KMAXNTRACKS];
float _tpc_track_length_trajectory[KMAXNTRACKS];
float _tpc_track_length_straight_line[KMAXNTRACKS];

float _tpc_drift_time_at_pmt[N_PMT];
float _tpc_drift_time_at_pmt2[N_PMT];

float _tpc_track_max_CBR_view[KMAXNTRACKS][2];


void initDPDvariables();

void make_dpd(TChain* t2, int runNum, int trigConf, double drift_field, double gains[N_PMT], string outfilename, int light_subrun);



void make_dpd(TChain* t2, int runNum, int trigConf, double drift_field, double gains[N_PMT], string outfilename, int light_subrun=-1)
{
	bool doCharge = (TString(outfilename).Contains("matched"))?true:false;
	
	cout << "doCharge = " << doCharge << ", runNum = " << runNum << endl;
	cout << "Drift field for this run: " << drift_field << " kV/cm" << endl;
	
	TFile ofile(outfilename.c_str(),"RECREATE");
	if (ofile.IsZombie())
	{
	       cout << "Error opening output file " << outfilename.c_str() << endl;
	       exit(-1);
	}

	TTree *dpd= new TTree("dpd","result tree");
	
	
	// Input charge variables
	int _runcharge;
	int _subruncharge;
	int _event;
    int _eventtime_seconds;
    int _eventtime_nanoseconds;
	
	int   _NumberOfHits;
	short _NumberOfTracks;
	//short _nclusters;
	
	short _Hit_TrackID[KMAXNHITS];
	float _Hit_ChargeIntegral[KMAXNHITS];
	float _Hit_ChargeSummedADC[KMAXNHITS];
	float _Hit_PeakTime[KMAXNHITS];
	short _Hit_Channel[KMAXNHITS];
	short _Hit_View[KMAXNHITS];
	short _Hit_TPC[KMAXNHITS];

	short _TrackID[KMAXNTRACKS];
	float _Track_TPC[KMAXNHITS];
	short _Track_NumberOfHits[KMAXNTRACKS];
	float _Track_PitchInViews[KMAXNTRACKS][2];
	short _Track_NumberOfHitsPerView[KMAXNTRACKS][2];
	float _Track_Momentum[KMAXNTRACKS];
	float _Track_Length_Trajectory[KMAXNTRACKS];
	float _Track_Length_StraightLine[KMAXNTRACKS];
	float _Track_StartDirection_Theta[KMAXNTRACKS];
	float _Track_StartDirection_Phi[KMAXNTRACKS];
	float _Track_EndDirection_Theta[KMAXNTRACKS];
	float _Track_EndDirection_Phi[KMAXNTRACKS];
	float _Track_StartX[KMAXNTRACKS];
	float _Track_StartY[KMAXNTRACKS];
	float _Track_StartZ[KMAXNTRACKS];
	float _Track_EndX[KMAXNTRACKS];
	float _Track_EndY[KMAXNTRACKS];
	float _Track_EndZ[KMAXNTRACKS];
	
	float _Track_HitX[KMAXNHITS];
	float _Track_HitY[KMAXNHITS];
	float _Track_HitZ[KMAXNHITS];
	short _Track_Hit_View[KMAXNHITS];
	float _Track_Hit_ChargeSummedADC[KMAXNHITS]; 
	float _Track_Hit_ChargeIntegral[KMAXNHITS];
	float _Track_Hit_PeakTime[KMAXNHITS];
	short _Track_Hit_Channel[KMAXNHITS];
	float _Track_Hit_ds_LocalTrackDirection[KMAXNHITS];
	float _Track_Hit_ds_3DPosition[KMAXNHITS];
	
	//float _deltatime;
	
	// Input light Variables
	int _PCTimeTag[3];
	int _nsamples;
	int _time_stamp;
    int _nevent;
    int _nchannels;
  	int _crt_daq_match;
  	int _crt_reco;
	int _crt_reco_sat;
    int _time_sample;
	int _runlight;
	short _adc_value[5][300000];
	//int	_crt_adc[4][32];

	int	  _crt_isFV;
	int   _crt_isFC;
	int   _crt_isclosestpoint[5];
	float _crt_pmt_dist[5];
	float _crt_track_param[5];
	float _crt_drift_len[5];
	float _crt_closest_coord[5][3];
	float _crt_track_door[3];
	float _crt_track_wall[3];
	float _crt_point_door_fc[3];
	float _crt_point_wall_fc[3];
	float _crt_point_door_fv[3];
	float _crt_point_wall_fv[3];
	float _crt_ToF; 

	

	if (doCharge)
	{
		 // Input charge branches
		
		TBranch * _b2_runcharge;	  			t2->SetBranchAddress("Run", 					&_runcharge, 						&_b2_runcharge);
		TBranch * _b2_subruncharge;	  			t2->SetBranchAddress("Subrun", 					&_subruncharge, 					&_b2_subruncharge);
		TBranch * _b2_event;	  				t2->SetBranchAddress("EventNumberInRun",		&_event , 							&_b2_event);
		TBranch * _b2_eventtime_seconds;		t2->SetBranchAddress("EventTimeSeconds", 		&_eventtime_seconds,				&_b2_eventtime_seconds);
		TBranch * _b2_eventtime_nanoseconds;	t2->SetBranchAddress("EventTimeNanoseconds",	&_eventtime_nanoseconds,			&_b2_eventtime_nanoseconds);

		TBranch * _b2_NumberOfHits;				t2->SetBranchAddress("NumberOfHits", 			&_NumberOfHits, 					&_b2_NumberOfHits);
		TBranch * _b2_NumberOfTracks;			t2->SetBranchAddress("NumberOfTracks", 			&_NumberOfTracks, 					&_b2_NumberOfTracks);
		//TBranch * _b2_nclusters;				t2->SetBranchAddress("NumberOfClusters", 		&_nclusters, 						&_b2_nclusters);
		
		TBranch * _b2_Hit_TrackID;				t2->SetBranchAddress("Hit_TrackID", 			_Hit_TrackID, 						&_b2_Hit_TrackID);
		TBranch * _b2_Hit_ChargeIntegral;		t2->SetBranchAddress("Hit_ChargeIntegral", 		_Hit_ChargeIntegral, 				&_b2_Hit_ChargeIntegral);
		TBranch * _b2_Hit_ChargeSummedADC;		t2->SetBranchAddress("Hit_ChargeSummedADC", 	_Hit_ChargeSummedADC, 				&_b2_Hit_ChargeSummedADC);
		TBranch * _b2_Hit_TPC;					t2->SetBranchAddress("Hit_TPC", 				_Hit_TPC, 							&_b2_Hit_TPC);
		TBranch * _b2_Hit_View;					t2->SetBranchAddress("Hit_View", 				_Hit_View, 							&_b2_Hit_View);
		TBranch * _b2_Hit_Channel;				t2->SetBranchAddress("Hit_Channel", 			_Hit_Channel, 						&_b2_Hit_Channel);
		TBranch * _b2_Hit_PeakTime;				t2->SetBranchAddress("Hit_PeakTime", 			_Hit_PeakTime, 						&_b2_Hit_PeakTime);

		TBranch * _b2_TrackID;					t2->SetBranchAddress("TrackID", 					_TrackID, 						&_b2_TrackID);
		TBranch * _b2_Track_TPC;				t2->SetBranchAddress("Track_Hit_TPC", 				_Track_TPC, 					&_b2_Track_TPC);
		TBranch * _b2_Track_Momentum;			t2->SetBranchAddress("Track_Momentum", 				_Track_Momentum, 				&_b2_Track_Momentum);
		TBranch * _b2_Track_NumberOfHits;		t2->SetBranchAddress("Track_NumberOfHits", 			_Track_NumberOfHits, 			&_b2_Track_NumberOfHits);
		TBranch * _b2_Track_NumberOfHitsPerView;t2->SetBranchAddress("Track_NumberOfHitsPerView",	_Track_NumberOfHitsPerView, 	&_b2_Track_NumberOfHitsPerView);
		TBranch * _b2_Track_Length_Trajectory;	t2->SetBranchAddress("Track_Length_Trajectory", 	_Track_Length_Trajectory, 		&_b2_Track_Length_Trajectory);
		TBranch * _b2_Track_Length_StraightLine;t2->SetBranchAddress("Track_Length_StraightLine",	_Track_Length_StraightLine, 	&_b2_Track_Length_StraightLine);
		TBranch * _b2_Track_PitchInViews;		t2->SetBranchAddress("Track_PitchInViews", 			_Track_PitchInViews, 			&_b2_Track_PitchInViews);
		TBranch * _b2_Track_StartDirection_Theta;t2->SetBranchAddress("Track_StartDirection_Theta", _Track_StartDirection_Theta,	&_b2_Track_StartDirection_Theta);
		TBranch * _b2_Track_StartDirection_Phi;	t2->SetBranchAddress("Track_StartDirection_Phi", 	_Track_StartDirection_Phi, 		&_b2_Track_StartDirection_Phi);
		TBranch * _b2_Track_EndDirection_Theta;	t2->SetBranchAddress("Track_EndDirection_Theta", 	_Track_EndDirection_Theta, 		&_b2_Track_EndDirection_Theta);
		TBranch * _b2_Track_EndDirection_Phi;	t2->SetBranchAddress("Track_EndDirection_Phi", 		_Track_EndDirection_Phi, 		&_b2_Track_EndDirection_Phi);
		TBranch * _b2_Track_StartX;				t2->SetBranchAddress("Track_StartPoint_X", 			_Track_StartX, 					&_b2_Track_StartX);
		TBranch * _b2_Track_StartY;				t2->SetBranchAddress("Track_StartPoint_Y", 			_Track_StartY, 					&_b2_Track_StartY);
		TBranch * _b2_Track_StartZ;				t2->SetBranchAddress("Track_StartPoint_Z", 			_Track_StartZ, 					&_b2_Track_StartZ);
		TBranch * _b2_Track_EndX;				t2->SetBranchAddress("Track_EndPoint_X", 			_Track_EndX, 					&_b2_Track_EndX);
		TBranch * _b2_Track_EndY;				t2->SetBranchAddress("Track_EndPoint_Y", 			_Track_EndY, 					&_b2_Track_EndY);
		TBranch * _b2_Track_EndZ;				t2->SetBranchAddress("Track_EndPoint_Z", 			_Track_EndZ, 					&_b2_Track_EndZ);
		
		TBranch * _b2_Track_Hit_View;			t2->SetBranchAddress("Track_Hit_View", 				_Track_Hit_View, 				&_b2_Track_Hit_View);
		TBranch * _b2_Track_Hit_PeakTime;		t2->SetBranchAddress("Track_Hit_PeakTime", 			_Track_Hit_PeakTime, 			&_b2_Track_Hit_PeakTime);
		TBranch * _b2_Track_Hit_ChargeSummedADC;t2->SetBranchAddress("Track_Hit_ChargeSummedADC", 	_Track_Hit_ChargeSummedADC, 	&_b2_Track_Hit_ChargeSummedADC);
		TBranch * _b2_Track_Hit_ChargeIntegral;	t2->SetBranchAddress("Track_Hit_ChargeIntegral", 	_Track_Hit_ChargeIntegral, 		&_b2_Track_Hit_ChargeIntegral);
		TBranch * _b2_Track_Hit_ds_LocalTrackDirection; t2->SetBranchAddress("Track_Hit_ds_LocalTrackDirection", _Track_Hit_ds_LocalTrackDirection, &_b2_Track_Hit_ChargeIntegral);
		TBranch * _b2_Track_Hit_ds_3DPosition;	t2->SetBranchAddress("Track_Hit_ds_3DPosition", 	_Track_Hit_ds_3DPosition, 		&_b2_Track_Hit_ds_3DPosition);
		TBranch * _b2_Track_Hit_Channel; 		t2->SetBranchAddress("Track_Hit_Channel", 			_Track_Hit_Channel, 			&_b2_Track_Hit_Channel);
		TBranch * _b2_Track_HitX;				t2->SetBranchAddress("Track_Hit_X", 				_Track_HitX, 					&_b2_Track_HitX);
		TBranch * _b2_Track_HitY;				t2->SetBranchAddress("Track_Hit_Y", 				_Track_HitY, 					&_b2_Track_HitY);
		TBranch * _b2_Track_HitZ;				t2->SetBranchAddress("Track_Hit_Z", 				_Track_HitZ, 					&_b2_Track_HitZ);

		// Light branches whose names were changed during matching
		TBranch * _b2_runlight;					t2->SetBranchAddress("light_run", 					&_runlight, 					&_b2_runlight);
		TBranch * _b2_light_event;				t2->SetBranchAddress("light_event", 				&_nevent, 						&_b2_light_event);
		
		// Light branches that were created during matching
		//TBranch * _b2_deltatime;   				t2->SetBranchAddress("deltatime",     		&_deltatime, 			&_b2_deltatime);
	}
	
	else 
	{ 
		// Use names from light branches in re-processed files 
		TBranch * _b2_light_event;	t2->SetBranchAddress("event"		, &_nevent    		, &_b2_light_event	);
	}
		
		
	// common light + CRT branches
	
	TBranch * _b2_time_sample;			t2->SetBranchAddress("TimeSample", 			&_time_sample, 			&_b2_time_sample);	
	TBranch * _b2_nchannels;			t2->SetBranchAddress("nchannels", 			&_nchannels, 			&_b2_nchannels);
	TBranch * _b2_PCTimeTag;			t2->SetBranchAddress("PCTimeTag", 			_PCTimeTag, 			&_b2_PCTimeTag);
	TBranch * _b2_nsamples;				t2->SetBranchAddress("nsamples", 			&_nsamples,  			&_b2_nsamples);	

	TBranch * _b2_adc_value_0; 			t2->SetBranchAddress("adc_value_0", 		_adc_value[0], 			&_b2_adc_value_0);
	TBranch * _b2_adc_value_1; 			t2->SetBranchAddress("adc_value_1", 		_adc_value[1], 			&_b2_adc_value_1);
	TBranch * _b2_adc_value_2; 			t2->SetBranchAddress("adc_value_2", 		_adc_value[2], 			&_b2_adc_value_2);
	TBranch * _b2_adc_value_3;  		t2->SetBranchAddress("adc_value_3", 		_adc_value[3], 			&_b2_adc_value_3);
	TBranch * _b2_adc_value_4;  		t2->SetBranchAddress("adc_value_4", 		_adc_value[4], 			&_b2_adc_value_4);

	TBranch * _b2_crt_daq_match;		t2->SetBranchAddress("crt_daq_match", 		&_crt_daq_match, 		&_b2_crt_daq_match);
	TBranch * _b2_crt_reco;				t2->SetBranchAddress("crt_reco", 			&_crt_reco, 			&_b2_crt_reco);
	TBranch * _b2_crt_reco_sat;			t2->SetBranchAddress("crt_reco_sat", 		&_crt_reco_sat, 		&_b2_crt_reco_sat);
	TBranch * _b2_crt_isclosestpoint;	t2->SetBranchAddress("crt_isclosestpoint", 	_crt_isclosestpoint, 	&_b2_crt_isclosestpoint);
	TBranch * _b2_crt_ToF;  			t2->SetBranchAddress("crt_ToF", 			&_crt_ToF, 				&_b2_crt_ToF);
	TBranch * _b2_crt_isFC; 	 		t2->SetBranchAddress("crt_isFC", 			&_crt_isFC, 			&_b2_crt_isFC);
	TBranch * _b2_crt_isFV;  			t2->SetBranchAddress("crt_isFV", 			&_crt_isFV, 			&_b2_crt_isFV);
	
	TBranch * _b2_crt_track_param;  t2->SetBranchAddress("crt_track_param", 		_crt_track_param, 		&_b2_crt_track_param);
	TBranch * _b2_crt_drift_len;  	t2->SetBranchAddress("crt_drift_len", 			_crt_drift_len, 		&_b2_crt_drift_len);
	TBranch * _b2_crt_pmt_dist;  	t2->SetBranchAddress("crt_pmt_dist", 			_crt_pmt_dist, 			&_b2_crt_pmt_dist);
	TBranch * _b2_crt_closest_coord;t2->SetBranchAddress("crt_closest_coord", 		_crt_closest_coord, 	&_b2_crt_closest_coord);
	TBranch * _b2_crt_track_door;  	t2->SetBranchAddress("crt_track_door", 			_crt_track_door, 		&_b2_crt_track_door);
	TBranch * _b2_crt_track_wall;  	t2->SetBranchAddress("crt_track_wall", 			_crt_track_wall, 		&_b2_crt_track_wall);
	TBranch * _b2_crt_point_door_fc;t2->SetBranchAddress("crt_point_door_fc", 		_crt_point_door_fc, 	&_b2_crt_point_door_fc);
	TBranch * _b2_crt_point_wall_fc;t2->SetBranchAddress("crt_point_wall_fc", 		_crt_point_wall_fc, 	&_b2_crt_point_wall_fc);
	TBranch * _b2_crt_point_door_fv;t2->SetBranchAddress("crt_point_door_fv", 		_crt_point_door_fv, 	&_b2_crt_point_door_fv);
	TBranch * _b2_crt_point_wall_fv;t2->SetBranchAddress("crt_point_wall_fv", 		_crt_point_wall_fv, 	&_b2_crt_point_wall_fv);
	
	

	
	//DPD branches
	
	if (doCharge)
	{
		TBranch * _bn_run_charge				= dpd->Branch("run_charge", 			&_run_charge_out, 		"run_charge/I");
		TBranch * _bn_event_charge					= dpd->Branch("event_charge", 				&_nevent_charge_out, 		"event_charge/I");
		TBranch * _bn_subrun_charge 			= dpd->Branch("subrun_charge", 			&_subrun_charge_out, 	"subrun_charge/I");
		TBranch * _bn_time_charge 				= dpd->Branch("time_charge", 			&_time_charge  , 		"time_charge/D");	
		//TBranch * _bn_deltatime     			= dpd->Branch("deltatime", 				&_deltatime, "deltatime/F");
	}
	
	TBranch * _bn_run_light						= dpd->Branch("run_light", 				&_run_light_out, 		"run_light/I");
	TBranch * _bn_trig_conf						= dpd->Branch("trig_conf", 				&_trig_conf, 	   		"trig_conf/I");
	TBranch * _bn_event_light						= dpd->Branch("event_light", 				&_nevent_light_out, 		"event_light/I");
	TBranch * _bn_nsamples 						= dpd->Branch("nsamples", 				&_nsamples, 			"nsamples/I");
	TBranch * _bn_time_light 					= dpd->Branch("time_light", 			&_time_light, 			"time_light/D");
	
	if (doCharge)
	{
		TBranch * _bn_ntracks						= dpd->Branch("ntracks",						&_ntracks, 						"ntracks/S");
		TBranch * _bn_n_sel_tracks					= dpd->Branch("n_sel_tracks",					&_n_sel_tracks, 				"n_sel_tracks/S");
		TBranch * _bn_tpc_totcharge 				= dpd->Branch("tpc_totcharge", 					&_tpc_totcharge, 				"tpc_totcharge/D");	
		TBranch * _bn_tpc_totrecocharge 			= dpd->Branch("tpc_totrecocharge", 				&_tpc_totrecocharge, 			"tpc_totrecocharge/D");	
		TBranch * _bn_tpc_totcharge_anode			= dpd->Branch("tpc_totcharge_anode", 			&_tpc_totcharge_anode, 			"tpc_totcharge_anode/D");	
		TBranch * _bn_tpc_totrecocharge_anode 		= dpd->Branch("tpc_totrecocharge_anode", 		&_tpc_totrecocharge_anode, 		"tpc_totrecocharge_anode/D");	
		TBranch * _bn_tpc_totcharge_anode_view 		= dpd->Branch("tpc_totcharge_anode_view",		&_tpc_totcharge_anode_view, 	"tpc_totcharge_anode_view[2]/D");	
		TBranch * _bn_tpc_totrecocharge_anode_view 	= dpd->Branch("tpc_totrecocharge_anode_view", 	&_tpc_totrecocharge_anode_view, "tpc_totrecocharge_anode_view[2]/D");
		TBranch * _bn_tpc_totcharge_LAr_view 		= dpd->Branch("tpc_totcharge_LAr_view", 		&_tpc_totcharge_LAr_view, 		"tpc_totcharge_LAr_view[2]/D");	
		TBranch * _bn_tpc_totrecocharge_LAr_view 	= dpd->Branch("tpc_totrecocharge_LAr_view", 	&_tpc_totrecocharge_LAr_view, 	"tpc_totrecocharge_LAr_view[2]/D");
	}
	 
	TBranch * _bn_crt_daq_match					= dpd->Branch("crt_daq_match", 			&_crt_daq_match, 		"crt_daq_match/I");
	TBranch * _bn_crt_reco						= dpd->Branch("crt_reco", 				&_crt_reco, 			"crt_reco/I");
	TBranch * _bn_crt_reco_sat					= dpd->Branch("crt_reco_sat", 			&_crt_reco_sat, 		"crt_reco_sat/I");
	TBranch * _bn_crt_matchreco					= dpd->Branch("crt_matchreco", 			&_crt_matchreco, 		"crt_matchreco/I");
	
	TBranch * _bn_crt_isFC						= dpd->Branch("crt_isFC",				&_crt_isFC, 			"crt_isFC/I");
	TBranch * _bn_crt_isFV						= dpd->Branch("crt_isFV",				&_crt_isFV, 			"crt_isFV/I");
	TBranch * _bn_crt_isclosestpoint 			= dpd->Branch("crt_isclosestpoint",		&_crt_isclosestpoint, 	"crt_isclosestpoint[5]/I");
	TBranch * _bn_crt_ToF						= dpd->Branch("crt_ToF",				&_crt_ToF, 				"crt_ToF/F");
	TBranch * _bn_crt_pmt_dist_n				= dpd->Branch("crt_pmt_dist",			&_crt_pmt_dist, 		"crt_pmt_dist[5]/F");
	TBranch * _bn_crt_track_param				= dpd->Branch("crt_track_param",		&_crt_track_param, 		"crt_track_param[5]/F");
	//TBranch * _bn_crt_track_param_LArSoft;  dpd->Branch("crt_track_param_LArSoft"    , _crt_track_param_LArSoft      , "crt_track_param_LArSoft[2]/F"    );
	TBranch * _bn_crt_drift_len					= dpd->Branch("crt_drift_len",			&_crt_drift_len, 		"crt_drift_len[5]/F");
	TBranch * _bn_crt_closest_coord				= dpd->Branch("crt_closest_coord",		&_crt_closest_coord, 	"crt_closest_coord[5][3]/F");
	
	TBranch * _bn_crt_track_door				= dpd->Branch("crt_track_door",			&_crt_track_door, 		"crt_track_door[3]/F");
	TBranch * _bn_crt_track_wall				= dpd->Branch("crt_track_wall",			&_crt_track_wall, 		"crt_track_wall[3]/F");
	TBranch * _bn_crt_point_door_fc				= dpd->Branch("crt_point_door_fc",		&_crt_point_door_fc, 	"crt_point_door_fc[3]/F");
	TBranch * _bn_crt_point_wall_fc				= dpd->Branch("crt_point_wall_fc",		&_crt_point_wall_fc, 	"crt_point_wall_fc[3]/F");
	TBranch * _bn_crt_point_door_fv				= dpd->Branch("crt_point_door_fv",		&_crt_point_door_fv, 	"crt_point_door_fv[3]/F");
	TBranch * _bn_crt_point_wall_fv				= dpd->Branch("crt_point_wall_fv",		&_crt_point_wall_fv, 	"crt_point_wall_fv[3]/F");

	TBranch * _bn_crt_fit_yx_m					= dpd->Branch("crt_fit_yx_m",			&_crt_fit_yx_m,			"crt_fit_yx_m/F");
	TBranch * _bn_crt_fit_yx_n					= dpd->Branch("crt_fit_yx_n",			&_crt_fit_yx_n,			"crt_fit_yx_n/F");
	TBranch * _bn_crt_fit_zx_m					= dpd->Branch("crt_fit_zx_m",			&_crt_fit_zx_m,			"crt_fit_zx_m/F");
	TBranch * _bn_crt_fit_zx_n					= dpd->Branch("crt_fit_zx_n",			&_crt_fit_zx_n,			"crt_fit_zx_n/F");
	
	TBranch * _bn_pmt_wvf_properties 			= dpd->Branch("pmt_wvf_properties", 	_pmt_wvf_properties, 	"pmt_wvf_properties[5][2]/D");
	TBranch * _bn_pmt_charge 					= dpd->Branch("pmt_charge", 			_pmt_charge, 			"pmt_charge[6]/D");	
	TBranch * _bn_pmt_npe						= dpd->Branch("pmt_npe",				_pmt_npe, 				"pmt_npe[6]/D");
	TBranch * _bn_pmt_ped       				= dpd->Branch("pmt_ped",  				_pmt_ped, 				"pmt_ped[5]/D");
	TBranch * _bn_pmt_pedrms    				= dpd->Branch("pmt_pedrms", 			_pmt_pedrms, 			"pmt_pedrms[5]/D");
	TBranch * _bn_pmt_ped_uncorr     			= dpd->Branch("pmt_ped_uncorr",  		_pmt_ped_uncorr, 		"pmt_ped_uncorr[5]/D");
	TBranch * _bn_pmt_pedrms_uncorr  			= dpd->Branch("pmt_pedrms_uncorr", 		_pmt_pedrms_uncorr, 	"pmt_pedrms_uncorr[5]/D");
	TBranch * _bn_pmt_ped2       				= dpd->Branch("pmt_ped2",  				_pmt_ped2, 				"pmt_ped2[5]/D");
	TBranch * _bn_pmt_pedrms2    				= dpd->Branch("pmt_pedrms2", 			_pmt_pedrms2, 			"pmt_pedrms2[5]/D");
	TBranch * _bn_pmt_ped_end   				= dpd->Branch("pmt_ped_end",  			_pmt_ped_end, 			"pmt_ped_end[5]/D");
	TBranch * _bn_pmt_ped_end_uncorr   			= dpd->Branch("pmt_ped_end_uncorr", 	_pmt_ped_end_uncorr, 	"pmt_ped_end_uncorr[5]/D");
	TBranch * _bn_pmt_pedrms_end   				= dpd->Branch("pmt_pedrms_end",  		_pmt_pedrms_end, 		"pmt_pedrms_end[5]/D");
	TBranch * _bn_pmt_pedrms_end_uncorr   		= dpd->Branch("pmt_pedrms_end_uncorr",	_pmt_pedrms_end_uncorr, "pmt_pedrms_end_uncorr[5]/D");
	TBranch * _bn_pmt_wvf_ADC_sat  				= dpd->Branch("pmt_wvf_ADC_sat",  		_pmt_wvf_ADC_sat, 		"pmt_wvf_ADC_sat[5]/I");
	TBranch * _bn_pmt_wvf_end_fit_c0   			= dpd->Branch("pmt_wvf_end_fit_c0",		_pmt_wvf_end_fit_c0, 	"pmt_wvf_end_fit_c0[5]/D");
	TBranch * _bn_pmt_wvf_end_fit_c1   			= dpd->Branch("pmt_wvf_end_fit_c1",		_pmt_wvf_end_fit_c1, 	"pmt_wvf_end_fit_c1[5]/D");
	TBranch * _bn_pmt_wvf_end_fit_c2   			= dpd->Branch("pmt_wvf_end_fit_c2",		_pmt_wvf_end_fit_c2, 	"pmt_wvf_end_fit_c2[5]/D");
	TBranch * _bn_pmt_wvf_end_fit_chi2  		= dpd->Branch("pmt_wvf_end_fit_chi2",	_pmt_wvf_end_fit_chi2, 	"pmt_wvf_end_fit_chi2[5]/D");
	TBranch * _bn_pmt_wvf_end_fit_ndof  		= dpd->Branch("pmt_wvf_end_fit_ndof",	_pmt_wvf_end_fit_ndof, 	"pmt_wvf_end_fit_ndof[5]/I");
	TBranch * _bn_pmt_wvf_corr_RC_eff     		= dpd->Branch("pmt_wvf_corr_RC_eff",	_pmt_wvf_corr_RC_eff, 	"pmt_wvf_corr_RC_eff[5]/D");
	TBranch * _bn_pmt_wvf_max_corr    			= dpd->Branch("pmt_wvf_max_corr",		_pmt_wvf_max_corr, 		"pmt_wvf_max_corr[5]/D");
	TBranch * _bn_pmt_wvf_max_uncorr    		= dpd->Branch("pmt_wvf_max_uncorr",		_pmt_wvf_max_uncorr, 	"pmt_wvf_max_uncorr[5]/D");		
	TBranch * _bn_pmt_npeaks	   			 	= dpd->Branch("pmt_npeaks",  			_pmt_npeaks, 			"pmt_npeaks[5]/I");
	TBranch * _bn_pmt_peaks_tau 				= dpd->Branch("pmt_peaks_tau", 			_pmt_peaks_tau,  		Form("pmt_peaks_tau[5][%d]/D",KMAXNPEAKS));
	
	TBranch * _bn_pmt_S1_charge_1us				= dpd->Branch("pmt_S1_charge_1us",		_pmt_S1_charge_1us, 	"pmt_S1_charge_1us[6]/D");
	TBranch * _bn_pmt_S1_charge_4us				= dpd->Branch("pmt_S1_charge_4us",		_pmt_S1_charge_4us, 	"pmt_S1_charge_4us[6]/D");
	TBranch * _bn_pmt_S1_charge_80ns			= dpd->Branch("pmt_S1_charge_80ns",		_pmt_S1_charge_80ns, 	"pmt_S1_charge_80ns[6]/D");
	TBranch * _bn_pmt_S1_charge_m2				= dpd->Branch("pmt_S1_charge_m2",		_pmt_S1_charge_m2, 		"pmt_S1_charge_m2[6]/D");
	TBranch * _bn_pmt_S1_npe_1us				= dpd->Branch("pmt_S1_npe_1us", 		_pmt_S1_npe_1us, 		"pmt_S1_npe_1us[6]/D");
	TBranch * _bn_pmt_S1_npe_4us				= dpd->Branch("pmt_S1_npe_4us", 		_pmt_S1_npe_4us, 		"pmt_S1_npe_4us[6]/D");
	TBranch * _bn_pmt_S1_npe_80ns				= dpd->Branch("pmt_S1_npe_80ns", 		_pmt_S1_npe_80ns, 		"pmt_S1_npe_80ns[6]/D");
	TBranch * _bn_pmt_S1_npe_m2					= dpd->Branch("pmt_S1_npe_m2", 			_pmt_S1_npe_m2 , 		"pmt_S1_npe_m2[6]/D");
	TBranch * _bn_pmt_S1_width					= dpd->Branch("pmt_S1_width", 			_pmt_S1_width,	 		"pmt_S1_width[5]/D");
	TBranch * _bn_pmt_S1_amp					= dpd->Branch("pmt_S1_amp", 			_pmt_S1_amp, 			"pmt_S1_amp[5]/D");
	TBranch * _bn_pmt_S1_tau					= dpd->Branch("pmt_S1_tau",				_pmt_S1_tau, 			"pmt_S1_tau[5]/D");
	TBranch * _bn_pmt_S1_tau_end				= dpd->Branch("pmt_S1_tau_end",			_pmt_S1_tau_end, 		"pmt_S1_tau_end[5]/D");

	TBranch * _bn_pmt_S2_charge					= dpd->Branch("pmt_S2_charge",			_pmt_S2_charge, 		"pmt_S2_charge[6]/D");
	TBranch * _bn_pmt_S2_charge_m2				= dpd->Branch("pmt_S2_charge_m2",		_pmt_S2_charge_m2, 		"pmt_S2_charge_m2[6]/D");
	//TBranch * _bn_pmt_S2_charge_corr 			= dpd->Branch("pmt_S2_charge_corr",		_pmt_S2_charge_corr, 	"pmt_S2_charge_corr[6]/D");
	//TBranch * _bn_pmt_S2_charge_corr_2 		= dpd->Branch("pmt_S2_charge_corr_2",	_pmt_S2_charge_corr_2, "pmt_S2_charge_corr_2[6]/D");
	TBranch * _bn_pmt_S2_charge_corr_p 			= dpd->Branch("pmt_S2_charge_corr_p",	_pmt_S2_charge_corr_p, 	"pmt_S2_charge_corr_p[6]/D");
	TBranch * _bn_pmt_S2_charge_corr_m 			= dpd->Branch("pmt_S2_charge_corr_m",	_pmt_S2_charge_corr_m, 	"pmt_S2_charge_corr_m[6]/D");
	TBranch * _bn_pmt_S2_npe					= dpd->Branch("pmt_S2_npe", 			_pmt_S2_npe, 			"pmt_S2_npe[6]/D");
	TBranch * _bn_pmt_S2_npe_m2					= dpd->Branch("pmt_S2_npe_m2", 			_pmt_S2_npe_m2, 		"pmt_S2_npe_m2[6]/D");
	//TBranch * _bn_pmt_S2_npe_corr				= dpd->Branch("pmt_S2_npe_corr", 		_pmt_S2_npe_corr, 		"pmt_S2_npe_corr[6]/D"   );
	TBranch * _bn_pmt_S2_npe_corr_p				= dpd->Branch("pmt_S2_npe_corr_p", 		_pmt_S2_npe_corr_p, 	"pmt_S2_npe_corr_p[6]/D");
	TBranch * _bn_pmt_S2_npe_corr_m				= dpd->Branch("pmt_S2_npe_corr_m", 		_pmt_S2_npe_corr_m, 	"pmt_S2_npe_corr_m[6]/D");
	TBranch * _bn_pmt_S2_width					= dpd->Branch("pmt_S2_width", 			_pmt_S2_width, 			"pmt_S2_width[5]/D");
	TBranch * _bn_pmt_S2_amp					= dpd->Branch("pmt_S2_amp", 			_pmt_S2_amp , 			"pmt_S2_amp[5]/D");
	TBranch * _bn_pmt_S2_tau					= dpd->Branch("pmt_S2_tau", 			_pmt_S2_tau , 			"pmt_S2_tau[5]/D");
	TBranch * _bn_pmt_S2_tau_avg 				= dpd->Branch("pmt_S2_tau_avg", 		_pmt_S2_tau_avg, 		"pmt_S2_tau_avg[5]/D");
	TBranch * _bn_pmt_S2_tau_start				= dpd->Branch("pmt_S2_tau_start", 		_pmt_S2_tau_start, 		"pmt_S2_tau_start[5]/D");
	TBranch * _bn_pmt_S2_tau_end				= dpd->Branch("pmt_S2_tau_end",			_pmt_S2_tau_end, 		"pmt_S2_tau_end[5]/D");
	
	if (doCharge)
	{
		TBranch * _bn_tpc_track_id 					= dpd->Branch("tpc_track_id", 					_tpc_track_id,   				"tpc_track_id[ntracks]/S");
		TBranch * _bn_tpc_track_nhits_view 			= dpd->Branch("tpc_track_nhits_view", 			_tpc_track_nhits_view, 			"tpc_track_nhits_view[ntracks][2]/S");
		
		TBranch * _bn_tpc_track_charge				= dpd->Branch("tpc_track_charge",				&_tpc_track_charge, 			"tpc_track_charge[ntracks]/D");
		TBranch * _bn_tpc_track_charge_anode		= dpd->Branch("tpc_track_charge_anode",			&_tpc_track_charge_anode, 		"tpc_track_charge_anode[ntracks]/D");
		TBranch * _bn_tpc_track_charge_anode_view	= dpd->Branch("tpc_track_charge_anode_view",	&_tpc_track_charge_anode_view, 	"tpc_track_charge_anode_view[ntracks][2]/D");
		TBranch * _bn_tpc_track_charge_LAr_view		= dpd->Branch("tpc_track_charge_LAr_view", 		&_tpc_track_charge_LAr_view , 	"tpc_track_charge_LAr_view[ntracks][2]/D");
		TBranch * _bn_tpc_track_dE_view				= dpd->Branch("tpc_track_dE_view",				&_tpc_track_dE_view , 			"tpc_track_dE_view[ntracks][2]/D");
		TBranch * _bn_tpc_track_avg_dEds_view		= dpd->Branch("tpc_track_avg_dEds_view",		&_tpc_track_avg_dEds_view, 		"tpc_track_avg_dEds_view[ntracks][2]/D");
		
		TBranch * _bn_tpc_track_momentum 			= dpd->Branch("tpc_track_momentum", 			_tpc_track_momentum, 			"tpc_track_momentum[ntracks]/F");
		TBranch * _bn_tpc_track_length_trajectory 	= dpd->Branch("tpc_track_length_trajectory",	_tpc_track_length_trajectory, 	"tpc_track_length_trajectory[ntracks]/F");
		TBranch * _bn_tpc_track_length_straight_line= dpd->Branch("tpc_track_length_straight_line", _tpc_track_length_straight_line,"tpc_track_length_straight_line[ntracks]/F");
		TBranch * _bn_tpc_track_max_CBR_view		= dpd->Branch("tpc_track_max_CBR_view",			_tpc_track_max_CBR_view, 		"tpc_track_max_CBR_view[ntracks][2]/F");
		
		TBranch * _bn_tpc_track_start_x				= dpd->Branch("tpc_track_start_x",				_tpc_track_start_x, 			"tpc_track_start_x[ntracks]/F");
		TBranch * _bn_tpc_track_start_y				= dpd->Branch("tpc_track_start_y",				_tpc_track_start_y, 			"tpc_track_start_y[ntracks]/F");
		TBranch * _bn_tpc_track_start_z				= dpd->Branch("tpc_track_start_z",				_tpc_track_start_z, 			"tpc_track_start_z[ntracks]/F");
		TBranch * _bn_tpc_track_end_x				= dpd->Branch("tpc_track_end_x",				_tpc_track_end_x, 				"tpc_track_end_x[ntracks]/F");
		TBranch * _bn_tpc_track_end_y				= dpd->Branch("tpc_track_end_y",				_tpc_track_end_y, 				"tpc_track_end_y[ntracks]/F");
		TBranch * _bn_tpc_track_end_z				= dpd->Branch("tpc_track_end_z",				_tpc_track_end_z, 				"tpc_track_end_z[ntracks]/F");
		
		TBranch * _bn_tpc_track_start_theta			= dpd->Branch("tpc_track_start_theta",			_tpc_track_start_theta, 		"tpc_track_start_theta[ntracks]/F");
		TBranch * _bn_tpc_track_start_phi			= dpd->Branch("tpc_track_start_phi",			_tpc_track_start_phi, 			"tpc_track_start_phi[ntracks]/F");
		TBranch * _bn_tpc_track_end_theta			= dpd->Branch("tpc_track_end_theta",			_tpc_track_end_theta,			"tpc_track_end_theta[ntracks]/F");
		TBranch * _bn_tpc_track_end_phi				= dpd->Branch("tpc_track_end_phi",				_tpc_track_end_phi,				"tpc_track_end_phi[ntracks]/F");

		TBranch * _bn_tpc_track_fit_yx_m			= dpd->Branch("tpc_track_fit_yx_m",				_tpc_track_fit_yx_m, 			"tpc_track_fit_yx_m[ntracks]/F");
		TBranch * _bn_tpc_track_fit_yx_n			= dpd->Branch("tpc_track_fit_yx_n",				_tpc_track_fit_yx_n, 			"tpc_track_fit_yx_n[ntracks]/F");
		TBranch * _bn_tpc_track_fit_yx_chi2			= dpd->Branch("tpc_track_fit_yx_chi2",			_tpc_track_fit_yx_chi2, 		"tpc_track_fit_yx_chi2[ntracks]/F");
		TBranch * _bn_tpc_track_fit_yx_ndof			= dpd->Branch("tpc_track_fit_yx_ndof",			_tpc_track_fit_yx_ndof, 		"tpc_track_fit_yx_ndof[ntracks]/I");
		
		TBranch * _bn_tpc_track_fit_zx_m			= dpd->Branch("tpc_track_fit_zx_m",				_tpc_track_fit_zx_m, 			"tpc_track_fit_zx_m[ntracks]/F");
		TBranch * _bn_tpc_track_fit_zx_n			= dpd->Branch("tpc_track_fit_zx_n",				_tpc_track_fit_zx_n, 			"tpc_track_fit_zx_n[ntracks]/F");
		TBranch * _bn_tpc_track_fit_zx_chi2			= dpd->Branch("tpc_track_fit_zx_chi2",			_tpc_track_fit_zx_chi2, 		"tpc_track_fit_zx_chi2[ntracks]/F");
		TBranch * _bn_tpc_track_fit_zx_ndof			= dpd->Branch("tpc_track_fit_zx_ndof",			_tpc_track_fit_zx_ndof, 		"tpc_track_fit_zx_ndof[ntracks]/I");
		
		TBranch * _bn_tpc_drift_time_at_pmt			= dpd->Branch("tpc_drift_time_at_pmt", 			_tpc_drift_time_at_pmt, 		"tpc_drift_time_at_pmt[5]/F"       );
		TBranch * _bn_tpc_drift_time_at_pmt2		= dpd->Branch("tpc_drift_time_at_pmt2", 		_tpc_drift_time_at_pmt2, 		"tpc_drift_time_at_pmt2[5]/F"       );
	}
	
	
	// initialize vector for highway algorithm CBR
	std::vector<std::vector<std::pair<float,float>>> vvMaxCBR;	

	if (debug) cout << "Total number of events in tree: \t" << t2->GetEntries()  << endl;

	// file splitting for large light runs
	int startEvent=0;
	int lastEvent=t2->GetEntries();
	if (light_subrun>=0)
	{
		startEvent=light_subrun*1000;
		lastEvent=startEvent+TMath::Min(1000,lastEvent-startEvent);
	}
	
	for(int ev=startEvent; ev<lastEvent ; ++ev)
	{	
		cout << "dpdMaker: Event " << ev << " over " << (lastEvent-startEvent) << endl;
		
		if (debug) cout << "... Initializing DPD variables ..." << endl;
		initDPDvariables();

 		t2->GetEntry(ev);
		
		_trig_conf=trigConf;
				
		if (!doCharge) _runlight=runNum;
		
		if (debug)
		{
			if (doCharge) printf("Charge (run,subrun,event) = %d, %d, %d\n",_runcharge,_subruncharge,_event); 
			printf("Light (run,event) = %d, %d\n",_runlight,_nevent);
		}
		
		// speed up
		/*
		if (_crt_daq_match!=1 || _crt_reco!=1) continue;
		if (_crt_ToF<=0. || _crt_isFV!=1) continue;
		if (_crt_track_lenFV<=3100.) continue;
		*/
		
		// calculate deltatime 
		double time_charge=0.0;
		double time_light=0.0;
		
		// check for light run 1290: if PCTimeTag[0]==1, then use PCTimeTag[1]. also, PCTimeTAg[2] might not work!
		// might be true for light runs < 1320
		time_light=(double)_PCTimeTag[0]+(double)_PCTimeTag[2]*1e-3;
		
		
		//
		// Charge stuff
		//
		
		if (doCharge)
		{
			time_charge=(double)_eventtime_seconds-35.0+(double)_eventtime_nanoseconds*1e-9;
			
			if (debug)
			{
				cout << "\tTime: C_" << _eventtime_seconds << " L_" << _PCTimeTag[0] << endl;
				cout << "\tTime C_ns_" << _eventtime_nanoseconds << " L_ms_" << _PCTimeTag[2] << endl;
				cout << "\tdt: " <<time_light-time_charge << endl<< endl;
			}
			
			TH2F * view[2];
			view[0] = new TH2F("view0","View 0;Channel;Time",320,0,320,1667,0,1667);
			view[1] = new TH2F("view1","View 1;Channel;Time",960,0,960,1667,0,1667);

			double totq=0;
			double totrecoq=0;
			double trackcharge[KMAXNTRACKS]={0};
			int number_reco_hits=0;
	
			if (debug) cout << "Charge Event "<<_event <<" - n.hits "<<_NumberOfHits<< " - n.tracks " << _NumberOfTracks << endl;
			
			// JOSE: Total charge, total reco charge, track charge calculation using Hit_ChargeIntegral
			
			for (int j=0; j<_NumberOfHits; j++)
			{
				//if (debug) cout << "\t hit "<< j << " over " << _NumberOfHits << endl;
				totq+=_Hit_ChargeIntegral[j];
				if(_Hit_TrackID[j]!=-1)
				{
					//if (debug) cout << "Hit belongs to the track " << j << " "  << _Hit_TrackID[j]   << endl;
					totrecoq+=_Hit_ChargeIntegral[j];
					number_reco_hits++;
					trackcharge[_Hit_TrackID[j]]+=_Hit_ChargeIntegral[j];
				}
				//if (debug) cout << "\t\t hitview "<< _Hit_View[j]  << endl;
				if( _Hit_View[j]==0) view[_Hit_View[j]]->Fill(_Hit_Channel[j],0.4*(1667-_Hit_PeakTime[j]),_Hit_ChargeIntegral[j]);
				else view[_Hit_View[j]]->Fill(_Hit_Channel[j]-320,0.4*(1667-_Hit_PeakTime[j]),_Hit_ChargeIntegral[j]);
			}
			

			// MICHAEL + CHRISTOPH: Total charge, total reco charge, track charge calculation using Hit_ChargeSummedADC
			
			// charge collected at the anode
			double totcharge_anode[2]={0};
			double totrecocharge_anode[2]={0};
			double mytrackcharge_anode[KMAXNTRACKS][2]={0};
			
			// charge deposited in liquid (dQ_l)
			// divide anode_charge by effective gain = extraction efficiency * LEM gain * collection efficiency
			float  G_eff = getEffectiveGain(_runcharge);
			if (ev==startEvent) cout << "Effective gain for charge run " << _runcharge << ": " << G_eff << endl;
			double totcharge_liquid[2]={0};
			double totrecocharge_liquid[2]={0};
			double mytrackcharge_liquid[KMAXNTRACKS][2]={0};
			
			// dQ_l / ds
			// divide by either Track_Hit_ds_LocalTrackDirection or Track_Hit_ds_3DPosition (both work similarly well)
			double totrecocharge_liquid_ds[2]={0};
			double mytrackcharge_liquid_ds[KMAXNTRACKS][2]={0};
			
			// dE/ds
			double mytrack_dE[KMAXNTRACKS][2]={0};
			double mytrack_avg_dEds[KMAXNTRACKS][2]={0};
			
			// calculation of total charge 
			for (int j=0; j<_NumberOfHits; j++)
			{
				totcharge_anode[_Hit_View[j]]+=_Hit_ChargeSummedADC[j];
				
				if (G_eff>0) totcharge_liquid[_Hit_View[j]]+=2.*_Hit_ChargeSummedADC[j]/G_eff;
			}

			// calculation of total reco charge, track charge 

			int counter=0;
			//cout << "Run " << _runcharge << ", subrun " << _subruncharge << ", event " << _event << endl;
			for (int j=0; j<_NumberOfTracks; j++)
			{
				for (int k=0; k<_Track_NumberOfHits[j]; k++)
				//for (int k=0; k<(_Track_NumberOfHitsPerView[j][0]+_Track_NumberOfHitsPerView[j][1]); k++)
				{
					double charge = _Track_Hit_ChargeSummedADC[counter];
					int view = _Track_Hit_View[counter];
					double ds = _Track_Hit_ds_LocalTrackDirection[counter];
					
					totrecocharge_anode[view]+=charge;
					mytrackcharge_anode[j][view]+=charge; // /chargeCalib[0];
					
					if (G_eff<=0) { counter++; continue; }
					
					totrecocharge_liquid[view]+=2.*charge/G_eff;
					mytrackcharge_liquid[j][view]+=2.*charge/G_eff;
					
					// smaller than channel pitch or larger than view length 
					if (ds < 0.3125 || (view==0 && ds>3000.) || (view==1 && ds>1000.)) { counter++; continue; }
					
					//if (_Track_Hit_ds_3DPosition[counter]<=0) { counter++; continue; }
					
					totrecocharge_liquid_ds[view]+=2.*charge/G_eff/ds;
					mytrackcharge_liquid_ds[j][view]+=2.*charge/G_eff/ds;
					
					double dQdsl = 2.*charge*1.E15/chargeCalib[view]/G_eff/ds;
					double dEds = 1./ ( (BirksA/(dQdsl*WorkFunctionMeVperfC)) - (BirksKDividedByLArDensity/drift_field) );
					
					if (dEds>0.) 
					{
						mytrack_dE[j][view] += ds*dEds;
						mytrack_avg_dEds[j][view] += dEds;
					}
					
					counter++;
				} // track loop
	
				// divide each view by charge calibration or total number of hits
				for (int v=0; v<2; v++)
				{
					mytrackcharge_anode[j][v]/=chargeCalib[v];
					mytrackcharge_liquid[j][v]/=chargeCalib[v];
					mytrackcharge_liquid_ds[j][v]/=chargeCalib[v];
					
					mytrack_avg_dEds[j][v]/=_Track_NumberOfHitsPerView[j][v];
				}
			}
			
			// divide each view by charge calibration
			for (int v=0; v<2; v++)
			{
				totcharge_anode[v]/=chargeCalib[v];
				totcharge_liquid[v]/=chargeCalib[v];
				totrecocharge_anode[v]/=chargeCalib[v];
				totrecocharge_liquid[v]/=chargeCalib[v];
				totrecocharge_liquid_ds[v]/=chargeCalib[v];
			}

			if (debug)
			{
				cout << "Charge Event "<<_event << " Reco proportion (" <<Form("%.0d%s",100*number_reco_hits/_NumberOfHits,"%") <<" hits) - "<< Form("%.2f%s",100*totrecoq/totq,"%") << "IADC"  <<  endl;
				
				cout << "anode view 0 = " << totcharge_anode[0] << ", anode view 1 = " << totcharge_anode[1] << endl;
				cout << "liquid view 0 = " << totcharge_liquid[0] << ", liquid view 1 = " << totcharge_liquid[1] << endl;
				cout << "reco anode view 0 = " << totrecocharge_anode[0] << ", reco anode view 1 = " << totrecocharge_anode[1] << endl;
				cout << "reco liquid view 0 = " << totrecocharge_liquid[0] << ", reco liquid view 1 = " << totrecocharge_liquid[1] << endl;
				cout << "reco liquid ds view 0 = " << totrecocharge_liquid_ds[0] << ", reco liquid ds view 1 = " << totrecocharge_liquid_ds[1] << endl;
				
				for (int j=0; j<_NumberOfTracks; j++)
				{
					cout << "Track " << j << " hits in view 0 "<< _Track_NumberOfHitsPerView[j][0] <<  endl;
					cout << "Track " << j << " hits in view 1 "<< _Track_NumberOfHitsPerView[j][1] <<  endl;
					cout << "Total Charge per track " << j << ": "<< trackcharge[j] <<  endl;
				}
			}



			// TRACK SORTING

			if (debug) cout << "\t... sorting 3d tracks by their charge." << endl;

			short sorted_tracks[KMAXNTRACKS]; // tracks sorted by total charge in the track

			std::size_t n(0);
			std::generate(std::begin(sorted_tracks), std::end(sorted_tracks), [&]{ return n++; });

			// JOSE: sort by track_charge
			//std::sort(  std::begin(sorted_tracks), std::end(sorted_tracks), [&](int i1, int i2) { return trackcharge[i1] > trackcharge[i2]; } );
			
			// Michael: sort by track_charge_anode
			std::sort(  std::begin(sorted_tracks), std::end(sorted_tracks), [&](int i1, int i2) { 
				return (mytrackcharge_anode[i1][0]+mytrackcharge_anode[i1][1])>(mytrackcharge_anode[i2][0]+mytrackcharge_anode[i2][1]); } );
				
			// Chiara: sort by track length
			//std::sort(  std::begin(sorted_tracks), std::end(sorted_tracks), [&](int i1, int i2) { return _Track_Length_Trajectory[i1]>_Track_Length_Trajectory[i2]; } );

			//for (int j=0; j<_NumberOfTracks; j++) cout << "\t\tSorted track " << j << " " << sorted_tracks[j] << " " << trackcharge[sorted_tracks[j]] <<"IAC"<< endl;
		
			if (debug) cout << "\t... 3d tracks sorted." << endl;
				

				
				
				
			
			// HIGHWAY ALGORITHM
			
			if (debug) cout << "\t... reading Highway algorithm information ..." << endl;
			
			float _Track_MaxCBR_view[KMAXNTRACKS][2]={0};
			if (vvMaxCBR.size()==0) vvMaxCBR = HighwayReader(_runcharge,_subruncharge);
			for (int j=0; j<_NumberOfTracks; j++)
			{
				if (vvMaxCBR.at(ev).size()>0 && j<vvMaxCBR.at(ev).size()) 
				{
					_Track_MaxCBR_view[j][0] = (vvMaxCBR.at(ev).at(j).first!=0)?vvMaxCBR.at(ev).at(j).first:ERRVAL;
					_Track_MaxCBR_view[j][1] = (vvMaxCBR.at(ev).at(j).second!=0)?vvMaxCBR.at(ev).at(j).second:ERRVAL;
				}
				else
				{
					_Track_MaxCBR_view[j][0] = ERRVAL;
					_Track_MaxCBR_view[j][1] = ERRVAL;
				}
				//cout << "Track " << j << ": " << _Track_MaxCBR_view[j][0] << ", " << _Track_MaxCBR_view[j][1] << endl;
			}
			
			if (debug) cout << "\t... Highway reader done ..." << endl;
		
		
		
			// TRACK SELECTION
			
			if (debug) cout << "\t ... selecting tracks ..." << endl;
			
			int isSelected[KMAXNTRACKS]={-1};
			short NumberOfSelectedTracks=0;
			for (int j=0; j<_NumberOfTracks; j++)
			{
				/*
				// track is greater than 100 cm and reconstructed end-points 
				// must be closer than 2 cm to anode and cathode
					if (_Track_Length_StraightLine[j]>100. && 
						((fabs(_Track_StartX[j]-50.)<2. && fabs(_Track_EndX[j]+50.)<2.) || 
						 (fabs(_Track_StartX[j]+50.)<2. && fabs(_Track_EndX[j]-50.)<2.)))
								NumberOfSelectedTracks++;
				*/
				if (_Track_Length_StraightLine[j]>100. && _Track_StartDirection_Theta[j]>92.)
					//&& fabs(_Track_StartDirection_Phi[j])>5. && fabs(_Track_StartDirection_Phi[j])<85. )
				{
					NumberOfSelectedTracks++;
					isSelected[j]=true;
					
					printf("Track %d (selected), view 0: dE = %f MeV, dL = %f cm, dE/dL = %f MeV/cm, <dE/ds> = %f\n",j,mytrack_dE[j][0],_Track_Length_Trajectory[j],mytrack_dE[j][0]/_Track_Length_Trajectory[j],mytrack_avg_dEds[j][0]);
					printf("Track %d (selected), view 1: dE = %f MeV, dL = %f cm, dE/dL = %f MeV/cm, <dE/ds> = %f\n",j,mytrack_dE[j][1],_Track_Length_Trajectory[j],mytrack_dE[j][1]/_Track_Length_Trajectory[j],mytrack_avg_dEds[j][1]);
				}
				else isSelected[j]=false;
				
			} 

			if (debug) cout << "\t ... done selecting tracks ..." << endl;



			// TRACK FITTING
			
			if (debug) cout << "\t ... Fitting tracks ..." << endl;
	
			TGraph *trackYX[KMAXNTRACKS], *trackZX[KMAXNTRACKS];
			double t_yx_m[KMAXNTRACKS], t_yx_n[KMAXNTRACKS], t_yx_chi2[KMAXNTRACKS];
			double t_zx_m[KMAXNTRACKS], t_zx_n[KMAXNTRACKS], t_zx_chi2[KMAXNTRACKS];
			int t_yx_ndof[KMAXNTRACKS], t_zx_ndof[KMAXNTRACKS];
			
			int b=0;
			for (int j=0; j<_NumberOfTracks; j++)
			{
				trackYX[j]= new TGraph();
				trackZX[j]= new TGraph();
				
				for (int l=b; l < (b+_Track_NumberOfHits[j]); l++) //Hit loop (for this track)  
				{
					trackYX[j]->SetPoint(trackYX[j]->GetN(),_Track_HitY[l],_Track_HitX[l]);
					trackZX[j]->SetPoint(trackZX[j]->GetN(),_Track_HitZ[l],_Track_HitX[l]);
				} 
				
				b+=_Track_NumberOfHits[j]; 
				
				//LeastSquareLinearFit(Int_t n, Double_t& constant, Double_t& slope, Int_t& ifail, Double_t xmin = 0, Double_t xmax = 0)
				//trackYX[j]->LeastSquareLinearFit(trackYX[j]->GetN()-1,t_yx_n[j],t_yx_m[j],fail_yx[j],-50,50);
				//trackZX[j]->LeastSquareLinearFit(trackZX[j]->GetN()-1,t_zx_n[j],t_zx_m[j],fail_zx[j],0,300);   
				
				
				TFitResultPtr p_yx = trackYX[j]->Fit("pol1","QRS","",-51,51);
				t_yx_n[j] = p_yx->Value(0);
				//t_yx_n_err[j] = p_yx->ParError(0);
				t_yx_m[j] = p_yx->Value(1);
				t_yx_chi2[j] = p_yx->Chi2();
				t_yx_ndof[j] = p_yx->Ndf();
				
				TFitResultPtr p_zx = trackZX[j]->Fit("pol1","QRS","",-1,301);
				t_zx_n[j] = p_zx->Value(0);
				t_zx_m[j] = p_zx->Value(1);
				t_zx_chi2[j] = p_zx->Chi2();
				t_zx_ndof[j] = p_zx->Ndf();
				
				
				if (display_tracks && isSelected[j])
				{
					cout << "Drawing Track " << j << ": " << endl;
					cout << "\tTrackYX m,n,chi2,ndof: "<< t_yx_m[j] << "  " << t_yx_n[j] << "  "  << t_yx_chi2[j] << "  "  << t_yx_ndof[j] << endl;
					cout << "\tTrackZX m,n,chi2,ndof: "<< t_zx_m[j] << "  " << t_zx_n[j] << "  "  << t_zx_chi2[j] << "  "  << t_zx_ndof[j] << endl;
					cout << "\tNumber of points in track " << j << ": "<< trackYX[j]->GetN() << ", " << trackZX[j]->GetN() << " totq = " <<trackcharge[j]<<endl;
					
					TCanvas* c5 = new TCanvas("c5","c5",1200,600);
					c5->Divide(2,1);
					c5->cd(1);
					trackYX[j]->Draw("AP");
					//TLine lin1 = TLine(-50.,t_yx_n[j]-50*t_yx_m[j],50.,t_yx_n[j]+50.*t_yx_m[j]);
					//lin1.SetLineColor(kRed);
					//lin1.Draw("same");
			
					c5->cd(2);
					trackZX[j]->Draw("AP");
					//TLine lin2 = TLine(0.,t_zx_n[j],300.,t_zx_n[j]+300.*t_zx_m[j]);
					//lin2.SetLineColor(kRed);
					//lin2.Draw("same");
			
					c5->Modified();
					c5->Update();
					lets_pause();
				
					delete c5;
				}
			}
			
			if (debug) cout << "\t ... done fitting tracks ..." << endl;
			
			
			// TPC DRIFT TIME AT PMT POSITION
			// for leading (highest-charge) track in the event
			
			if (debug) cout << "\t ... calculating TPC drift time ..." << endl;
			
			double S2_width_tolerance_channel=5.;//cm
			
			// JOSE's calculation (normal average)
			
			float tpc_drift_time_at_pmt[N_PMT]={0};
			int tpc_drift_time_at_pmt_counter[N_PMT]={0};
			
			for (int j=0; j<_NumberOfHits; j++)
			{			
				if(_Hit_ChargeIntegral[j]==_Hit_ChargeIntegral[j])
				{
					for (int k=0; k<N_PMT; k++)
					{
						float hitpos = 300.0*(_Hit_Channel[j]-320)/960.0;

						if (_Hit_TrackID[j] == sorted_tracks[0] && _Hit_View[j]==1 && 
							fabs(hitpos - pmt_z_pos[k]) < S2_width_tolerance_channel )
						{
							tpc_drift_time_at_pmt[k]+=0.4*_Hit_PeakTime[j];
							tpc_drift_time_at_pmt_counter[k]++;
							//cout << "channel " <<  _Hit_Channel[j]<< " in cm " << (_Hit_Channel[j]-320)*300.0/960.0<< " time " << 0.4*_Hit_PeakTime[j] << " assing to pmt " << k << endl;
						}
					}
				}
			}
			
			
			// MICHAEL's calculation (charge-weighted average)
			
			float tpc_drift_time_at_pmt2[N_PMT]={0};
			float tpc_drift_time_at_pmt_counter2[N_PMT]={0};
			
			b=0;
			for (int j=0; j<_NumberOfTracks; j++)
			{	
				if (j==sorted_tracks[0])
				{
					for (int l=b; l < (b+_Track_NumberOfHits[j]); l++) //Hit loop (for this track)  
					{
						if (_Track_Hit_View[l]==1)
						{
							for (int k=0; k<N_PMT; k++)
							{
								if (fabs(_Track_HitZ[l]-pmt_z_pos[k]) < S2_width_tolerance_channel)
								{
									tpc_drift_time_at_pmt2[k]+=0.4*_Track_Hit_ChargeSummedADC[l]*_Track_Hit_PeakTime[l];
									tpc_drift_time_at_pmt_counter2[k]+=_Track_Hit_ChargeSummedADC[l];
									//cout << "hit " << counter << ", z pos = " << z_hit << " cm, time = " << time_hit << " us, assigned to PMT " << k << endl;
								}
							}
						}
					} 
				}
				b+=_Track_NumberOfHits[j]; 
			}

			for (int k=0; k<N_PMT; k++)
			{
				if (tpc_drift_time_at_pmt_counter[k]!=0) tpc_drift_time_at_pmt[k]/=tpc_drift_time_at_pmt_counter[k];
				else tpc_drift_time_at_pmt[k]=ERRVAL;
				
				if (tpc_drift_time_at_pmt_counter2[k]!=0) tpc_drift_time_at_pmt2[k]/=tpc_drift_time_at_pmt_counter2[k];
				else tpc_drift_time_at_pmt2[k]=ERRVAL;
			}
			
			if (debug) cout << "\t ... done TPC drift time ..." << endl;
		
			if (debug)
			{
				cout << "Time: charge: " << time_charge << ", light: " << time_light << ", delta: " << time_charge - time_light << endl;
				cout << "Charge (Run,SubRun,ev): (" <<_runcharge << ","<< _subruncharge<<","<<_event<<") /t totcharge: "<<  totq <<endl; 
				cout << "   Length: (";
				for (Int_t i=0;i<_NumberOfTracks;i++) cout << _Track_Length_StraightLine[i] << ", ";
				cout << ")"<< endl;
				cout << "   Charge: (";
				for (Int_t i=0;i<_NumberOfTracks;i++) cout << trackcharge[i] << ", ";
				cout << ")"<< endl;
			}
		
			if (debug) cout << "\t ... preparing charge DPD variables ..." << endl;
			// prepare all charge DPD variables for filling later
		
			_run_charge_out=_runcharge;
			_subrun_charge_out=_subruncharge;
			_nevent_charge_out=_event;
		
			_ntracks=_NumberOfTracks;
			_n_sel_tracks=NumberOfSelectedTracks;

			_time_charge = time_charge;
		
			_tpc_totcharge=totq;
			_tpc_totrecocharge=totrecoq;
			
			_tpc_totcharge_anode=0;
			_tpc_totrecocharge_anode=0;
			for (int j=0; j<_NumberOfTracks; j++) _tpc_track_charge_anode[j]=0;
			
			for (int v=0; v<2; v++)
			{
				_tpc_totcharge_anode += totcharge_anode[v];
				_tpc_totrecocharge_anode += totrecocharge_anode[v];
				
				_tpc_totcharge_anode_view[v] = totcharge_anode[v];
				_tpc_totrecocharge_anode_view[v] = totrecocharge_anode[v];
				_tpc_totcharge_LAr_view[v] = totcharge_liquid[v];
				_tpc_totrecocharge_LAr_view[v] = totrecocharge_liquid[v];
				
				for (int j=0; j<_NumberOfTracks; j++)
				{
					_tpc_track_charge_anode[j] += mytrackcharge_anode[sorted_tracks[j]][v];
					
					_tpc_track_nhits_view[j][v] = _Track_NumberOfHitsPerView[sorted_tracks[j]][v];
					_tpc_track_charge_anode_view[j][v] = mytrackcharge_anode[sorted_tracks[j]][v];
					_tpc_track_charge_LAr_view[j][v] = mytrackcharge_liquid[sorted_tracks[j]][v];
					_tpc_track_dE_view[j][v] = mytrack_dE[sorted_tracks[j]][v];
					_tpc_track_avg_dEds_view[j][v] = mytrack_avg_dEds[sorted_tracks[j]][v];
					_tpc_track_max_CBR_view[j][v] = _Track_MaxCBR_view[sorted_tracks[j]][v];
				}
			}

			for (int j=0; j<_NumberOfTracks; j++) 
			{
				_tpc_track_charge[j] = trackcharge[sorted_tracks[j]];
				_tpc_track_fit_yx_m[j] = t_yx_m[sorted_tracks[j]];
				_tpc_track_fit_yx_n[j] = t_yx_n[sorted_tracks[j]];
				_tpc_track_fit_yx_chi2[j]=  t_yx_chi2[sorted_tracks[j]];
				_tpc_track_fit_yx_ndof[j]=  t_yx_ndof[sorted_tracks[j]];
				_tpc_track_fit_zx_m[j] = t_zx_m[sorted_tracks[j]];
				_tpc_track_fit_zx_n[j] = t_zx_n[sorted_tracks[j]];
				_tpc_track_fit_zx_chi2[j]=  t_zx_chi2[sorted_tracks[j]];
				_tpc_track_fit_zx_ndof[j]=  t_zx_ndof[sorted_tracks[j]];
				_tpc_track_start_x[j] = _Track_StartX[sorted_tracks[j]];
				_tpc_track_start_y[j] = _Track_StartY[sorted_tracks[j]];
				_tpc_track_start_z[j] = _Track_StartZ[sorted_tracks[j]];
				_tpc_track_end_x[j] = _Track_EndX[sorted_tracks[j]];
				_tpc_track_end_y[j] = _Track_EndY[sorted_tracks[j]];
				_tpc_track_end_z[j] = _Track_EndZ[sorted_tracks[j]];
				_tpc_track_start_theta[j] = _Track_StartDirection_Theta[sorted_tracks[j]];
				_tpc_track_start_phi[j] = _Track_StartDirection_Phi[sorted_tracks[j]];
				_tpc_track_end_theta[j] = _Track_EndDirection_Theta[sorted_tracks[j]];
				_tpc_track_end_phi[j] = _Track_EndDirection_Phi[sorted_tracks[j]];
				_tpc_track_id[j] = _TrackID[sorted_tracks[j]];
				_tpc_track_momentum[j] = _Track_Momentum[sorted_tracks[j]];
				_tpc_track_length_trajectory[j] = _Track_Length_Trajectory[sorted_tracks[j]];
				_tpc_track_length_straight_line[j] = _Track_Length_StraightLine[sorted_tracks[j]];
			}			
		
			for (int k=0; k<N_PMT; k++)
			{
				_tpc_drift_time_at_pmt[k] = tpc_drift_time_at_pmt[k];
				_tpc_drift_time_at_pmt2[k] = tpc_drift_time_at_pmt2[k];
			}
		
			if (view[0]) view[0]->Delete();
			if (view[1]) view[1]->Delete();
		
			for (int j=0; j<_NumberOfTracks; j++)
			{
				if (trackYX[j]) trackYX[j]->Delete();
				if (trackZX[j]) trackZX[j]->Delete();
			}
		}
		
		
		
		//
		// Light stuff
		//
		
 		TH1F *h[N_PMT];
		TH1F *h_plot[N_PMT];
		//TH1F *h_corr[N_PMT];
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
		//double ped_corr[N_PMT]={0};
		//double pedrms_corr[N_PMT]={0};
		double ped_uncorr[N_PMT]={0};
		double pedrms_uncorr[N_PMT]={0};
		double ped2[N_PMT]={0};
		double pedrms2[N_PMT]={0};
		double ped_end[N_PMT]={0};
		double pedrms_end[N_PMT]={0};
		double ped_end_uncorr[N_PMT]={0};
		double pedrms_end_uncorr[N_PMT]={0};
		double wvf_end_fit_c0[N_PMT]={0};
		double wvf_end_fit_c1[N_PMT]={0};
		double wvf_end_fit_c2[N_PMT]={0};
		double wvf_end_fit_chi2[N_PMT]={0};
		int    wvf_end_fit_ndof[N_PMT]={0};
		double wvf_corr_RC_eff[N_PMT]={0};
		double wvf_max_corr[N_PMT]={0};
		double wvf_max_uncorr[N_PMT]={0};

		
		//so let's do 0-0.4 for 4 us and 228-228.4 for 1 ms
		
		double ped_start = 228.0; // 218.;
		double ped_stop  = 228.4; // 228.;
		double ped2_start = 228.6;
		double ped2_stop  = 229.0;
		double ped_end_start = 910.;
		double ped_end_stop  = 912.; // 920.;
		
		int rebinfactor = 64;
		double S1_mintime = 228;
		double S1_maxtime = 231;
		double S2_maxtime = 1000.; // 900.;	
		
		if (_nsamples<=1000)
		{
			ped_start = 0.0;
			ped_stop  = 0.4;
			ped2_start = 0.0;
			ped2_stop  = 0.4;
			ped_end_start = 3.6;
			ped_end_stop  = 4.0;
		
			rebinfactor = 40;
			S1_mintime = 0.3;
			S1_maxtime = 0.85;
			S2_maxtime = 4.0;	
		}
		
		vector<int> pmt_valleys[N_PMT];
		double pmt_valleys_tau[N_PMT][KMAXNPEAKS]={0};
		
		for (int k=0; k<N_PMT; k++)	
		{	
			if (debug) cout << "\t... calculating pedestals"<< endl;
			
			ped[k] = WaveformAnalysis::baseline(h[k],pedrms[k],h[k]->FindBin(ped_start),h[k]->FindBin(ped_stop)-1); 
			
			ped2[k] = WaveformAnalysis::baseline(h[k],pedrms2[k],h[k]->FindBin(ped2_start),h[k]->FindBin(ped2_stop)-1); 
			
			ped_end[k] = WaveformAnalysis::baseline(h[k],pedrms_end[k],h[k]->FindBin(ped_end_start),h[k]->FindBin(ped_end_stop)-1); 
			
			double ped_diff = ped_end[k]-ped[k];
			
			
			// rebin histo for calculating binmax and fitting
			TH1F* htmp_rebin = dynamic_cast<TH1F*>(h[k]->Rebin(rebinfactor,Form("%s%i","htmp_rebin",k)));
			htmp_rebin->Scale(1.0/rebinfactor);
	
			// find bin maximum of rebinned original histo
			int orig_binmax = WaveformAnalysis::find_S2_binmax(htmp_rebin, S1_mintime, S2_maxtime);
			wvf_max_uncorr[k] = hget(htmp_rebin,orig_binmax)-ped[k];


			// WF fitting and correction for positive-base PMTs
			if (_nsamples>1000 && (k==2 || k==3)) 
			{
				if (debug) cout << "\t... fitting end of wf"<< endl;
				
				float RC_c_max = 1200.; 
			
				TF1* fit = new TF1("fit","[0]+[1]*exp((900-x)/[2])",850,1100);
				fit->SetParameter(0,ped[k]);
				fit->SetParameter(1,ped_diff);
				fit->SetParameter(2,250.);
				fit->SetParLimits(0,0.9*ped[k],1.1*ped[k]);
				fit->SetParLimits(1,0.9*ped_diff,1.1*ped_diff);
				fit->SetParLimits(2,10.,RC_c_max);
				htmp_rebin->Fit(fit,"QRMN","",900,1040);
				wvf_end_fit_c0[k] = fit->GetParameter(0);
				wvf_end_fit_c1[k] = fit->GetParameter(1);
				wvf_end_fit_c2[k] = fit->GetParameter(2);
				wvf_end_fit_chi2[k] = fit->GetChisquare();
				wvf_end_fit_ndof[k] = fit->GetNDF();
			
				if (debug) printf("Fit results, chan %d: c0 = %f, c1 = %f, c2 = %f, chi2 = %f, pedend_rms = %f\n",k,wvf_end_fit_c0[k],wvf_end_fit_c1[k],wvf_end_fit_c2[k],wvf_end_fit_chi2[k],pedrms_end[k]);
			
				// default case
				double RC_c_def=284.7;
				double RC_d_def=254.5;
				double Cap=200.E-9; // F
							
				double RC_c_cap = RC_c_def;
				double RC_d_cap = RC_d_def;
				
				if (_runlight < 609) // 6.8 nF
				{
					Cap = 6.8E-9;
					RC_c_cap = 82.9;
					RC_d_cap = 77.0;
				}
				else if (_runlight < 1064) // 150 nF
				{
					// note that ch 3 splitter disconnected for light runs 935-1063
					Cap = 150.E-9;
					RC_c_cap = 310.7;
					RC_d_cap = 297.8;
				}
				else if (k==3) // 300 nF
				{
					Cap = 300.E-9;
					RC_c_cap = 307.2;
					RC_d_cap = 294.6;
				}
				else if (_runlight < 1717) // 150 nF
				{
					Cap = 150.E-9;
					RC_c_cap = 310.7;
					RC_d_cap = 297.8;
				}
				else // 200 nF
				{
					Cap = 200.E-9;
					RC_c_cap = 260.4;
					RC_d_cap = 249.3;
				}
				
				if (debug) printf("Chan %d: Cap = %e, RC_c_cap = %f, RC_d_cap = %f\n",k,Cap,RC_c_cap,RC_d_cap);
				
				/*
				TH1F* htmp_rb = dynamic_cast<TH1F*>(h[k]->Rebin(rebinfactor,"htmp_rb"));
				htmp_rb->Draw("hist");
				lets_pause();
				delete htmp_rb;
				*/
				
				// set h_corr, h_corr_p and h_corr_m
				TH1F* h_corr = (TH1F*)h[k]->Clone("hcorr");
				h_corr_p[k] = (TH1F*)h[k]->Clone("hcorr_p");
				h_corr_m[k] = (TH1F*)h[k]->Clone("hcorr_m");
				
				// first do the easy variations on correction
				WaveformAnalysis::correct_wvf_histo(h[k],h_corr_p[k],ped[k],RC_c_cap,RC_d_cap);
				WaveformAnalysis::correct_wvf_histo(h[k],h_corr_m[k],ped[k],RC_c_def,RC_d_def);	
				
				// Now for the nominal correction
				// Scan RC_c parameter to find value that minimizes the maximum value of waveform during S2
				
				// set values to default in case there is error
				double RC_c = RC_c_def;
				double RC_d = RC_d_def;
				
				double MyBestMax=-9999;
				double d_factor = 1.0; // 0.925;
				
				// coarse search
				RC_c = WaveformAnalysis::calc_optimal_RC(h[k],ped[k],rebinfactor,orig_binmax,d_factor,50.,RC_c_max,50.,MyBestMax);
				
				if (debug || 1) cout << "chan " << k << ", best result step 1: RC_c = " << RC_c << ", RC_d = " << d_factor*RC_c << ", bestmax = " << MyBestMax << endl;
	
				// fine search
				RC_c = WaveformAnalysis::calc_optimal_RC(h[k],ped[k],rebinfactor,orig_binmax,d_factor,TMath::Max(0.,RC_c-40),TMath::Min((double)RC_c_max,RC_c+40),10.,MyBestMax);
				
				if (debug || 1) cout << "chan " << k << ", best result step 2: RC_c = " << RC_c << ", RC_d = " << d_factor*RC_c << ", bestmax = " << MyBestMax << endl;
				
				wvf_corr_RC_eff[k] = RC_c;
				wvf_max_corr[k] = MyBestMax;
				
				if (debug) cout << "\t... correcting wf histogram and re-calculating pedestal" << endl;
				
				
				// OVERWRITE original histos if wvf correction succeeds
				if (RC_c>=100. && RC_c<=700.) 
				{
					if (debug) cout << "Applying wvf correction for channel " << k << endl;
					WaveformAnalysis::correct_wvf_histo(h[k],h_corr,ped[k],RC_c,d_factor*RC_c);	
					
					// save uncorrected pedestal values
					ped_uncorr[k]=ped[k];
					pedrms_uncorr[k]=pedrms[k];
					ped_end_uncorr[k]=ped_end[k];
					pedrms_end_uncorr[k]=pedrms_end[k];
					
					h[k]->Reset();
					h[k]->Clear();
					h[k]->Add(h_corr);
					h[k]->SetTitle(Form("%s%i%s","Channel ",k," corrected;Time [#mus]; Voltage [ADC counts]"));
					h[k]->SetLineColor(kGreen+2);
					
					// recalculate pedestals
					ped[k] = WaveformAnalysis::baseline(h[k],pedrms[k],h[k]->FindBin(ped_start),h[k]->FindBin(ped_stop)-1); 
					ped2[k] = WaveformAnalysis::baseline(h[k],pedrms2[k],h[k]->FindBin(ped2_start),h[k]->FindBin(ped2_stop)-1); 
					ped_end[k] = WaveformAnalysis::baseline(h[k],pedrms_end[k],h[k]->FindBin(ped_end_start),h[k]->FindBin(ped_end_stop)-1); 
				}
				else 
				{
					if (debug) cout << "Not applying wvf correction for channel " << k << endl;
					h_corr->Reset();
					h_corr->Clear();
					h_corr->Add(h[k]);
				}
				
				
				/*
				double DC_c,DC_d;
				if (k==2 || k==3) WaveformAnalysis::calc_wvf_DC(h[k],ped[k],pedrms[k],DC_c,DC_d);
				*/	
				
				//pedrms_end[k]<2.0 && !pmt_wvf_ADC_sat[k] && wvf_end_fit_chi2[k]<4.0 && wvf_end_fit_c2[k]>20. && wvf_end_fit_c2[k]<1000
					
					
				if (0 && pedrms[k]<1.4) //  && _tpc_track_length_straight_line[0]>50.)
				{
					printf("evt %d: ped[%i] = %0.2f +/- %0.2f\n", ev, k, ped[k], pedrms[k]); 
					if (RC_c>10. && RC_c<RC_c_max) printf("ped_uncorr = %f, ped_end_uncorr = %f, ped_uncorr_diff = %f\n",ped_uncorr[k],ped_end_uncorr[k],ped_end_uncorr[k]-ped_uncorr[k]);
					printf("ped = %f, ped_diff = %f\n",ped[k],ped_end[k]-ped[k]);
					
					cout << endl;
				
					gStyle->SetOptStat(0);
					TCanvas *c1 = new TCanvas("c1");
					
					//h_plot[k]->GetXaxis()->SetRange(h_plot[k]->FindBin(650),h_plot[k]->FindBin(1048));
					htmp_rebin->Draw("hist");
					h_corr->SetLineColor(kGreen+2);
					h_corr->Rebin(rebinfactor);
					h_corr->Scale(1.0/rebinfactor);
					h_corr->Draw("hist.same");
					fit->SetLineColor(kRed);
					fit->Draw("same");
					
					float x = 0.55;
					float y = 0.45;
					TLatex l;
					l.SetNDC();
					l.DrawLatex(x,y+0.06,Form("Run %d, Event %d",_run_charge_out,_nevent_charge_out));
					//l.DrawLatex(x,y,Form("ped = %0.2f +/- %0.2f",ped[k],pedrms[k]));
					l.DrawLatex(x,y-0.06,Form("RC_fit = %0.1f #mus",wvf_end_fit_c2[k]));
					l.DrawLatex(x,y-0.12,Form("RC_eff = %0.1f #mus",RC_c));					
					
					c1->Modified();
					c1->Update();
					
					lets_pause();
				}

				delete h_corr;
				delete fit;				
			}
			
			delete htmp_rebin;
			
			if (debug) cout << "\t... rebinning wf histogram for plotting and fitting"<< endl;
			
			h_plot[k] = dynamic_cast<TH1F*>(h[k]->Rebin(rebinfactor,Form("%s%i","h_plot",k)));
			h_plot[k]->Scale(1.0/rebinfactor);
			
			
			if (debug) cout << "\t... finding peaks"<< endl;
			
			/*
			double thresh = ped[k]-15.*pedrms[k];
			pmt_valleys[k] = WaveformAnalysis::valleys(h_plot[k],thresh,1,h_plot[k]->FindBin(S1_maxtime));
			for (int i=0; i<pmt_valleys[k].size(); i++) pmt_valleys_tau[k][i]=hcenter(h_plot[k],pmt_valleys[k].at(i));
			*/
			
			pmt_valleys[k] = WaveformAnalysis::valleys(h_plot[k],1,h_plot[k]->FindBin(S2_maxtime));
			for (int i=0; i<TMath::Min(KMAXNPEAKS,(int)pmt_valleys[k].size()); i++) pmt_valleys_tau[k][i]=hcenter(h_plot[k],pmt_valleys[k].at(i));

			
			if (debug) 
			{ 
				printf("ped[%i] = %0.2f +/- %0.2f\n", k, ped[k], pedrms[k]); 
				printf("ped_end[%i] = %0.2f +/- %0.2f\n", k, ped_end[k], pedrms_end[k]); 
				printf("ch: %d: No of valleys: %lu\n",k,pmt_valleys[k].size());
				for (int i=0; i<TMath::Min(KMAXNPEAKS,(int)pmt_valleys[k].size()); i++) printf("\t%d: %d, %f\n",i,pmt_valleys[k].at(i),pmt_valleys_tau[k][i]);
				
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
		//double q_S1_4us_corr[N_PMT]={0.0};
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
			q_S1_1us[k]   = q*WaveformAnalysis::calc_S1_charge(h[k],ped[k],tau_S1[k]-0.04,tau_S1[k]+1.0);
			q_S1_4us[k]   = q*WaveformAnalysis::calc_S1_charge(h[k],ped[k],tau_S1[k]-0.04,tau_S1[k]+4.0);
			q_S1_80ns[k]  = q*WaveformAnalysis::calc_S1_charge(h[k],ped[k],tau_S1[k]-0.04,tau_S1[k]+0.08);
			q_S1_m2[k]	  = q*WaveformAnalysis::calc_S1_charge_m2(h[k],ped[k],binpeak_S1[k],endbin_S1[k]);
			
			tau_S1_end[k] = hcenter(h[k],endbin_S1[k]);
	
			//printf("S1: q_1us = %f, q_4us = %f, q_m2 = %f, width = %f, amp = %f, tau = %f, tau_end = %f\n", 
			//1e6*q_S1_1us[k],1e6*q_S1_4us[k],1e6*q_S1_m2[k],width_S1[k],amp_S1[k],tau_S1[k],tau_S1_end[k]);
		}

		// S2 calculations
		
		int    binpeak_S2[N_PMT]={0};
		int    binavg_S2[N_PMT]={0};
		double q_S2[N_PMT]={0.0};
		//double q_S2_corr[N_PMT]={0.0};
		//double q_S2_corr_2[N_PMT]={0.0};
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
			
			if (debug) cout << "\t... calculating S2 parameters for all PMTs."<< endl;
			
			for (int k=0; k<N_PMT; k++)	
			{
				double S2_mintime = tau_S1[k]+4.0;
				
				// charge integral calculations - old method 
				
				q_S2[k] 	    = q*WaveformAnalysis::calc_S2_parameters(h[k],ped[k],S2_mintime,S2_maxtime,binavg_S2[k]);
				q_S2_corr_p[k]  = (h_corr_p[k])?q*WaveformAnalysis::integral_S2(h_corr_p[k],S2_mintime,S2_maxtime,ped[k]):0;
				q_S2_corr_m[k]  = (h_corr_m[k])?q*WaveformAnalysis::integral_S2(h_corr_m[k],S2_mintime,S2_maxtime,ped[k]):0;
				
				//binavg_S2[k]    = WaveformAnalysis::calc_S2_binavg(h[k],tau_S1[k]+4.0,S2_maxtime,ped[k],pb); // S2_maxtime was 1048
				tau_S2_avg[k] 	= hcenter(h[k],binavg_S2[k]);
				
				double t_bpc    = WaveformAnalysis::find_S2_peak_coarse(h[k],S2_mintime,S2_maxtime);
				q_S2_m2[k]		= q*WaveformAnalysis::calc_S2_parameters_m2(h[k],h_plot[k],ped[k],pedrms[k],S2_mintime,S2_maxtime,t_bpc,tau_S2_start[k],tau_S2_end[k],width_S2[k]);
				binpeak_S2[k]   = WaveformAnalysis::find_S2_binpeak(h_plot[k],tau_S2_start[k],tau_S2_end[k]);
				tau_S2[k] 		= hcenter(h_plot[k],binpeak_S2[k]);
				amp_S2[k] 		= ADC_to_volts*(ped[k]-hget(h_plot[k],binpeak_S2[k]));
				
				// now for 

				/*
				double q1=0;
				double q2=0;
				q_S2_corr_2[k]  = q*WaveformAnalysis::integral_S2_corr2(h[k],ped[k],tau_S1[k],q1,q2);
				
				q1=q1*q;
				q2=q2*q;
				double vmax = h[k]->GetMaximum()-ped[k];
				q_S2_corr_2[k]+=vmax*200.E-9*ADC_to_volts;
				double q3=vmax*200.E-9*ADC_to_volts;
				*/
			}
		}
		
		if (debug) 
		{
			cout << "\t... done."<< endl << endl;
			
			cout << "Light: (Run,ev): ("<<_runlight <<","<<_nevent <<") /ttotlight: " << totlight[N_PMT] << " " << _nsamples << endl;
			cout << "S1 light 1us: " << q_S1_1us[0]+q_S1_1us[1]+q_S1_1us[2]+q_S1_1us[3]+q_S1_1us[4] << " S2 light: " << q_S2[0]+q_S2[1]+q_S2[2]+q_S2[3]+q_S2[4]<< endl;
		
			//cout <<"Pos1-STAIRS_CRT = "<< _crt_track_pos0[0] << " "<<_crt_track_pos0[1] << " "<<_crt_track_pos0[2]  << endl;
			//cout <<"Pos2-WALL_CRT   = "<< _crt_track_pos1[0] << " "<<_crt_track_pos1[1] << " "<<_crt_track_pos1[2]  << endl;
		}
		
		if (debug) cout << "... preparing light DPD variables ... " << endl;
		
		// prepare DPD variables for filling
		
		_run_light_out=_runlight;
		_nevent_light_out=_nevent;
		_time_light = time_light;
		
		for (int k=0; k<N_PMT; k++) 
		{
			_pmt_wvf_properties[k][0]=ped[k];//ADC Counts
			_pmt_wvf_properties[k][1]=pedrms[k]; //ADC Counts
		
			_pmt_ped[k]=ped[k];
			_pmt_pedrms[k]=pedrms[k];
			_pmt_ped_uncorr[k]=ped_uncorr[k];
			_pmt_pedrms_uncorr[k]=pedrms_uncorr[k];
			_pmt_ped2[k]=ped2[k];
			_pmt_pedrms2[k]=pedrms2[k];
			_pmt_ped_end[k]=ped_end[k];
			_pmt_ped_end_uncorr[k]=ped_end_uncorr[k];
			_pmt_pedrms_end[k]=pedrms_end[k];
			_pmt_pedrms_end_uncorr[k]=pedrms_end_uncorr[k];
			_pmt_wvf_ADC_sat[k]=pmt_wvf_ADC_sat[k];
			_pmt_wvf_end_fit_c0[k]=wvf_end_fit_c0[k];
			_pmt_wvf_end_fit_c1[k]=wvf_end_fit_c1[k];
			_pmt_wvf_end_fit_c2[k]=wvf_end_fit_c2[k];
			_pmt_wvf_end_fit_chi2[k]=wvf_end_fit_chi2[k];
			_pmt_wvf_end_fit_ndof[k]=wvf_end_fit_ndof[k];
			_pmt_wvf_corr_RC_eff[k]=wvf_corr_RC_eff[k];
			_pmt_wvf_max_corr[k]=wvf_max_corr[k];
			_pmt_wvf_max_uncorr[k]=wvf_max_uncorr[k];
			
			_pmt_npeaks[k] = TMath::Min(KMAXNPEAKS,(int)pmt_valleys[k].size());
			//for (int i=0; i<TMath::Min(_pmt_npeaks[k],KMAXNPEAKS); i++) _pmt_peaks_tau[k][i]=pmt_valleys_tau[k][i]; 
			for (int i=0; i<KMAXNPEAKS; i++) _pmt_peaks_tau[k][i]=pmt_valleys_tau[k][i];
	
			_pmt_charge[k]=(double)q*totlight[k];
			_pmt_npe[k]=(double)q*totlight[k]/(charge_e)/gains[k];

			_pmt_S1_charge_1us[k]=q_S1_1us[k]; // in volts
			_pmt_S1_charge_4us[k]=q_S1_4us[k]; // in volts
			_pmt_S1_charge_80ns[k]=q_S1_80ns[k]; // in volts
			_pmt_S1_charge_m2[k]=q_S1_m2[k]; // in volts
			
			_pmt_S1_npe_1us[k]=q_S1_1us[k]/(charge_e)/gains[k];
			_pmt_S1_npe_4us[k]=q_S1_4us[k]/(charge_e)/gains[k];
			_pmt_S1_npe_80ns[k]=q_S1_80ns[k]/(charge_e)/gains[k];
			_pmt_S1_npe_m2[k]=q_S1_m2[k]/(charge_e)/gains[k];
			
			_pmt_S1_width[k]=width_S1[k]; // in ns
			_pmt_S1_amp[k]=amp_S1[k]; // in volts
			_pmt_S1_tau[k]=tau_S1[k]; // in us
			_pmt_S1_tau_end[k]=tau_S1_end[k]; // in us

			_pmt_S2_charge[k]=q_S2[k];
			_pmt_S2_charge_m2[k]=q_S2_m2[k];
			//_pmt_S2_charge_corr[k]=q_S2_corr[k];
			//_pmt_S2_charge_corr_2[k]=q_S2_corr_2[k];
			_pmt_S2_charge_corr_p[k]=q_S2_corr_p[k];
			_pmt_S2_charge_corr_m[k]=q_S2_corr_m[k];
			
			_pmt_S2_npe[k]=q_S2[k]/(charge_e)/gains[k];
			_pmt_S2_npe_m2[k]=q_S2_m2[k]/(charge_e)/gains[k];
			//_pmt_S2_npe_corr[k]=q_S2_corr[k]/(charge_e)/gains[k];
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
		
		_pmt_S2_charge[N_PMT]=0;
		_pmt_S2_npe[N_PMT]=0;
		
		_pmt_S2_charge_m2[N_PMT]=0;
		_pmt_S2_npe_m2[N_PMT]=0;
		
		//_pmt_S2_charge_corr[N_PMT]=0;
		//_pmt_S2_npe_corr[N_PMT]=0;
		
		//_pmt_S2_charge_corr_2[N_PMT]=0;
		//_pmt_S2_npe_corr_2[N_PMT]=0;
		
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

			_pmt_S2_charge[N_PMT]+=_pmt_S2_charge[k]; // in volts
			_pmt_S2_npe[N_PMT]+=_pmt_S2_npe[k];
			
			_pmt_S2_charge_m2[N_PMT]+=_pmt_S2_charge_m2[k]; // in volts
			_pmt_S2_npe_m2[N_PMT]+=_pmt_S2_npe_m2[k];
			
			//_pmt_S2_charge_corr[N_PMT]+=_pmt_S2_charge_corr[k]; // in volts
			//_pmt_S2_npe_corr[N_PMT]+=_pmt_S2_npe_corr[k];
			
			//_pmt_S2_charge_corr_2[N_PMT]+=_pmt_S2_charge_corr_2[k]; // in volts
			
			_pmt_S2_charge_corr_p[N_PMT]+=_pmt_S2_charge_corr_p[k]; // in volts
			_pmt_S2_npe_corr_p[N_PMT]+=_pmt_S2_npe_corr_p[k];
			
			_pmt_S2_charge_corr_m[N_PMT]+=_pmt_S2_charge_corr_m[k]; // in volts
			_pmt_S2_npe_corr_m[N_PMT]+=_pmt_S2_npe_corr_m[k];
		}
		
		
		
		// CRT stuff
	

		if (_crt_daq_match==1 && _crt_reco_sat==1) _crt_matchreco=1;
		else _crt_matchreco=0;
		
		
		// CRT FITTING
		
		float c_yx_m=ERRVAL;
		float c_yx_n=ERRVAL;
		float c_zx_m=ERRVAL;
		float c_zx_n=ERRVAL;
		
		if (_crt_track_param[0]!=ERRVAL && _crt_track_param[1]!=ERRVAL)
		{	
			// convert from CRT coordinate system to LArSoft coordinates
			double TrackIn[3]  = {0.1*_crt_track_door[2],-0.1*_crt_track_door[0],0.1*(_crt_track_door[1]+Z_CENTER)}; // mm to cm
			double TrackOut[3] = {0.1*_crt_track_wall[2],-0.1*_crt_track_wall[0],0.1*(_crt_track_wall[1]+Z_CENTER)}; // mm to cm
			
			// 'converstion' of theta and phi to LArSoft coordinates
			/*
			float dx=(TrackIn[0] - TrackOut[0]);
			float dy=(TrackIn[1] - TrackOut[1]);
			float dz=(TrackIn[2] - TrackOut[2]);
			
			float dyz = sqrt(dy*dy+dz*dz);

			float TanTheta = dx/dyz;
			float TanPhi = -dz/dy;
			//_crt_track_param_LArSoft[0] = TanTheta;
			//_crt_track_param_LArSoft[1] = TanPhi;
			
			printf("_crt_track_door[0] = %f, _crt_track_door[1] = %f, _crt_track_door[2] = %f\n",_crt_track_door[0],_crt_track_door[1],_crt_track_door[2]);
			printf("_crt_track_wall[0] = %f, _crt_track_wall[1] = %f, _crt_track_wall[2] = %f\n",_crt_track_wall[0],_crt_track_wall[1],_crt_track_wall[2]);
			printf("TrackIn[0] = %f, TrackIn[1] = %f, TrackIn[2] = %f\n",TrackIn[0],TrackIn[1],TrackIn[2]);
			printf("TrackOut[0] = %f, TrackOut[1] = %f, TrackOut[2] = %f\n",TrackOut[0],TrackOut[1],TrackOut[2]);
			printf("_crt_track_param[0] = %f, _crt_track_param_LArSoft[0] = %f\n",_crt_track_param[0],_crt_track_param_LArSoft[0]);
			printf("_crt_track_param[1] = %f, _crt_track_param_LArSoft[1] = %f\n",_crt_track_param[1],_crt_track_param_LArSoft[1]);
			*/
			
			
			c_yx_m = (TrackOut[0]-TrackIn[0])/(TrackOut[1]-TrackIn[1]);
			c_yx_n = TrackOut[0]-c_yx_m*TrackOut[1];
			
			c_zx_m = (TrackOut[0]-TrackIn[0])/(TrackOut[2]-TrackIn[2]);
			c_zx_n = TrackOut[0]-c_zx_m*TrackOut[2];
			
			/*
			if (c_yx_n<-3000)
			{
				cout << "WARNING: c_yx_n < -3000!" << endl;
				printf("TrackIn[0] = %f, TrackIn[1] = %f, TrackIn[2] = %f\n",TrackIn[0],TrackIn[1],TrackIn[2]);
				printf("TrackOut[0] = %f, TrackOut[1] = %f, TrackOut[2] = %f\n",TrackOut[0],TrackOut[1],TrackOut[2]);
				printf("c_yx_m = %f, c_yx_n = %f\n",c_yx_m,c_yx_n);	
			}
			*/
			
			if (_crt_matchreco && display_crt)
			{
				cout << "Drawing CRT:" << endl;
				cout << "\tcrtYX m,n: "<< c_yx_m << "  " << c_yx_n << "  " << endl;
				cout << "\tcrtZX m,n: "<< c_zx_m << "  " << c_zx_n << "  " << endl;
			
				TF1* crtYX_line = new TF1("YX_crt_track",Form("%.4f%s%.4f",c_yx_m,"*x+",c_yx_n),-50,50);
				TF1* crtZX_line = new TF1("ZX_crt_track",Form("%.4f%s%.4f",c_zx_m,"*x+",c_zx_n),0,300);
			
				TCanvas* c6 = new TCanvas("c6","c6",1200,600);
				c6->Divide(2,1);
				c6->cd(1);
				crtYX_line->SetLineColor(kGreen+2);
				crtYX_line->Draw();

				c6->cd(2);
				crtZX_line->SetLineColor(kGreen+2);
				crtZX_line->Draw();

				c6->Modified();
				c6->Update();
				lets_pause();
		
				delete c6;
				delete crtYX_line;
				delete crtZX_line;
			}
		}
		
		_crt_fit_yx_m = c_yx_m;
		_crt_fit_yx_n = c_yx_n;
		_crt_fit_zx_m = c_zx_m;
		_crt_fit_zx_n = c_zx_n;
		
		
		if (debug) cout << "... Filling DPD ... " << endl;
		
		dpd->Fill();
		
		/*
		if (_crt_daq_match!=1 || _crt_reco!=1) continue;
		if (_crt_ToF<=0.) continue;
		if (_crt_track_lenFV<=3100.) continue;
		*/
		
		bool myEvent=true;
		for (int k=0; k<N_PMT; k++)
		{
			if (_crt_matchreco!=1 || _crt_isFV!=1 || _crt_ToF<=0.)
				//|| _crt_track_lenFV<=3100) 
				myEvent=false; // double check
			if (_pmt_wvf_ADC_sat[k]!=0) myEvent=false;
			if (_pmt_S1_amp[k]<=40*ADC_to_volts) myEvent=false;
			if (_pmt_S2_amp[k] <= _pmt_pedrms[k]*ADC_to_volts) myEvent=false;
			//if ((_pmt_S2_tau_start[k]-_pmt_S1_tau[k])<=20.) myEvent=false;
		}
		if (myEvent)
		{
			for (int k=0; k<N_PMT; k++)
			{
	   			/*
	   			 ev, 
	   			 Pedestal_ADC[ch], 
	   			 RMS_Ped[ch], 
	   			 Charge_aroundS1_100ns[ch], 
	   			 nPEs_aroundS1_100ns[ch], 
	   			 TMinS2[ch],  
	   			 Rebinned_TMinStartS2[ch], 
	   			 TMin_StartS2[ch], 
	   			 Rebinned_TMinEndS2[ch], 
	   			 TMin_EndS2[ch], 
	   			 Charge_aroundS2[ch], 
	   			 nPEs_S2[ch], 
	   			 Charge_aroundS2_NB[ch], 
	   			 nPEs_S2_NB[ch]
	   			*/
				
	   			cout << _nevent << " " << k << " " << _pmt_ped[k] << " " << _pmt_pedrms[k] << " " << 
	   				_pmt_S1_charge_80ns[k]/q << " " << _pmt_S1_npe_80ns[k] << " " << 
	   				_pmt_S2_tau[k]*1.E3 << " " << _pmt_S2_tau_start[k]*1.E3 << " 0 " << _pmt_S2_tau_end[k]*1.E3 << " 0 " << 
	   				_pmt_S2_charge[k]/q << " " << _pmt_S2_npe[k] << " " << _pmt_S2_charge_m2[k]/q << " " << _pmt_S2_npe_m2[k] << endl;
			}
		}
			
		// Wf displays	
		
		
		//if (_pmt_S2_tau_end[3]<500. && fabs(_pmt_ped[3]-_pmt_ped_end[3])>5) display_waveforms=true;
		
		/*
		for (int k=0; k<N_PMT; k++)
		{
			if (k!=10) continue;
			
			if (_pmt_npeaks[k]>1)
			{
				gStyle->SetOptStat(0);
			    gStyle->SetPadTickX(1);
			    gStyle->SetPadTickY(1);
			
				c1 = new TCanvas("c1","c1");
			
				//h_plot[k]->GetYaxis()->SetTitleOffset(1.25);
				//h_plot[k]->GetYaxis()->SetTitleSize(0.05);
				//h_plot[k]->GetXaxis()->SetTitleSize(0.05);
				
				h_plot[k]->Draw("HIST");
			
				double x = 0.6;
				double y = 0.59;
				TLatex l;
				l.SetNDC();
				l.DrawLatex(x,y+0.06,Form("Run %d, Event %d",_run_charge_out,_nevent_charge_out));
				l.DrawLatex(x,y,Form("No. of peaks: %d",_pmt_npeaks[k]));
				for (int i=0; i<_pmt_npeaks[k]; i++) l.DrawLatex(x+0.05,y-0.06*(i+1),Form("%f",pmt_valleys_tau[k][i]));
	
				c1->Update();
				c1->Modified();
				//c1->Print(Form("/Volumes/data/Dropbox/physics/311/git/Light_analysis/ana_michael/working/S1peaks/run%d_ev%d_ch%d.png",_run_charge_out,_nevent_charge_out,k));	
			
				//lets_pause();
			
			}
		}
		*/
		
		
		if (display_waveforms)		
		//if (_tpc_track_length_trajectory[0]>50. && _tpc_track_max_CBR[0]>0. && _tpc_track_max_CBR[0]<0.2 && display_waveforms)
		{
			TCanvas *c1;
			
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
				//if (k==2 || k==3) h_corr[k]->Draw("hist.same");
			
				double x = 0.15;
				double y = 0.43;
				TLatex l;
				l.SetNDC();
				if (doCharge) l.DrawLatex(x,y+0.06,Form("Run %d, Event %d",_run_charge_out,_nevent_charge_out));
				else l.DrawLatex(x,y+0.06,Form("Run %d, Event %d",_run_light_out,_nevent));
				l.DrawLatex(x,y,Form("ped = %0.2f +/- %0.2f",ped[k],pedrms[k]));
				l.DrawLatex(x,y-0.06,Form("S1: %0.1f ADC, %0.1f PE",_pmt_S1_charge_80ns[k]/q,_pmt_S1_npe_80ns[k]));
				l.DrawLatex(x,y-0.12,Form("S2: %0.1f #mus, %0.1f #mus",_pmt_S2_tau[k],_pmt_S2_tau_avg[k]));
				l.DrawLatex(x,y-0.18,Form("S2_m2: %0.2f #mus, %0.1f #mus, %0.1f #mus",_pmt_S2_tau_start[k],_pmt_S2_tau_end[k],_pmt_S2_width[k]));
				l.DrawLatex(x,y-0.24,Form("S2: %0.1f ADC, %0.1f PE",_pmt_S2_charge[k]/q,_pmt_S2_npe[k]));
				l.DrawLatex(x,y-0.30,Form("S2_m2: %0.1f ADC, %0.1f PE",_pmt_S2_charge_m2[k]/q,_pmt_S2_npe_m2[k]));
				cout << "k = " << k << ", S2 amp = " << _pmt_S2_amp[k] << ", npeaks = " << _pmt_npeaks[k] << endl;
				/*
				l.DrawLatex(x,y-0.06,Form("ped_end = %0.2f +/- %0.2f",ped_end[k],pedrms_end[k]));
				l.DrawLatex(x,y-0.12,Form("S1: %0.2f nC, %0.2f V, %0.1f ns, %0.1f #mus",1e9*_pmt_S1_charge_1us[k],_pmt_S1_amp[k],_pmt_S1_width[k],_pmt_S1_tau[k]));
				if (_pmt_S2_amp[k]>0.) 
				{
					l.DrawLatex(x,y-0.18,Form("S2: %0.2f uC, %0.2f mV, %0.1f #mus",1e6*_pmt_S2_charge[k],_pmt_S2_amp[k]*1000.,_pmt_S2_tau[k]));
					l.DrawLatex(x,y-0.24,Form("S2_m2: %0.2f #mus, %0.1f #mus, %0.1f #mus",_pmt_S2_tau_start[k],_pmt_S2_tau_end[k],_pmt_S2_width[k]));
				}	
				*/
				//l.DrawLatex(x,y-0.30,Form("ped_diff = %0.2f",(ped_end[k]-ped[k])/pedrms[k]));
			}
			c1->Update();
			c1->Modified();
			//c1->Print(Form("chiara/run%d_ev%d.png",_run_light_out,_nevent));	
			
			lets_pause();
			
			//lets_pause();
		}
		
		//if (debug && !display_waveforms) lets_pause();
		
		//if (c1) delete c1;
		
		//debug=false;
		//display_waveforms=false;

		for (int k=0; k<N_PMT; k++) if (h[k]) h[k]->Delete(); 
		for (int k=0; k<N_PMT; k++) if (h_plot[k]) h_plot[k]->Delete(); 
		//for (int k=0; k<N_PMT; k++) if (h_corr[k]) h_corr[k]->Delete(); 
		for (int k=0; k<N_PMT; k++) if (h_corr_p[k]) h_corr_p[k]->Delete(); 
		for (int k=0; k<N_PMT; k++) if (h_corr_m[k]) h_corr_m[k]->Delete(); 

	} //event loop

	dpd->Write();
	ofile.Close();
}


void initDPDvariables()
{
	_run_charge_out=-1;
	_subrun_charge_out=-1;
	_nevent_charge_out=-1;
	_trig_conf=-2;
	_run_light_out=-1;
	_nevent_light_out=-1;

	_time_charge=ERRVAL;
	_time_light=ERRVAL;

	for (int i=0; i<N_PMT; i++)
	{
		_pmt_ped[i]=ERRVAL;
		_pmt_pedrms[i]=ERRVAL;
		_pmt_ped_uncorr[i]=ERRVAL;
		_pmt_pedrms_uncorr[i]=ERRVAL;
		_pmt_ped2[i]=ERRVAL;
		_pmt_pedrms2[i]=ERRVAL;
		_pmt_ped_end[i]=ERRVAL;
		_pmt_ped_end_uncorr[i]=ERRVAL;
		_pmt_pedrms_end[i]=ERRVAL;
		_pmt_pedrms_end_uncorr[i]=ERRVAL;
		_pmt_wvf_ADC_sat[i]=-1;
		_pmt_wvf_end_fit_c0[i]=ERRVAL;
		_pmt_wvf_end_fit_c1[i]=ERRVAL;
		_pmt_wvf_end_fit_c2[i]=ERRVAL;
		_pmt_wvf_end_fit_chi2[i]=ERRVAL;
		_pmt_wvf_end_fit_ndof[i]=-1;
		_pmt_wvf_corr_RC_eff[i]=ERRVAL;
		_pmt_wvf_max_corr[i]=ERRVAL;
		_pmt_wvf_max_uncorr[i]=ERRVAL;
		_pmt_wvf_properties[i][0]=ERRVAL;
		_pmt_wvf_properties[i][1]=ERRVAL;
		_pmt_npeaks[i]=-1;
	
		_pmt_S1_width[i]=ERRVAL;
		_pmt_S1_amp[i]=ERRVAL;
		_pmt_S1_tau[i]=ERRVAL;
		_pmt_S1_tau_end[i]=ERRVAL;

		_pmt_S2_width[i]=ERRVAL;
		_pmt_S2_amp[i]=ERRVAL;
		_pmt_S2_tau[i]=ERRVAL;
		_pmt_S2_tau_start[i]=ERRVAL;
		_pmt_S2_tau_end[i]=ERRVAL;
		_pmt_S2_tau_avg[i]=ERRVAL;
		
		_tpc_drift_time_at_pmt[i]=ERRVAL;
		_tpc_drift_time_at_pmt2[i]=ERRVAL;
		
		for (int p=0; p<KMAXNPEAKS; p++) _pmt_peaks_tau[i][p]=ERRVAL;
	}
	
	for (int i=0; i<N_PMT+1; i++)
	{
		_pmt_S1_charge_1us[i]=ERRVAL;
		_pmt_S1_charge_4us[i]=ERRVAL;
		_pmt_S1_charge_80ns[i]=ERRVAL;
		_pmt_S1_charge_m2[i]=ERRVAL;
		_pmt_S1_npe_1us[i]=ERRVAL;
		_pmt_S1_npe_4us[i]=ERRVAL;
		_pmt_S1_npe_80ns[i]=ERRVAL;
		_pmt_S1_npe_m2[i]=ERRVAL;
	
		_pmt_charge[i]=ERRVAL;
		_pmt_npe[i]=ERRVAL;
	
		_pmt_S2_charge[i]=ERRVAL;
		_pmt_S2_charge_m2[i]=ERRVAL;
		_pmt_S2_charge_corr_p[i]=ERRVAL;
		_pmt_S2_charge_corr_m[i]=ERRVAL;
		_pmt_S2_npe[i]=ERRVAL;
		_pmt_S2_npe_m2[i]=ERRVAL;
		_pmt_S2_npe_corr_p[i]=ERRVAL;
		_pmt_S2_npe_corr_m[i]=ERRVAL;
	}
	
	_tpc_totcharge=ERRVAL;
	_tpc_totrecocharge=ERRVAL;
	_tpc_totcharge_anode=ERRVAL;
	_tpc_totrecocharge_anode=ERRVAL;
	
	for (int v=0; v<2; v++)
	{
		_tpc_totcharge_anode_view[v]=ERRVAL;
		_tpc_totrecocharge_anode_view[v]=ERRVAL;
		_tpc_totcharge_LAr_view[v]=ERRVAL;
		_tpc_totrecocharge_LAr_view[v]=ERRVAL;
	}
	
	_crt_matchreco=-1;
	_crt_fit_yx_m=ERRVAL;
	_crt_fit_yx_n=ERRVAL;
	_crt_fit_zx_m=ERRVAL;
	_crt_fit_zx_n=ERRVAL;

	_ntracks=-1;
	_n_sel_tracks=-1;
	
	for (int i=0; i<KMAXNTRACKS; i++)
	{
		_tpc_track_charge[i]=ERRVAL;
		_tpc_track_charge_anode[i]=ERRVAL;
		
		_tpc_track_start_theta[i]=ERRVAL;
		_tpc_track_start_phi[i]=ERRVAL;
		_tpc_track_end_theta[i]=ERRVAL;
		_tpc_track_end_phi[i]=ERRVAL;

		_tpc_track_start_x[i]=ERRVAL;
		_tpc_track_start_y[i]=ERRVAL;
		_tpc_track_start_z[i]=ERRVAL;
		_tpc_track_end_x[i]=ERRVAL;
		_tpc_track_end_y[i]=ERRVAL;
		_tpc_track_end_z[i]=ERRVAL;

		_tpc_track_fit_yx_m[i]=ERRVAL;
		_tpc_track_fit_yx_n[i]=ERRVAL;
		_tpc_track_fit_yx_chi2[i]=ERRVAL;
		_tpc_track_fit_yx_ndof[i]=ERRVAL;

		_tpc_track_fit_zx_m[i]=ERRVAL;
		_tpc_track_fit_zx_n[i]=ERRVAL;
		_tpc_track_fit_zx_chi2[i]=ERRVAL;
		_tpc_track_fit_zx_ndof[i]=ERRVAL;
		
		_tpc_track_id[i]=ERRVAL;
		
		_tpc_track_momentum[i]=ERRVAL;
		_tpc_track_length_trajectory[i]=ERRVAL;
		_tpc_track_length_straight_line[i]=ERRVAL;

		for (int v=0; v<2; v++)
		{
			_tpc_track_charge_anode_view[i][v]=ERRVAL;
			_tpc_track_charge_LAr_view[i][v]=ERRVAL;
			_tpc_track_dE_view[i][v]=ERRVAL;
			_tpc_track_avg_dEds_view[i][v]=ERRVAL;
			_tpc_track_nhits_view[i][v]=ERRVAL;
			_tpc_track_max_CBR_view[i][v]=ERRVAL;
		}
	}
}

