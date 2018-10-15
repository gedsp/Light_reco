#include "../LIB/WaveformAnalysis.cc"
#include "../config_reco.h"
//#include "THDet.C"

bool debug=false;
bool display_tracks=false;
bool display_waveforms=false;

double ADC_to_volts = 2./4096.;
double charge_e = 1.602E-19;

static const int NMAXPEAKS=100;

// in this library, the S2 will be integrated over the whole waveform.

void make_dpd(int run, int subrun, double gains[N_PMT], string outfilename);

void make_dpd(int run, int subrun, double gains[N_PMT], string outfilename)
{
	TString infilename = Form("%s/%d-%d-Parser.root",matched_data_dir.c_str(),run,subrun);
	
	cout << "Matched data infile = " << infilename.Data() << endl;
	
	TChain * t2 = new TChain("analysistree/anatree");
	t2->Add(Form("%s/%d-%d-Parser.root",matched_data_dir.c_str(),run,subrun));

	//const int kmaxnhits=100000;
	//const int kmaxntracks=2000;
	
	#include "../LIB/branches_2018Feb05.h"
	
	// open pedestal correction file to retrieve correction functions
	TFile *fpc = new TFile(pc_file.c_str());
	static const int nsteps=8;
	TF1* ped_diff_fits[N_PMT][nsteps]={0};
	for (int i=0; i<N_PMT; i++) for (int j=0; j<nsteps; j++) ped_diff_fits[i][j]=(TF1*)fpc->Get(Form("ped_diff_fit_ch%d_%d",i,j));

	double time_charge=0.0;
	double time_light=0.0;
	
	TFile ofile(outfilename.c_str(),"RECREATE");
	if (ofile.IsZombie())
	{
	       cout << "Error opening output file " << outfilename.c_str() << endl;
	       exit(-1);
	}

	TTree *dpd= new TTree("dpd","result tree");

	int _run_out;
	int _subrun_out;
	int _nev_out;

	double _time_charge;
	double _time_light;
	
	double _pmt_charge[N_PMT+1];
	double _pmt_npe[N_PMT+1];
	double _pmt_ped[N_PMT];
	double _pmt_pedrms[N_PMT];
	double _pmt_ped2[N_PMT];
	double _pmt_pedrms2[N_PMT];
	double _pmt_ped_end[N_PMT];
	double _pmt_pedrms_end[N_PMT];
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
	double _pmt_S2_charge_corr[N_PMT+1];
	double _pmt_S2_npe[N_PMT+1];
	double _pmt_S2_npe_m2[N_PMT+1];
	double _pmt_S2_npe_corr[N_PMT+1];
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
	double _tpc_track_mYX[kmaxntracks]={0};
	double _tpc_track_mZX[kmaxntracks]={0};
	double _tpc_track_nYX[kmaxntracks]={0};
	double _tpc_track_nZX[kmaxntracks]={0};
	double _tpc_track_startX[kmaxntracks];
	double _tpc_track_endX[kmaxntracks];
	
	double _tpc_track_start_theta[kmaxntracks]={0};
	double _tpc_track_start_phi[kmaxntracks]={0};
	
	int _tpc_track_fitresult_yx[kmaxntracks]={-1};
	int _tpc_track_fitresult_zx[kmaxntracks]={-1};
	int _crt_matchreco;
	double _tpc_track_length_n[kmaxntracks]={0};

	double _tpc_drift_time_at_pmt_pos[N_PMT]={0};

	TBranch * _bn_ev			= dpd->Branch("ev"      , &_nev_out         , "ev/I"      );
	TBranch * _bn_run			= dpd->Branch("run"      , &_run_out         , "run/I"      );
	TBranch * _bn_subrun		= dpd->Branch("subrun"      , &_subrun_out         , "subrun/I"      );
	TBranch * _bn_nsamples 		= dpd->Branch("nsamples", &_nsamples , "nsamples/I");
	TBranch * _bn_time_charge 	= dpd->Branch("time_charge", &_time_charge  , "time_charge/D");
	TBranch * _bn_time_light 	= dpd->Branch("time_light", &_time_light  , "time_light/D");
	
	TBranch * _bn_ntracks			= dpd->Branch("ntracks"   	,&_ntracks      	, "ntracks/S"   );
	TBranch * _bn_n_sel_tracks		= dpd->Branch("n_sel_tracks"   	,&_n_sel_tracks      	, "n_sel_tracks/S"   );
	TBranch * _bn_tpc_totcharge 	= dpd->Branch("tpc_totcharge"   , &_tpc_totcharge       , "tpc_totcharge/D"    );	
	TBranch * _bn_tpc_totrecocharge = dpd->Branch("tpc_totrecocharge"   , &_tpc_totrecocharge       , "tpc_totrecocharge/D"    );	
	
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
	TBranch * _bn_pmt_pedrms_end   	= dpd->Branch("pmt_pedrms_end"   ,  _pmt_pedrms_end      , "pmt_pedrms_end[5]/D"   );
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
	TBranch * _bn_pmt_S1_npe_1us	= dpd->Branch("pmt_S1_npe_1us"   , _pmt_S1_npe_1us      , "pmt_S1_npe_1us[6]/D"   );
	TBranch * _bn_pmt_S1_npe_4us	= dpd->Branch("pmt_S1_npe_4us"   , _pmt_S1_npe_4us      , "pmt_S1_npe_4us[6]/D"   );
	TBranch * _bn_pmt_S1_npe_80ns	= dpd->Branch("pmt_S1_npe_80ns"   , _pmt_S1_npe_80ns      , "pmt_S1_npe_80ns[6]/D"   );
	TBranch * _bn_pmt_S1_npe_m2		= dpd->Branch("pmt_S1_npe_m2"   , _pmt_S1_npe_m2      , "pmt_S1_npe_m2[6]/D"   );
	TBranch * _bn_pmt_S1_width		= dpd->Branch("pmt_S1_width"   , _pmt_S1_width      , "pmt_S1_width[5]/D"   );
	TBranch * _bn_pmt_S1_amp		= dpd->Branch("pmt_S1_amp"   , _pmt_S1_amp      , "pmt_S1_amp[5]/D"   );
	TBranch * _bn_pmt_S1_tau		= dpd->Branch("pmt_S1_tau"   , _pmt_S1_tau      , "pmt_S1_tau[5]/D"   );
	TBranch * _bn_pmt_S1_tau_end	= dpd->Branch("pmt_S1_tau_end"   , _pmt_S1_tau_end      , "pmt_S1_tau_end[5]/D"   );

	TBranch * _bn_pmt_S2_charge		= dpd->Branch("pmt_S2_charge"   ,_pmt_S2_charge      , "pmt_S2_charge[6]/D"   );
	TBranch * _bn_pmt_S2_charge_m2	= dpd->Branch("pmt_S2_charge_m2"   ,_pmt_S2_charge_m2      , "pmt_S2_charge_m2[6]/D"   );
	TBranch * _bn_pmt_S2_charge_corr = dpd->Branch("pmt_S2_charge_corr"   ,_pmt_S2_charge_corr      , "pmt_S2_charge_corr[6]/D"   );
	TBranch * _bn_pmt_S2_npe		= dpd->Branch("pmt_S2_npe"   , _pmt_S2_npe      , "pmt_S2_npe[6]/D"   );
	TBranch * _bn_pmt_S2_npe_m2		= dpd->Branch("pmt_S2_npe_m2"   , _pmt_S2_npe_m2      , "pmt_S2_npe_m2[6]/D"   );
	TBranch * _bn_pmt_S2_npe_corr	= dpd->Branch("pmt_S2_npe_corr"   , _pmt_S2_npe_corr      , "pmt_S2_npe_corr[6]/D"   );
	TBranch * _bn_pmt_S2_width		= dpd->Branch("pmt_S2_width"   , _pmt_S2_width      , "pmt_S2_width[5]/D"   );
	TBranch * _bn_pmt_S2_amp		= dpd->Branch("pmt_S2_amp"   , _pmt_S2_amp      , "pmt_S2_amp[5]/D"   );
	TBranch * _bn_pmt_S2_tau		= dpd->Branch("pmt_S2_tau"   , _pmt_S2_tau      , "pmt_S2_tau[5]/D"   );
	TBranch * _bn_pmt_S2_tau_avg 	= dpd->Branch("pmt_S2_tau_avg"   , _pmt_S2_tau_avg      , "pmt_S2_tau_avg[5]/D"   );
	TBranch * _bn_pmt_S2_tau_start	= dpd->Branch("pmt_S2_tau_start"   , _pmt_S2_tau_start      , "pmt_S2_tau_start[5]/D"   );
	TBranch * _bn_pmt_S2_tau_end	= dpd->Branch("pmt_S2_tau_end"   , _pmt_S2_tau_end      , "pmt_S2_tau_end[5]/D"   );
	
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


	if (debug) cout << "Number of events: \t" << t2->GetEntries()  << endl;

	for(int ev=0; ev < t2->GetEntries() ; ++ev)
	{		
		cout << "dpdMaker: Event = " << ev << endl;
		
		if (debug) cout << "lets check event "<< ev << " - ("<< _runcharge<<"," <<_subruncharge<< ","<< _event <<")" << endl;

		if(ev%100==0) cout << ev << " over " <<t2->GetEntries()<<" events!" << endl;

 		t2->GetEntry(ev);
		
		
		// calculate delta t
		time_charge=(double)_eventtime_seconds-35.0+(double)_eventtime_nanoseconds*1e-9;
		time_light=(double)_PCTimeTag[0]+(double)_PCTimeTag[2]*1e-3;
		
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
	
		TGraph *trackYX[kmaxntracks], *trackZX[kmaxntracks];
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

		short sorted_tracks[kmaxntracks]; // tracks sorted by total charge in the track

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

		double t_yx_m[kmaxntracks], t_zx_m[kmaxntracks], t_yx_n[kmaxntracks], t_zx_n[kmaxntracks];
		int fail_yx[kmaxntracks], fail_zx[kmaxntracks];
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


		//
		// Light stuff
		//
		
 		TH1F *h[N_PMT];
		TH1F *h_plot[N_PMT];
		double totlight[N_PMT+1]={0};
		
		if (debug) cout << "\tLight variables: Nsamples " << _nsamples  << "\t" << _time_sample<< endl;

		for (int k=0; k<N_PMT; k++) 
		{
			h[k] = new TH1F(Form("%s%i","h",k),Form("%s%i%s","Channel ",k," ;Time [#mus]; Voltage [ADC counts]"),_nsamples,0.0,0.001*_nsamples*_time_sample);

			for (int j=0; j<_nsamples; j++) 
			{
				h[k]->SetBinContent(j+1, _adc_value[k][j]); 
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
		double wvf_end_fit_c0[N_PMT]={0};
		double wvf_end_fit_c1[N_PMT]={0};
		double wvf_end_fit_c2[N_PMT]={0};
		double wvf_end_fit_chi2[N_PMT]={0};
		int    wvf_end_fit_ndof[N_PMT]={0};
		double pedrms_end[N_PMT]={0};
		
		double ped_start = _nsamples==1000?0.0:218.;
		double ped_stop  = _nsamples==1000?0.5:228.;
		double ped2_start = _nsamples==1000?0.0:228.5;
		double ped2_stop  = _nsamples==1000?0.5:229.;
		double ped_end_start = _nsamples==1000?3.5:910.;
		double ped_end_stop  = _nsamples==1000?4.0:920.;
		
		int pedminbin = 10;		
		int pedmaxbin =_nsamples==1000?100:1000;
		
		
		int rebinfactor =_nsamples==1000?1:64;
		double S1_mintime = _nsamples==1000?0.3:228;
		double S1_maxtime = _nsamples==1000?0.65:231;	
		
		vector<int> pmt_valleys[N_PMT];
		double pmt_valleys_tau[N_PMT][NMAXPEAKS]={0};
		
		for (int k=0; k<N_PMT; k++)	
		{
			if (debug) cout << "\t... calculating pedestals"<< endl;
			
			ped[k] = WaveformAnalysis::baseline(h[k],pedrms[k],h[k]->FindBin(ped_start),h[k]->FindBin(ped_stop)-1); 
			
			ped2[k] = WaveformAnalysis::baseline(h[k],pedrms2[k],h[k]->FindBin(ped2_start),h[k]->FindBin(ped2_stop)-1); 
			
			ped_end[k] = WaveformAnalysis::baseline(h[k],pedrms_end[k],h[k]->FindBin(ped_end_start),h[k]->FindBin(ped_end_stop)-1); 
						
			
			if (debug) cout << "\t... rebinning wf histogram for plotting and fitting"<< endl;
			
			h_plot[k] = dynamic_cast<TH1F*>(h[k]->Rebin(rebinfactor,Form("%s%i","h",k)));
			h_plot[k]->Scale(1.0/rebinfactor);
			
			if (debug) cout << "\t... finding peaks"<< endl;
			
			/*
			double thresh = ped[k]-15.*pedrms[k];
			pmt_valleys[k] = WaveformAnalysis::valleys(h_plot[k],thresh,1,h_plot[k]->FindBin(S1_maxtime));
			for (int i=0; i<pmt_valleys[k].size(); i++) pmt_valleys_tau[k][i]=hcenter(h_plot[k],pmt_valleys[k].at(i));
			*/
			
			pmt_valleys[k] = WaveformAnalysis::valleys(h_plot[k],1,h_plot[k]->FindBin(900.));
			for (int i=0; i<pmt_valleys[k].size(); i++) pmt_valleys_tau[k][i]=hcenter(h_plot[k],pmt_valleys[k].at(i));
			
			if (debug) cout << "\t... fitting end of wf"<< endl;
			
			if (_nsamples>1000 && (k==2 || k==3)) // only for positive-base PMTs
			{
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
				
				TF1* fit = new TF1("fit","[0]+[1]*exp((900-x)/[2])",850,1100);
				fit->SetParameter(0,ped[k]);
				fit->SetParameter(1,ped_end[k]-ped[k]);
				fit->SetParameter(2,250.);
				fit->SetParLimits(0,3700,4100);
				fit->SetParLimits(1,0,200);
				fit->SetParLimits(2,0,2000);
				h_plot[k]->Fit("fit","QRMN","",900,1040);
				wvf_end_fit_c0[k] = fit->GetParameter(0);
				wvf_end_fit_c1[k] = fit->GetParameter(1);
				wvf_end_fit_c2[k] = fit->GetParameter(2);
				wvf_end_fit_chi2[k] = fit->GetChisquare();
				wvf_end_fit_ndof[k] = fit->GetNDF();
				
				/*
				TCanvas *c1 = new TCanvas("c1");
				
				h_plot[k]->GetXaxis()->SetRange(h_plot[k]->FindBin(650),h_plot[k]->FindBin(1048));
				h_plot[k]->Draw("hist");
				fit->SetLineColor(kRed);
				fit->Draw("same");
				c1->Modified();
				c1->Update();
				lets_pause();
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
				for (int i=0; i<pmt_valleys[k].size(); i++) printf("\t%d: %d, %f\n",i,pmt_valleys[k].at(i),pmt_valleys_tau[k][i]);
				
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
			tau_S1_end[k] = hcenter(h[k],endbin_S1[k]);
	
			//printf("S1: q_1us = %f, q_4us = %f, q_m2 = %f, width = %f, amp = %f, tau = %f, tau_end = %f\n", 
			//1e6*q_S1_1us[k],1e6*q_S1_4us[k],1e6*q_S1_m2[k],width_S1[k],amp_S1[k],tau_S1[k],tau_S1_end[k]);
		}

		// S2 calculations
		
		int    binpeak_S2[N_PMT]={0};
		int    binavg_S2[N_PMT]={0};
		double q_S2[N_PMT]={0.0};
		double q_S2_corr[N_PMT]={0.0};
		double q_S2_m2[N_PMT]={0.0};
		double amp_S2[N_PMT]={0.0};
		double width_S2[N_PMT]={0};
		double tau_S2[N_PMT]={0}; 
		double tau_S2_start[N_PMT]={0}; 
		double tau_S2_end[N_PMT]={0}; 
		double tau_S2_avg[N_PMT]={0};
		
		// only do calculation if sufficient nsamples
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
				
				
				int p;
				if (tau_S2_avg[k]<400) p=0;
				else if (tau_S2_avg[k]>800) p=nsteps-1;
				else p=(int)((tau_S2_avg[k]-400.)/50.);
				TF1* myfit = ped_diff_fits[k][p];
				
				q_S2_corr[k]    = q*WaveformAnalysis::calc_S2_corrected_charge(h[k],myfit,ped[k],q,tau_S1[k]);
				//printf("k=%d, q_S2 = %f, q_S2_corr = %f\n",k,1.E6*q_S2[k],1.E6*q_S2_corr[k]);
			}
		}
		
	
		if (debug) 
		{
			cout << "\t... done."<< endl << endl;
			
			cout << "Light: (Run,ev): ("<<_run <<","<<_nevent <<") /ttotlight: " << totlight[N_PMT] << " " << _nsamples << endl;
			cout << "S1 light 1us: " << q_S1_1us[0]+q_S1_1us[1]+q_S1_1us[2]+q_S1_1us[3]+q_S1_1us[4] << " S2 light: " << q_S2[0]+q_S2[1]+q_S2[2]+q_S2[3]+q_S2[4]<< endl;
			cout << "Time: charge: " << time_charge << ", light: " << time_light << ", delta: " << time_charge - time_light << endl;
			cout << "Charge (Run,SubRun,ev): (" <<_runcharge << ","<< _subruncharge<<","<<_event<<") /t totcharge: "<<  totq <<endl; 
			cout << "   Length: (";
			for (Int_t i=0;i<_NumberOfTracks_pmtrack;i++)  cout << _Track_Length_pmtrack[i] << ", ";
			cout << ")"<< endl;
			cout << "   Charge: (";
			for (Int_t i=0;i<_NumberOfTracks_pmtrack;i++)cout << trackcharge[i] << ", ";
			cout << ")"<< endl;

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
		//view[0]->Draw("colz");
		//lets_pause();
		//ana_crt_evt(runnum_C,subrun_C,_event,1, gains);
		//view[1]->Draw("colz");
		//lets_pause();
		
		//if (1)//totrecoq/totq > 0.3 && _NumberOfTracks_pmtrack==1)// && _crt_ToF<0) _crt_daq_match==1&&_crt_reco==1 
			
		_run_out=_runcharge;
		_subrun_out=_subruncharge;
		_nev_out=_event;
		_ntracks=_NumberOfTracks_pmtrack;
		_n_sel_tracks=NumberOfSelectedTracks;
		//cout << "_ntracks" << _ntracks << endl;

		_time_charge = time_charge;
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
			_pmt_pedrms_end[k]=pedrms_end[k];
			_pmt_wvf_end_fit_c0[k]=wvf_end_fit_c0[k];
			_pmt_wvf_end_fit_c1[k]=wvf_end_fit_c1[k];
			_pmt_wvf_end_fit_c2[k]=wvf_end_fit_c2[k];
			_pmt_wvf_end_fit_chi2[k]=wvf_end_fit_chi2[k];
			_pmt_wvf_end_fit_ndof[k]=wvf_end_fit_ndof[k];
			
			_pmt_npeaks[k] = (pmt_valleys[k]).size();
			for (int i=0; i<TMath::Min(_pmt_npeaks[k],NMAXPEAKS); i++) _pmt_peaks_tau[k][i]=pmt_valleys_tau[k][i]; 
	
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
			_pmt_S2_charge_corr[k]=q_S2_corr[k];
			_pmt_S2_npe[k]=q_S2[k]/(charge_e)/gains[k];
			_pmt_S2_npe_m2[k]=q_S2_m2[k]/(charge_e)/gains[k];
			_pmt_S2_npe_corr[k]=q_S2_corr[k]/(charge_e)/gains[k];
			_pmt_S2_width[k]=width_S2[k]; // in us
			_pmt_S2_amp[k]=amp_S2[k]; // in volts
			_pmt_S2_tau[k]=tau_S2[k]; // in us
			_pmt_S2_tau_start[k]=tau_S2_start[k]; // in us
			_pmt_S2_tau_end[k]=tau_S2_end[k]; // in us
			_pmt_S2_tau_avg[k]=tau_S2_avg[k]; // in us
			
			_tpc_drift_time_at_pmt_pos[k]=tpc_drift_time_at_pmt_pos[k];
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
		
		_pmt_S2_charge_corr[N_PMT]=0;
		_pmt_S2_npe_corr[N_PMT]=0;
		
		
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
			
			_pmt_S2_charge_corr[N_PMT]+=_pmt_S2_charge_corr[k]; // in volts
			_pmt_S2_npe_corr[N_PMT]+=_pmt_S2_npe_corr[k];
		}

		_crt_mYX = crtYX_m;
		_crt_mZX = crtZX_m;
		_crt_nYX = crtYX_n;
		_crt_nZX = crtZX_n;
		_crt_ToF_n = _crt_ToF;
		_crt_track_lenFV_n = _crt_track_lenFV;
		for (int k=0; k<3; k++) _crt_pmt_dist_n[k] = _crt_pmt_dist[k];
		_crt_isFV_n = _crt_isFV;
		
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

		//Dor (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "TL: " <<  _tpc_track_nZX[j] << " " << t_zx_n[sorted_tracks[j]] << " " <<t_zx_n[j] << endl;

		//Dor (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "TL: " <<  sorted_tracks[j] << " " << _tpc_track_length_pmtrack[sorted_tracks[j]] << " " <<_tpc_track_length_pmtrack[j] << endl;


//			std::copy(std::begin(trackcharge), std::end(trackcharge), _tpc_track_charge);
//			std::copy(std::begin(t_yx_m), std::end(t_yx_m), _tpc_track_mYX);
//			std::copy(std::begin(t_zx_m), std::end(t_zx_m), _tpc_track_mZX);
//			std::copy(std::begin(t_yx_n), std::end(t_yx_n), _tpc_track_nYX);
//			std::copy(std::begin(t_zx_n), std::end(t_zx_n), _tpc_track_nZX);
//			std::copy(std::begin(fail_yx), std::end(fail_yx), _tpc_track_fitresult_yx);
//			std::copy(std::begin(fail_zx), std::end(fail_zx), _tpc_track_fitresult_zx);

//			_tpc_track_mYX = t_yx_m[0];
//			_tpc_track_mZX = t_zx_n[0];
//			_tpc_track_nYX = t_yx_m[0];
//			_tpc_track_nZX = t_zx_n[0];

		dpd->Fill();
		
	
		TCanvas *c1;
		
		// Wf displays	
		
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
				l.DrawLatex(x,y+0.06,Form("Run %d, Event %d",_run_out,ev));
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
		}

		if (debug || display_waveforms) lets_pause();
		
		//if (c1) delete c1;
		
		//debug=false;
		//display_waveforms=false;

		for (int k=0; k<N_PMT; k++) if (h[k]) h[k]->Delete(); 
		for (int k=0; k<N_PMT; k++) if (h_plot[k]) h_plot[k]->Delete(); 

		if (view[0]) view[0]->Delete();
		if (view[1]) view[1]->Delete();
		
		for (int j=0; j<_NumberOfTracks_pmtrack; j++) if (trackYX[j]) trackYX[j]->Delete();
		for (int j=0; j<_NumberOfTracks_pmtrack; j++) if (trackZX[j]) trackZX[j]->Delete();

	} //event loop

	dpd->Write();
	ofile.Close();
}

