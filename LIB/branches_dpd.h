
	TBranch * _bn_ev;	int _ev_out;		dpd->SetBranchAddress("ev"   	,&_ev_out      	, &_bn_ev   );
	TBranch * _bn_run;	int _run_out;		dpd->SetBranchAddress("run"   	,&_run_out      	, &_bn_run   );
	TBranch * _bn_subrun;	int _subrun_out;		dpd->SetBranchAddress("subrun"   	,&_subrun_out      	, &_bn_subrun   );
	TBranch * _bn_ntracks;	short _ntracks;		dpd->SetBranchAddress("ntracks"   	,&_ntracks      	, &_bn_ntracks   );
	TBranch * _bn_nsamples;  int _nsamples;     dpd->SetBranchAddress("nsamples"    ,&_nsamples       , &_bn_nsamples   );

	TBranch * _bn_time_charge; float _time_charge;    dpd->SetBranchAddress("time_charge",    &_time_charge   , &_bn_time_charge  );
	TBranch * _bn_time_light; float _time_light;      dpd->SetBranchAddress("time_light",    &_time_light   , &_bn_time_light  );

	TBranch * _bn_tpc_totcharge; float _tpc_totcharge; dpd->SetBranchAddress("tpc_totcharge"   , &_tpc_totcharge       , &_bn_tpc_totcharge    );
	TBranch * _bn_tpc_totrecocharge; float _tpc_totrecocharge; dpd->SetBranchAddress("tpc_totrecocharge"   , &_tpc_totcharge       , &_bn_tpc_totcharge    );	
	TBranch * _bn_tpc_track_charge; float _tpc_track_charge[500];	dpd->SetBranchAddress("tpc_track_charge",_tpc_track_charge     , &_bn_tpc_track_charge   );


	TBranch * _bn_crt_matchreco;	int _crt_matchreco;		dpd->SetBranchAddress("crt_matchreco"   	,&_crt_matchreco      	, &_bn_crt_matchreco   );

	TBranch * _bn_crt_mYX;	float _crt_mYX;	dpd->SetBranchAddress("crt_mYX"   	,&_crt_mYX      	, &_bn_crt_mYX   );
	TBranch * _bn_crt_mZX;	float _crt_mZX;	dpd->SetBranchAddress("crt_mZX"   	,&_crt_mZX      	, &_bn_crt_mZX  );
	TBranch * _bn_crt_nYX;	float _crt_nYX;	dpd->SetBranchAddress("crt_nYX"   	,&_crt_nYX      	, &_bn_crt_nYX   );
	TBranch * _bn_crt_nZX;	float _crt_nZX;	dpd->SetBranchAddress("crt_nZX"   	,&_crt_nZX      	, &_bn_crt_nZX   );

	TBranch * _bn_tpc_mYX;	float _tpc_mYX[500];	dpd->SetBranchAddress("tpc_mYX"   	,_tpc_mYX      	, &_bn_tpc_mYX   );
	TBranch * _bn_tpc_mZX;	float _tpc_mZX[500];	dpd->SetBranchAddress("tpc_mZX"   	,_tpc_mZX      	, &_bn_tpc_mZX   );
	TBranch * _bn_tpc_nYX;	float _tpc_nYX[500];	dpd->SetBranchAddress("tpc_nYX"   	,_tpc_nYX      	, &_bn_tpc_nYX   );
	TBranch * _bn_tpc_nZX;	float _tpc_nZX[500];	dpd->SetBranchAddress("tpc_nZX"   	,_tpc_nZX      	, &_bn_tpc_nZX   );

	TBranch * _bn_tpc_startX;float _tpc_startX[500];	dpd->SetBranchAddress("tpc_startX"   	,_tpc_startX      	, &_bn_tpc_startX   );
	TBranch * _bn_tpc_endX;	float _tpc_endX[500];	dpd->SetBranchAddress("tpc_endX"   	,_tpc_endX      	, &_bn_tpc_endX   );

	float _Track_Length_pmtrack[500];
	TBranch * _bn_Track_Length; dpd->SetBranchAddress("Track_Length"		, _Track_Length_pmtrack		,  &_bn_Track_Length      );

	TBranch * _bn_tpc_drift_time_at_pmt_pos; float _tpc_drift_time_at_pmt_pos[5]; dpd->SetBranchAddress("tpc_drift_time_at_pmt_pos", _tpc_drift_time_at_pmt_pos		, &_bn_tpc_drift_time_at_pmt_pos       );


	TBranch * _bn_pmt_charge;	float _pmt_charge[6];	dpd->SetBranchAddress("pmt_charge"   , _pmt_charge   , &_bn_pmt_charge );
	TBranch * _bn_pmt_npe; 		float _pmt_npe[6]; 	dpd->SetBranchAddress("pmt_npe"	, _pmt_npe      , &_bn_pmt_npe   	);
	TBranch * _bn_pmt_wvf_prop; 	float _pmt_wvf_properties[5][2]; 	dpd->SetBranchAddress("pmt_wvf_properties"	, _pmt_wvf_properties      , &_bn_pmt_wvf_prop 	);

    TBranch * _bn_pmt_ped; 			float _pmt_ped[5];     	dpd->SetBranchAddress("pmt_ped", _pmt_ped, &_bn_pmt_ped);
	TBranch * _bn_pmt_pedrms; 		float _pmt_pedrms[5];  	dpd->SetBranchAddress("pmt_pedrms", _pmt_pedrms, &_bn_pmt_pedrms);
    TBranch * _bn_pmt_ped_end; 		float _pmt_ped[5];     	dpd->SetBranchAddress("pmt_ped", _pmt_ped, &_bn_pmt_ped);
	TBranch * _bn_pmt_pedrms_end; 	float _pmt_pedrms_end[5];  dpd->SetBranchAddress("pmt_pedrms_end", _pmt_pedrms_end, &_bn_pmt_pedrms_end);

	TBranch * _bn_pmt_S1_charge;	float _pmt_S1_charge[6];dpd->SetBranchAddress("pmt_S1_charge"   , _pmt_S1_charge   , &_bn_pmt_S1_charge );
	TBranch * _bn_pmt_S1_npe; 		float _pmt_S1_npe[6]; 	dpd->SetBranchAddress("pmt_S1_npe"	, _pmt_S1_npe      , &_bn_pmt_S1_npe   	);
	TBranch * _bn_pmt_S1_width; 	float _pmt_S1_width[5]; dpd->SetBranchAddress("pmt_S1_width"	, _pmt_S1_width    , &_bn_pmt_S1_width  );
	TBranch * _bn_pmt_S1_amp; 		float _pmt_S1_amp[5]; 	dpd->SetBranchAddress("pmt_S1_amp"	, _pmt_S1_amp      , &_bn_pmt_S1_amp   	);
	TBranch * _bn_pmt_S1_tau; 		float _pmt_S1_tau[5]; 	dpd->SetBranchAddress("pmt_S1_tau"	, _pmt_S1_tau      , &_bn_pmt_S1_tau   	);

	TBranch * _bn_pmt_S2_charge;		float _pmt_S2_charge[6];	dpd->SetBranchAddress("pmt_S2_charge"   , _pmt_S2_charge   , &_bn_pmt_S2_charge );
	TBranch * _bn_pmt_S2_npe; 			float _pmt_S2_npe[6]; 		dpd->SetBranchAddress("pmt_S2_npe"	, _pmt_S2_npe      , &_bn_pmt_S2_npe   	);
	TBranch * _bn_pmt_S2_width; 		float _pmt_S2_width[5]; 	dpd->SetBranchAddress("pmt_S2_width"	, _pmt_S2_width    , &_bn_pmt_S2_width  );
	TBranch * _bn_pmt_S2_amp; 			float _pmt_S2_amp[5]; 		dpd->SetBranchAddress("pmt_S2_amp"	, _pmt_S2_amp      , &_bn_pmt_S2_amp   	);
	TBranch * _bn_pmt_S2_tau; 			float _pmt_S2_tau[5]; 		dpd->SetBranchAddress("pmt_S2_tau"	, _pmt_S2_tau      , &_bn_pmt_S2_tau   	);
	
	TBranch * _bn_pmt_S2_gaus_charge;	float _pmt_S2_gaus_charge[6];dpd->SetBranchAddress("pmt_S2_gaus_charge"   , _pmt_S2_gaus_charge   , &_bn_pmt_S2_gaus_charge );
	TBranch * _bn_pmt_S2_gaus_npe; 		float _pmt_S2_gaus_npe[6]; 	dpd->SetBranchAddress("pmt_S2_gaus_npe"	, _pmt_S2_gaus_npe      , &_bn_pmt_S2_gaus_npe   	);
	TBranch * _bn_pmt_S2_gaus_width; 	float _pmt_S2_gaus_width[5]; dpd->SetBranchAddress("pmt_S2_gaus_width"	, _pmt_S2_gaus_width    , &_bn_pmt_S2_gaus_width  );
	TBranch * _bn_pmt_S2_gaus_amp; 		float _pmt_S2_gaus_amp[5]; 	dpd->SetBranchAddress("pmt_S2_gaus_amp"	, _pmt_S2_gaus_amp      , &_bn_pmt_S2_gaus_amp   	);
	TBranch * _bn_pmt_S2_gaus_tau; 		float _pmt_S2_gaus_tau[5]; 	dpd->SetBranchAddress("pmt_S2_gaus_tau"	, _pmt_S2_gaus_tau      , &_bn_pmt_S2_gaus_tau   	);

