#include "../LIB/WaveUtils.h"
#include "../config_reco.h"
//#include "THDet.C"

static const int NSAMP_CUTOFF = 10000;

bool debug=false;
bool display_tracks=false;
bool display_waveforms=false;

// in this library, the S2 will be integrated over the whole waveform.

void make_dpd(int run, int subrun, float gains[5], string outfilename);

void make_dpd(int run, int subrun, float gains[5], string outfilename)
{
	TString infilename = Form("%s/%d-%d-Parser.root",matched_data_dir.c_str(),run,subrun);
	
	cout << "Matched data infile = " << infilename.Data() << endl;
	
	TChain * t2 = new TChain("analysistree/anatree");
	t2->Add(Form("%s/%d-%d-Parser.root",matched_data_dir.c_str(),run,subrun));

	//const int kmaxnhits=100000;
	//const int kmaxntracks=2000;
	
	#include "../LIB/branches_2018Feb05.h"
	
	TH1F * htime = new TH1F("delta_time_per_event","TimeStamp_{Charge} - TimeStamp_{Light}",200,-1,1);
	TH2F * hcorr = new TH2F("mean_deltat","Correlation of totlight vs totcharge;Charge;Light",1000,0,0,1000,0,0);

	float time_charge=0.0;
	float time_light=0.0;

	TFile ofile(outfilename.c_str(),"RECREATE");
	TTree *dpd= new TTree("dpd","result tree");

	int _run_out;
	int _subrun_out;
	int _nev_out;

	float _time_charge;
	float _time_light;
	
	float _pmt_charge[6];
	float _pmt_npe[6];
	float _pmt_ped[5];
	float _pmt_pedrms[5];
	float _pmt_ped_end[5];
	float _pmt_pedrms_end[5];

	float _pmt_S1_charge[6];
	float _pmt_S1_npe[6];
	float _pmt_S1_width[5];
	float _pmt_S1_amp[5];
	float _pmt_S1_tau[5];

	float _pmt_S2_charge[6];
	float _pmt_S2_npe[6];
	float _pmt_S2_width[5];
	float _pmt_S2_amp[5];
	float _pmt_S2_tau[5];
	
	float _pmt_S2_gaus_charge[6];
	float _pmt_S2_gaus_npe[6];
	float _pmt_S2_gaus_width[5];
	float _pmt_S2_gaus_amp[5];
	float _pmt_S2_gaus_tau[5];

	float _tpc_totcharge;
	float _tpc_totrecocharge;
	float _tpc_track_charge[1000]={0};

	float _crt_mYX;
	float _crt_mZX;
	float _crt_nYX;
	float _crt_nZX;

	float _tpc_mYX[kmaxntracks]={0};
	float _tpc_mZX[kmaxntracks]={0};
	float _tpc_nYX[kmaxntracks]={0};
	float _tpc_nZX[kmaxntracks]={0};
	short _ntracks;
	int _tpc_track_fitresult_yx[kmaxntracks]={-1};
	int _tpc_track_fitresult_zx[kmaxntracks]={-1};
	int _crt_matchreco;
	float _Track_Length_n[kmaxntracks]={0};

	float _tpc_drift_time_at_pmt_pos[5]={0};

	TBranch * _bn_ev=dpd->Branch("ev"      , &_nev_out         , "ev/I"      );
	TBranch * _bn_run=dpd->Branch("run"      , &_run_out         , "run/I"      );
	TBranch * _bn_subrun=dpd->Branch("subrun"      , &_subrun_out         , "subrun/I"      );
	TBranch * _bn_nsamples = dpd->Branch("nsamples", &_nsamples , "nsamples/I");
	TBranch * _bn_time_charge = dpd->Branch("time_charge", &_time_charge  , "time_charge/F");
	TBranch * _bn_time_light = dpd->Branch("time_light", &_time_light  , "time_light/F");

	TBranch * _bn_crt_matchreco=dpd->Branch("crt_matchreco"      , &_crt_matchreco         , "crt_matchreco/I"      );
	float _pmt_wvf_properties[5][2];
	TBranch * _bn_pmt_wvf_properties = dpd->Branch("pmt_wvf_properties"    , _pmt_wvf_properties       , "pmt_wvf_properties[5][2]/F"    );

	TBranch * _bn_pmt_charge 	= dpd->Branch("pmt_charge"    , _pmt_charge       , "pmt_charge[6]/F"    );	
	TBranch * _bn_pmt_npe		= dpd->Branch("pmt_npe"   , _pmt_npe      , "pmt_npe[6]/F"   );
	TBranch * _bn_pmt_ped       = dpd->Branch("pmt_ped"   ,  _pmt_ped      , "pmt_ped[5]/F"   );
	TBranch * _bn_pmt_pedrms    = dpd->Branch("pmt_pedrms", _pmt_pedrms   , "pmt_pedrms[5]/F"   );
	TBranch * _bn_pmt_ped_end   = dpd->Branch("pmt_ped_end"   ,  _pmt_ped_end      , "pmt_ped_end[5]/F"   );
	TBranch * _bn_pmt_pedrms_end   = dpd->Branch("pmt_pedrms_end"   ,  _pmt_pedrms_end      , "pmt_pedrms_end[5]/F"   );

	TBranch * _bn_pmt_S1_charge	= dpd->Branch("pmt_S1_charge"   ,_pmt_S1_charge      , "pmt_S1_charge[6]/F"   );
	TBranch * _bn_pmt_S1_npe	= dpd->Branch("pmt_S1_npe"   , _pmt_S1_npe      , "pmt_S1_npe[6]/F"   );
	TBranch * _bn_pmt_S1_width	= dpd->Branch("pmt_S1_width"   , _pmt_S1_width      , "pmt_S1_width[5]/F"   );
	TBranch * _bn_pmt_S1_amp	= dpd->Branch("pmt_S1_amp"   , _pmt_S1_amp      , "pmt_S1_amp[5]/F"   );
	TBranch * _bn_pmt_S1_tau	= dpd->Branch("pmt_S1_tau"   , _pmt_S1_tau      , "pmt_S1_tau[5]/F"   );

	TBranch * _bn_pmt_S2_charge	= dpd->Branch("pmt_S2_charge"   ,_pmt_S2_charge      , "pmt_S2_charge[6]/F"   );
	TBranch * _bn_pmt_S2_npe	= dpd->Branch("pmt_S2_npe"   , _pmt_S2_npe      , "pmt_S2_npe[6]/F"   );
	TBranch * _bn_pmt_S2_width	= dpd->Branch("pmt_S2_width"   , _pmt_S2_width      , "pmt_S2_width[5]/F"   );
	TBranch * _bn_pmt_S2_amp	= dpd->Branch("pmt_S2_amp"   , _pmt_S2_amp      , "pmt_S2_amp[5]/F"   );
	TBranch * _bn_pmt_S2_tau	= dpd->Branch("pmt_S2_tau"   , _pmt_S2_tau      , "pmt_S2_tau[5]/F"   );
	
	TBranch * _bn_pmt_S2_gaus_charge	= dpd->Branch("pmt_S2_gaus_charge"   ,_pmt_S2_gaus_charge      , "pmt_S2_gaus_charge[6]/F"   );
	TBranch * _bn_pmt_S2_gaus_npe	= dpd->Branch("pmt_S2_gaus_npe"   , _pmt_S2_gaus_npe      , "pmt_S2_gaus_npe[6]/F"   );
	TBranch * _bn_pmt_S2_gaus_width	= dpd->Branch("pmt_S2_gaus_width"   , _pmt_S2_gaus_width      , "pmt_S2_gaus_width[5]/F"   );
	TBranch * _bn_pmt_S2_gaus_amp	= dpd->Branch("pmt_S2_gaus_amp"   , _pmt_S2_gaus_amp      , "pmt_S2_gaus_amp[5]/F"   );
	TBranch * _bn_pmt_S2_gaus_tau	= dpd->Branch("pmt_S2_gaus_tau"   , _pmt_S2_gaus_tau      , "pmt_S2_gaus_tau[5]/F"   );


	TBranch * _bn_ntracks		= dpd->Branch("ntracks"   	,&_ntracks      	, "ntracks/S"   );
	TBranch * _bn_tpc_totcharge 	= dpd->Branch("tpc_totcharge"   , &_tpc_totcharge       , "tpc_totcharge/F"    );	
	TBranch * _bn_tpc_totrecocharge 	= dpd->Branch("tpc_totrecocharge"   , &_tpc_totrecocharge       , "tpc_totrecocharge/F"    );	
	TBranch * _bn_tpc_track_charge	= dpd->Branch("tpc_track_charge",&_tpc_track_charge     , "tpc_track_charge[ntracks]/F"   );

	TBranch * _bn_crt_mYX		= dpd->Branch("crt_mYX"   	,&_crt_mYX      	, "crt_mYX/F"   );
	TBranch * _bn_crt_mZX		= dpd->Branch("crt_mZX"   	,&_crt_mZX      	, "crt_mZX/F"   );
	TBranch * _bn_crt_nYX		= dpd->Branch("crt_nYX"   	,&_crt_nYX      	, "crt_nYX/F"   );
	TBranch * _bn_crt_nZX		= dpd->Branch("crt_nZX"   	,&_crt_nZX      	, "crt_nZX/F"   );

	TBranch * _bn_tpc_mYX		= dpd->Branch("tpc_mYX"   	,_tpc_mYX      	, "tpc_mYX[ntracks]/F"   );
	TBranch * _bn_tpc_mZX		= dpd->Branch("tpc_mZX"   	,_tpc_mZX      	, "tpc_mZX[ntracks]/F"   );
	TBranch * _bn_tpc_nYX		= dpd->Branch("tpc_nYX"   	,_tpc_nYX      	, "tpc_nYX[ntracks]/F"   );
	TBranch * _bn_tpc_nZX		= dpd->Branch("tpc_nZX"   	,_tpc_nZX      	, "tpc_nZX[ntracks]/F"   );

	float _tpc_startX[kmaxntracks];
	float _tpc_endX[kmaxntracks];
	TBranch * _bn_tpc_startX	= dpd->Branch("tpc_startX"   	,_tpc_startX      	, "tpc_startX[ntracks]/F"   );
	TBranch * _bn_tpc_endX		= dpd->Branch("tpc_endX"   	,_tpc_endX      	, "tpc_endX[ntracks]/F"   );

	TBranch * _bn_tpc_track_fitresult_yx= dpd->Branch("tpc_track_fitresult_yx"   	,_tpc_track_fitresult_yx      	, "tpc_track_fitresult_yx[ntracks]/I"   );
	TBranch * _bn_tpc_track_fitresult_zx= dpd->Branch("tpc_track_fitresult_zx"   	,_tpc_track_fitresult_zx      	, "tpc_track_fitresult_zx[ntracks]/I"   );

	TBranch * _bn_Track_Length_n	= dpd->Branch("Track_Length"		, _Track_Length_n		, "Track_Length[ntracks]/F"       );
	TBranch * _bn_tpc_drift_time_at_pmt_pos	= dpd->Branch("tpc_drift_time_at_pmt_pos"		, _tpc_drift_time_at_pmt_pos		, "tpc_drift_time_at_pmt_pos[5]/F"       );



	if (debug) cout << "Number of events: \t" << t2->GetEntries()  << endl;
	

	for(int ev=0; ev < t2->GetEntries() ; ++ev)
	{		
		cout << "dpdMaker: Event = " << ev << endl;
		
		if (debug) cout << "lets check event "<< ev << " - ("<< _runcharge<<"," <<_subruncharge<< ","<< _event <<")" << endl;

		if(ev%100==0) cout << ev << " over " <<t2->GetEntries()<<" events!" << endl;

 		t2->GetEntry(ev);

 		TH1F *h[5];
		if (debug) cout << "\tLight variables: Nsamples " << _nsamples  << "\t" << _time_sample<< endl;

		for (int k=0; k<5; k++) h[k] = new TH1F(Form("%s%i","h",k),Form("%s%i%s","Channel ",k," ;Time [#mus]; Voltage [ADC counts]"),_nsamples,0.0,0.001*_nsamples*_time_sample);

		if (debug) cout << "pmt histograms created " << ev << endl;

		//
		// Charge stuff
		//
	
		TH2F * view[2];
		view[0] = new TH2F("view0","View 0;Channel;Time",320,0,320,1667,0,1667);
		view[1] = new TH2F("view1","View 1;Channel;Time",960,0,960,1667,0,1667);

		float pmt_pos[5]={0};

		pmt_pos[0]=150.0-92.246; // pmt position in Y coordinate (larsoft coordinate systems) in cm
		pmt_pos[1]=150.0-46.123;
		pmt_pos[2]=150.0-0;
		pmt_pos[3]=150.0+46.123;
		pmt_pos[4]=150.0+92.246;


		float S2_width_tolerance_channel=5;//cm

		float tpc_drift_time_at_pmt_pos[5]={0};
		int tpc_drift_time_at_pmt_pos_counter[5]={0};
		//if (debug) cout << _no_hits << endl;

		float totq=0;
		float totrecoq=0;
		float trackcharge[1000]={0};
		int number_reco_hits=0;
	
		TGraph *trackYX[kmaxntracks], *trackZX[kmaxntracks];
		for (int j=0; j<_NumberOfTracks_pmtrack; j++) trackYX[j]= new TGraph();
		for (int j=0; j<_NumberOfTracks_pmtrack; j++) trackZX[j]= new TGraph();
	
		//for (int j=0; j<_no_hits; j++)	if (debug) cout << "Ev "<<ev << " hit: " << j <<"/"<<_no_hits<< "  _charge:" <<  _Hit_ChargeIntegral[j] <<  endl;

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

		//for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "\t\tSorted track " << j << " " << sorted_tracks[j] << " " << trackcharge[sorted_tracks[j]] <<"IAC"<< endl;

		//cout << "\t... 3d tracks sorted." << endl;

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
			}
		}

		for (int j=0; j<_no_hits; j++)
		{
			if(_Hit_ChargeIntegral[j]==_Hit_ChargeIntegral[j])
			{
				for (int k=0; k<5; k++) if (_hit_trkid[j]==sorted_tracks[0]&&_Hit_View[j]==1 && (_Hit_Channel[j]-320)*300.0/960.0>pmt_pos[k]-S2_width_tolerance_channel && (_Hit_Channel[j]-320)*300.0/960.0<pmt_pos[k]+S2_width_tolerance_channel ) 
				{
					tpc_drift_time_at_pmt_pos[k]+=0.4*_Hit_PeakTime[j];
					tpc_drift_time_at_pmt_pos_counter[k]++;
					//cout << "channel " <<  _Hit_Channel[j]<< " in cm " << _Hit_Channel[j]*300.0/960.0<< " time " << 0.4*_Hit_PeakTime[j] << " assing to pmt " << k << endl;
				}
			}
		}

		for (int k=0; k<5; k++)
		{
			if (tpc_drift_time_at_pmt_pos_counter[k]!=0) tpc_drift_time_at_pmt_pos[k]/=tpc_drift_time_at_pmt_pos_counter[k];
			else tpc_drift_time_at_pmt_pos[k]=-1;
		}


		// calculate delta t
		// MAL: This doesn't seem to work TODO
		time_charge=(double)_eventtime_seconds-35.0+(double)_eventtime_nanoseconds*1e-9;
		time_light=(double)_PCTimeTag[0]+(double)_PCTimeTag[2]*1e-3;

		if (debug)
		{
			cout << "\tTime: C_" << _eventtime_seconds << " L_" << _PCTimeTag[0] << endl;
			cout << "\tdt: " <<time_charge-time_light << endl<< endl;
		}
		
		float ped[5]={0};
		float pedrms[5]={0};
		float ped_end[5]={0};
		float pedrms_end[5]={0};
		
		double totlight[6]={0};
		//for (int k=0; k<5; k++)for (int j=0; j<_nsamples_2; j++) {cout << k << " " <<  j << " " <<  _adc_value_2[k][j] <<  endl;}

		for (int k=0; k<5; k++) for (int j=0; j<_nsamples; j++) {h[k]->SetBinContent(j+1, _adc_value[k][j]); totlight[k]+= _adc_value[k][j];}
		
		for (int k=0; k<5; k++)	
		{
			WaveUtils wu;
			if (_nsamples<=NSAMP_CUTOFF) 
			{
				ped[k] = wu.FindPedestal(*h[k], 10, 100); 
				pedrms[k]=wu.GetPedestalRMS();
				ped_end[k] = wu.FindPedestal(*h[k],h[k]->GetSize()-100,h[k]->GetSize()-10);
				pedrms_end[k]=wu.GetPedestalRMS();
			}
			else 
			{
				ped[k] = wu.FindPedestal(*h[k], 10, 1000); 
				pedrms[k]=wu.GetPedestalRMS();
				ped_end[k] = wu.FindPedestal(*h[k],h[k]->GetSize()-1000,h[k]->GetSize()-10);
				pedrms_end[k]=wu.GetPedestalRMS();
			}
			
			if (debug) 
			{ 
				printf("ped[%i] = %0.2f +/- %0.2f\n", k, ped[k], pedrms[k]); 
				printf("ped_end[%i] = %0.2f +/- %0.2f\n", k, ped_end[k], pedrms_end[k]); 
			}

			totlight[k]=-totlight[k]+ped[k]*_nsamples;
		}

		for (int k=0; k<5; k++)	{ totlight[5]+=totlight[k];}

		//cout << "\t... integrating S1 for all PMTs."<< endl;
		
		float q;
		q=(2.0/4096.0)*_time_sample*1e-9/(50.); // from iADC to Coulombs
		float s1_width[5]={0.0};
		float q_S1[5]={0.0};
		int binpeak[5]={0};
		float tau_S1[5]={0}; 
		float tau_S2[5]={0}; 
		float tau_S2_gaus[5]={0};
		float tau_S2_error[5]={0};
		float tau_S2_gaus_error[5]={0};
		float amp_S1[5]={0.0};
		for (int k=0; k<5; k++) q_S1[k]=q*get_S1full(h[k],s1_width[k], binpeak[k], ped[k]); // s1 in coulombs
		for (int k=0; k<5; k++) tau_S1[k]=h[k]->GetBinCenter(binpeak[k]); // in microseconds
		for (int k=0; k<5; k++) amp_S1[k]=(ped[k]-h[k]->GetBinContent(binpeak[k]))*2/4096; // in volts


		// S2 calculation
		
		float q_S2[5]={0.0};
		float q_S2_error[5]={0.0};
		float amp_S2[5]={0.0};
		float s2_width[5]={0};
				
		float q_S2_gaus[5]={0.0};
		float q_S2_gaus_error[5]={0.0};
		float amp_S2_gaus[5]={0.0};
		float S2_gaus_width[5]={0};
		
		TH1F *h_plot[5];
		
		// only do calculation if sufficient nsamples
		if (_nsamples>NSAMP_CUTOFF)
		{
			if (debug) cout << "\t... calculating S2 parameteres for all PMTs."<< endl;
			
			//for (int j=0; j<5; j++) q_S2[j]= ((2.0/4096.0)*1e-6/50)*get_S2_gaus(h_plot[j],s2_width[j], tau_S2[j], amp_S2[j], ped[j])
			for (int j=0; j<5; j++)	q_S2[j]= q*get_S2(h[j],s2_width[j], tau_S2[j], amp_S2[j],binpeak[j],s1_width[j],ped[j]);
			
			if (debug) cout << "\t... rebinning wf histogram for plotting / S2 calculation."<< endl;
			
			for (int k=0; k<5; k++)
			{
				h[k]->Sumw2();
				h_plot[k] = dynamic_cast<TH1F*>(h[k]->Rebin(128,Form("%s%i","h",k)));
				//h_plot[k]->Sumw2();
				h_plot[k]->Scale(1.0/128.0);
				h_plot[k]->GetYaxis()->SetTitleOffset(1.2);
				h_plot[k]->GetYaxis()->SetTitleSize(0.04);
				h_plot[k]->GetXaxis()->SetTitleSize(0.04);
			}
			
			if (debug) cout << "\t... calculating S2_gaus parameteres for all PMTs."<< endl;
		
			for (int j=0; j<5; j++)	q_S2_gaus[j]= ((2.0/4096.0)*1e-6/50.)*get_S2_gaus(h_plot[j],S2_gaus_width[j], tau_S2_gaus[j], amp_S2_gaus[j],ped[j]);
		}
		
	
		if (debug) 
		{
			cout << "\t... done."<< endl << endl;
			
			cout << "Light: (Run,ev): ("<<_run <<","<<_nevent <<") /ttotlight: " << totlight[5] << " " << _nsamples << endl;
			cout << "S1 light: " << q_S1[0]+q_S1[1]+q_S1[2]+q_S1[3]+q_S1[4] << " S2 light: " << q_S2[0]+q_S2[1]+q_S2[2]+q_S2[3]+q_S2[4]<< endl;
			cout << "Time: charge: " << time_charge << ", light: " << time_light << ", delta: " << time_charge - time_light << endl;
			printf("Double time: charge = %0.12e, light = %0.12e, delta = %0.6e\n",time_charge/1.E9,time_light/1.E9,time_charge/1.E9-time_light/1.E9);
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

		float panel_space=100+20;
		float det_dy = 7302.;
		float y_center = (det_dy + panel_space*2)/2.;
		y_center=1500;
    	float z_center = (1910.+1930.)/2 + 112.*8;
		float fc_zmin = 2129;
		float fc_zmax = 2129+950;
		z_center = (fc_zmin+fc_zmax)/2 ;
		float fv_dz   = 950 + (618.08-418.08);
		float fv_zmin = 2129 -(618.08-418.08);
		float fv_zmax = fv_zmin + fv_dz;
		z_center = (fv_zmin+fv_zmax)/2;

		float x1,y1,z1,x2,y2,z2;
		float crtYX_m=0;
		float crtYX_n=0;

		float crtZX_m=0;
		float crtZX_n=0;
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
		}
		else _crt_matchreco=0;
		//view[0]->Draw("colz");
		//lets_pause();
		//ana_crt_evt(runnum_C,subrun_C,_event,1, gains);
		//view[1]->Draw("colz");
		//lets_pause();
		
		if (1)//totrecoq/totq > 0.3 && _NumberOfTracks_pmtrack==1)// && _crt_ToF<0) _crt_daq_match==1&&_crt_reco==1 
		{
			_run_out=_runcharge;
			_subrun_out=_subruncharge;
			_nev_out=_event;
			_ntracks=_NumberOfTracks_pmtrack;
			//cout << "_ntracks" << _ntracks << endl;

			_time_charge = time_charge;
			_time_light = time_light;

			for (int k=0; k<5; k++) _pmt_wvf_properties[k][0]=ped[k];//ADC Counts
			for (int k=0; k<5; k++) _pmt_wvf_properties[k][1]=pedrms[k]; //ADC Counts
			
			for (int k=0; k<5; k++) _pmt_ped[k]=ped[k];
			for (int k=0; k<5; k++) _pmt_ped_end[k]=ped_end[k];
			for (int k=0; k<5; k++) _pmt_pedrms[k]=pedrms[k];
			for (int k=0; k<5; k++) _pmt_pedrms_end[k]=pedrms_end[k];
		
			for (int k=0; k<6; k++) _pmt_charge[k]=(float)q*totlight[k];
			for (int k=0; k<6; k++) _pmt_npe[k]=(float)q*totlight[k]/(1.602e-19)/gains[k];

			for (int k=0; k<5; k++) _pmt_S1_charge[k]=q_S1[k]; // in volts
			for (int k=0; k<5; k++) _pmt_S1_npe[k]=q_S1[k]/(1.602e-19)/gains[k];

			_pmt_S1_charge[5]=0;
			_pmt_S1_npe[5]=0;

			for (int k=0; k<5; k++) _pmt_S1_charge[5]+=_pmt_S1_charge[k]; // in volts
			for (int k=0; k<5; k++) _pmt_S1_npe[5]+=_pmt_S1_npe[k];

			for (int k=0; k<5; k++) _pmt_S1_width[k]=s1_width[k]*_time_sample; // in ns
			for (int k=0; k<5; k++) _pmt_S1_amp[k]=amp_S1[k]; // in volts
			for (int k=0; k<5; k++) _pmt_S1_tau[k]=tau_S1[k]; // in us

			for (int k=0; k<5; k++) _pmt_S2_charge[k]=q_S2[k];
			for (int k=0; k<5; k++) _pmt_S2_npe[k]=q_S2[k]/(1.602e-19)/gains[k];
			for (int k=0; k<5; k++) _pmt_S2_width[k]=s2_width[k]; // in us
			for (int k=0; k<5; k++) _pmt_S2_amp[k]=amp_S2[k]; // in volts
			for (int k=0; k<5; k++) _pmt_S2_tau[k]=tau_S2[k]; // in us
			
			for (int k=0; k<5; k++) _pmt_S2_gaus_charge[k]=q_S2_gaus[k];
			for (int k=0; k<5; k++) _pmt_S2_gaus_npe[k]=q_S2_gaus[k]/(1.602e-19)/gains[k];
			for (int k=0; k<5; k++) _pmt_S2_gaus_width[k]=S2_gaus_width[k]; // in us
			for (int k=0; k<5; k++) _pmt_S2_gaus_amp[k]=amp_S2_gaus[k]; // in volts
			for (int k=0; k<5; k++) _pmt_S2_gaus_tau[k]=tau_S2_gaus[k]; // in us
			

			_pmt_S2_charge[5]=0;
			_pmt_S2_npe[5]=0;
			for (int k=0; k<5; k++) _pmt_S2_charge[5]+=_pmt_S2_charge[k]; // in volts
			for (int k=0; k<5; k++) _pmt_S2_npe[5]+=_pmt_S2_npe[k];
			
			_pmt_S2_gaus_charge[5]=0;
			_pmt_S2_gaus_npe[5]=0;
			for (int k=0; k<5; k++) _pmt_S2_gaus_charge[5]+=_pmt_S2_gaus_charge[k]; // in volts
			for (int k=0; k<5; k++) _pmt_S2_gaus_npe[5]+=_pmt_S2_gaus_npe[k];

			_crt_mYX = crtYX_m;
			_crt_mZX = crtZX_m;
			_crt_nYX = crtYX_n;
			_crt_nZX = crtZX_n;

			for (int k=0; k<5; k++)_tpc_drift_time_at_pmt_pos[k]=tpc_drift_time_at_pmt_pos[k];
			_tpc_totcharge=totq;
			_tpc_totrecocharge=totrecoq;
//			_tpc_track_charge=trackcharge[0];

			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _tpc_track_charge[j] = trackcharge[sorted_tracks[j]];
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _tpc_track_fitresult_yx[j]=fail_yx[sorted_tracks[j]];
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _tpc_track_fitresult_zx[j]=fail_zx[sorted_tracks[j]];
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _tpc_mYX[j] = t_yx_m[sorted_tracks[j]];
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _tpc_mZX[j] = t_zx_m[sorted_tracks[j]];
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _tpc_nYX[j] = t_yx_n[sorted_tracks[j]];
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _tpc_nZX[j] = t_zx_n[sorted_tracks[j]];

			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _tpc_startX[j] = _Track_StartX_pmtrack[sorted_tracks[j]];
			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _tpc_endX[j] = _Track_EndX_pmtrack[sorted_tracks[j]];


			//for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "TL: " <<  _tpc_nZX[j] << " " << t_zx_n[sorted_tracks[j]] << " " <<t_zx_n[j] << endl;

			for (int j=0; j<_NumberOfTracks_pmtrack; j++) _Track_Length_n[j] = _Track_Length_pmtrack[sorted_tracks[j]];

			//for (int j=0; j<_NumberOfTracks_pmtrack; j++) cout << "TL: " <<  sorted_tracks[j] << " " << _Track_Length_pmtrack[sorted_tracks[j]] << " " <<_Track_Length_pmtrack[j] << endl;


//			std::copy(std::begin(trackcharge), std::end(trackcharge), _tpc_track_charge);
//			std::copy(std::begin(t_yx_m), std::end(t_yx_m), _tpc_mYX);
//			std::copy(std::begin(t_zx_m), std::end(t_zx_m), _tpc_mZX);
//			std::copy(std::begin(t_yx_n), std::end(t_yx_n), _tpc_nYX);
//			std::copy(std::begin(t_zx_n), std::end(t_zx_n), _tpc_nZX);
//			std::copy(std::begin(fail_yx), std::end(fail_yx), _tpc_track_fitresult_yx);
//			std::copy(std::begin(fail_zx), std::end(fail_zx), _tpc_track_fitresult_zx);

//			_tpc_mYX = t_yx_m[0];
//			_tpc_mZX = t_zx_n[0];
//			_tpc_nYX = t_yx_m[0];
//			_tpc_nZX = t_zx_n[0];

			dpd->Fill();
			
		
			// Wf displays				
			if (display_waveforms)
			{
				gStyle->SetOptStat(0);
			    gStyle->SetPadTickX(1);
			    gStyle->SetPadTickY(1);
				
				TCanvas *c1 = new TCanvas("c1","c1",1200,800);
				c1->Divide(3,2);
				for (int k=0; k<5; k++)
				{
					c1->cd(k+1);
					if (_nsamples<=NSAMP_CUTOFF) h[k]->Draw("HIST");
					else h_plot[k]->Draw("HIST");
				
					float x = 0.35;
					float y = 0.45;
					TLatex l;
					l.SetNDC();
					l.DrawLatex(x,y+0.06,Form("Run %d, Event %d",_run_out,ev));
					l.DrawLatex(x,y,Form("ped = %0.2f +/- %0.2f",ped[k],pedrms[k]));
					l.DrawLatex(x,y-0.06,Form("ped_end = %0.2f +/- %0.2f",ped_end[k],pedrms_end[k]));
					l.DrawLatex(x,y-0.12,Form("S1: %0.2f V, %0.1f ns, %0.1f #mus",_pmt_S1_amp[k],_pmt_S1_width[k],_pmt_S1_tau[k]));
					if (_pmt_S2_amp[k]>0.) 
					{
						l.DrawLatex(x,y-0.18,Form("S2: %0.2f mV, %0.1f #mus, %0.1f #mus",_pmt_S2_amp[k]*1000.,_pmt_S2_width[k],_pmt_S2_tau[k]));
						l.DrawLatex(x,y-0.24,Form("S2 gaus: %0.2f mV, %0.1f #mus, %0.1f #mus",_pmt_S2_gaus_amp[k]*1000.,_pmt_S2_gaus_width[k],_pmt_S2_gaus_tau[k]));
					}	
					//l.DrawLatex(x,y-0.30,Form("ped_diff = %0.2f",(ped_end[k]-ped[k])/pedrms[k]));
				}
				c1->Update();
				c1->Modified();
				//c1->Print(Form("working_v2/S1_amp_1/run%d_%d_ev%d.png",_run_out,_subrun_out,_event));	
			}

			if (debug || display_waveforms) lets_pause();
			
		}
		
		htime->Fill(time_charge-time_light);
		hcorr->Fill(totq,totlight[5]);
		for (int k=0; k<5; k++) 
		{ 
			h[k]->Delete(); 
			if (h_plot[k]) h_plot[k]->Delete(); 
		}
		view[0]->Delete();
		view[1]->Delete();	

	} //event loop

	dpd->Write();

	//htime->Draw("HIST");
	//lets_pause();
	//hcorr->Draw("colz");
	//lets_pause();
	//dpd->Draw("pmt_charge[5]:tpc_totcharge");
	//lets_pause();
}

