#ifndef CONFIG_RECO_H
#define CONFIG_RECO_H

static const int N_PMT = 5;
static const int N_LEM = 12;

static const int KMAXNSEG=10;
static const int KMAXNPEAKS=50;
static const int KMAXNTRACKS=200;
static const int KMAXNEVENTS = 350; // max events in a charge subrun (>335)
static const int KMAXNHITS=50000;
static const float ERRVAL=-9999;

// PMT positions in LArSoft coordinate system (in cm)
static const float pmt_x_pos = -698.0;
static const float pmt_y_pos = 0.0;
static const float pmt_z_pos[N_PMT]={50.0,100.0,150.0,200.0,250.0};

// set up the input folder where the light-charge matched data is stored.
const string matched_data_dir = "/eos/user/l/leyton/dpd_20190419/matched/"; // local
//const string matched_data_dir = "/eos/user/l/leyton/2018Aug24/"; // on lxplus

// set up the input folder where the light-only data is stored.
const string light_data_dir = "/eos/experiment/wa105/data/311_PMT/data/root/reprocessed_5apr19/"; // local
//const string light_data_dir = "/eos/user/c/chlastor/ReprocessedData/21Feb/"; // on lxplus

// set up the input folder where the Highway Algorithm data is stored 
const string highway_data_dir ="/afs/cern.ch/user/g/gldesape/WA105_mine/Analysis/Event-track-selection/HighwayAlgorithm"; // local
//const string highway_data_dir = "/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_August_24/Highway_VBW/"; // on lxplus

// full path to Light_reco directory                        
const string reco_dir = "/afs/cern.ch/user/g/gldesape/WA105_mine/Light_reco"; // local
//const string reco_dir = "/afs/cern.ch/user/l/leyton/git/Light_reco"; // on lxplus 

// full path to DPD storage
//const string dpd_dir = "/afs/cern.ch/user/l/leyton/work/dpd_20190403"; // on lxplus
const string dpd_dir = "/afs/cern.ch/user/g/gldesape/WA105_mine/Light_reco/dpdMaker/dpd"; // local

// full path to db files
const string db_charge_file = reco_dir+"/dbVoltages_charge.root";
const string db_light_file = reco_dir+"/dbVoltages_light.root";

#endif
