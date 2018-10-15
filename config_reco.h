#ifndef CONFIG_RECO_H
#define CONFIG_RECO_H

static const int N_PMT = 5;
static const int N_LEM = 12;

//set up the input folder where the light-charge matched data is stored.
const string matched_data_dir = "/Volumes/data/root_files/2018Feb05/"; // local 
//const string matched_data_dir = "/eos/user/j/jsotooto/root_files/2018Feb05/"; // on lxplus

// full path to Light_reco directory                        
const string reco_dir = "/Volumes/data/Dropbox/physics/311/git/Light_reco"; // local
//const string reco_dir = "/afs/cern.ch/user/l/leyton/git/Light_reco"; // on lxplus 

// full path to DPD storage
//const string dpd_dir = "/afs/cern.ch/user/l/leyton/public/dpd_20180914"; // on lxplus
const string dpd_dir = reco_dir+"/dpdMaker/dpd"; // local

// full path to db file
const string db_file = reco_dir+"/dbVoltages.root";

// full path to pedestal correction file
const string pc_file = reco_dir+"/PedestalCorrection.root";


#endif