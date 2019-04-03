#!/usr/bin/env python
#
import sys
import os
import glob
#import ROOT

if __name__ == "__main__":
        
    datafile = file("../config_reco.h")
    filedir=''
    for line in datafile:
        if line.startswith('//'):
                continue
        if 'matched_data_dir' in line:
                filedir=line.split('=',2)[1]
                filedir=filedir.split(';',2)[0]
                filedir=filedir.replace(' ','')
                filedir=filedir.replace('\"','')
                
    RunPath=os.getcwd()
    InputPath=filedir

    #RunPath = "/afs/cern.ch/user/l/leyton/311/macros/macro_analysis"
    #InputPath = "/eos/user/j/jsotooto/root_files/2018Feb05/"
    
    print "RunPath %s" % ( RunPath )
    print "InputPath %s" % ( InputPath )
    print ""
    
    myfiles=glob.glob(InputPath+"/*.root")
    
    file = open("mychargeruns.txt", 'w')
    sys.stdout = file
    
    for i in myfiles:
        i=i.split(InputPath)[1]
        run=int(i.split("-")[0])
        subRun=int(i.split("-")[1])
        print run,subRun

    sys.stdout = sys.__stdout__
    file.close()
