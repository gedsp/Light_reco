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
    
    os.system("mkdir -p run")
    os.system("mkdir -p run/scripts")
    os.system("mkdir -p run/log")
    os.system("mkdir -p run/out")
    os.system("mkdir -p dpd")
    
    myfiles=glob.glob(InputPath+"/*.root")
    
    for i in myfiles:
        i=i.split(InputPath)[1]
        run=int(i.split("-")[0])
        subRun=int(i.split("-")[1])
        #print run,subRun
        #print run

        #if run!=840:
        #    continue
            
        # write batch script
        runName = "run-"+str(run)+"-"+str(subRun)
        file = open("run/scripts/"+runName+".sh", 'w')
        sys.stdout = file

        print """
hostname
date
cd """+RunPath+"""
pwd
source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.06/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
root -l -q -b dpdMaker.C'("""+str(run)+""","""+str(subRun)+""")'
"""
        
        sys.stdout = sys.__stdout__
        file.close()

        # make bash script executable
        os.system("chmod +x run/scripts/" + runName+".sh")

        # cmd to run and preserve LSF output
        #bsubcmd= "bsub -q 1nh "+RunPath+"/run/scripts/"+runName+".sh -e "+RunPath+"/run/log/"+runName+".log -o "+RunPath+"run/out/"+runName+".out"
    
        # cmd to run and throw away output 
        bsubcmd= "bsub -q 1nh -o /dev/null -e /dev/null "+RunPath+"/run/scripts/"+runName+".sh" 

        print bsubcmd
        os.system(bsubcmd)
