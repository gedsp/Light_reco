#!/usr/bin/env python                                                                                                                                                                                     
#                                                                                                                                                                                                         
import sys
import os
import glob
import ROOT                                                                                                                                                                                              

if __name__ == "__main__":

    datafile = file("../config_reco.h")
    filedir=''
    dpddir=''
    for line in datafile:
        if line.startswith('//'):
                continue
        if 'light_data_dir' in line:
                filedir=line.split('=',2)[1]
                filedir=filedir.split(';',2)[0]
                filedir=filedir.replace(' ','')
                filedir=filedir.replace('\"','')
        if 'dpd_dir' in line:
                dpddir=line.split('=',2)[1]
                dpddir=dpddir.split(';',2)[0]
                dpddir=dpddir.replace(' ','')
                dpddir=dpddir.replace('\"','')

    RunPath=os.getcwd()
    InputPath=filedir 
    DpdPath=dpddir                                  

    print "RunPath %s" % ( RunPath )
    print "InputPath %s" % ( InputPath )
    print "DpdPath %s" % ( DpdPath )
    print ""

    myfiles=glob.glob(InputPath+"/output0000*_reprocessed.root")

    fileDone = open("light_runs_done.txt", 'w')
    fileTBC  = open("light_runs_tbc.txt", 'w')

    for i in myfiles:
        fullpathname=i
        fsize=os.path.getsize(fullpathname)
        i=i.split(InputPath)[1]
        blah=i.split("_reprocessed.root")[0]
        blah2=blah.split("0000")[1]
        run=-1
        subRun=-1
        if (blah2.isdigit()):
            run=int(blah2)
        if run<0:
            continue
        #if run>=1414:
        #    continue
        
        chain = ROOT.TChain("midas_data")
        chain.Add(fullpathname)
        entries = chain.GetEntries()
        if entries<1000 or fsize/entries<100000:
            file=DpdPath+'/light/dpd-light-'+str(run)+'.root'
            exists = os.path.isfile(file)
            if exists:
                sys.stdout = fileDone                                           
            else:
                sys.stdout = fileTBC
            print run,subRun,"1"
            continue
        
        for subRun in range (0,int(entries/1000)+1):
            file=DpdPath+'/light/dpd-light-'+str(run)+'-'+str(subRun)+'.root'
            exists = os.path.isfile(file)
            if exists:
                sys.stdout = fileDone                                           
            else:
                sys.stdout = fileTBC
                                                                                                                                                                                          
            print run,subRun,"1"

    sys.stdout = sys.__stdout__
    fileDone.close()
    fileTBC.close()

