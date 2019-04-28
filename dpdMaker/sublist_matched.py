#!/usr/bin/env python                                                                                                                                                                                     
#                                                                                                                                                                                                         
import sys
import os
import glob
#import ROOT                                                                                                                                                                                              

if __name__ == "__main__":

    datafile = file("../config_reco.h")
    filedir=''
    dpddir=''
    for line in datafile:
        if line.startswith('//'):
                continue
        if 'matched_data_dir' in line:
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

    myfiles=glob.glob(InputPath+"/*.root")

    fileDone = open("matched_runs_done.txt", 'w')
    fileTBC  = open("matched_runs_tbc.txt", 'w')

    for i in myfiles:
        i=i.split(InputPath)[1]
        run=int(i.split("-")[0])
        subRun=int(i.split("-")[1])

        file=DpdPath+'/matched/dpd-matched-'+str(run)+'-'+str(subRun)+'.root'
        exists = os.path.isfile(file)
        if exists:
            sys.stdout = fileDone                                           
        else:
            sys.stdout = fileTBC
                                                                                                                                                                                          
        print run,subRun,"0"

    sys.stdout = sys.__stdout__
    fileDone.close()
    fileTBC.close()

