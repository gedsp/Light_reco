#!/bin/bash

run=$1
subrun=$2
dolight=$3

source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.04.18/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh

cd /afs/cern.ch/user/l/leyton/git/Light_reco/dpdMaker

echo "ls before macro: "
ls --color -lh -c

root -l -q -b "dpdMaker.C($run,$subrun,$dolight)"

#'("""+str(run)+""","""+str(subRun)+""")'      
#root -l -b -q "macrogrid.C+g($i,$j,$k)"

echo "ls after running macro: "
ls --color -lh -c

echo "done"
