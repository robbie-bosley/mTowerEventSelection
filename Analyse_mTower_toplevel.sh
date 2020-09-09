#!/bin/bash

#run with run number as argument, e.g. ./Analyse_mTower.sh 1250

#to submit this script in the quark queue
#        qsub -cwd -V ./Analyse_mTower.sh 1250
# -cwd uses the current working directory
# -V keeps the environment variables

echo "====================================="
echo "starting job Analyse_mTower($1)"
echo "====================================="

root -l -b <<SHELL
.!pwd
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionPerEvent/classes/mTowerHit.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionPerEvent/classes/mTowerEvent.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionPerEvent/classes/mTowerCluster.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionPerEvent/classes/mTowerChip.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionPerEvent/subprocessors/EventSelection.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionPerEvent/Analyse_mTower_toplevel.cxx+
Analyse_mTower(1336)
.q
SHELL
