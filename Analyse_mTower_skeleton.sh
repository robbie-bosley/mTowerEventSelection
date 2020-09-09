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
.L ./classes/mTowerHit.cxx+
.L ./classes/mTowerEvent.cxx+
.L ./classes/mTowerCluster.cxx+
.L ./classes/mTowerChip.cxx+
.L ./subprocessors/EventSelection.cxx+
.L ./Analyse_mTower_skeleton.cxx+
Analyse_mTower(1336)
.q
SHELL
