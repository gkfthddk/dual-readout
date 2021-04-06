#!/bin/bash
CMND=$1
BOX=$3
SEED=$(( $2 + $BOX ))
filename=$4
maxE=$5
WORKPLACE=${PWD}
mkdir box
date
#cd /pad/yulee/dream/build
cd build
export LD_LIBRARY_PATH=${PWD}/lib$LD_LIBRARY_PATH
cd bin
echo "Gen"
./P8ptcgen ${CMND}.cmnd ${SEED} ${WORKPLACE}/${filename} ${maxE}
echo "filename ${filename}"
echo "./P8ptcgen ${CMND}.cmnd ${SEED} ${WORKPLACE}/${filename} ${maxE}"
echo "DRsim"
./DRsim run_hepmc.mac ${SEED} ${WORKPLACE}/${filename}
echo "Reco"
./Reco ${SEED} ${WORKPLACE}/${filename}
cd ${WORKPLACE}
ls *
date
