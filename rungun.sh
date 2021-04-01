#!/bin/bash
CMND=$1
BOX=$3
SEED=$(( $2 + $BOX ))
filename=$4
WORKPLACE=${PWD}
mkdir box
date
hostname
#ls /pad
#cd /pad/yulee/dream/build
cd build
echo "pwd"
pwd
export LD_LIBRARY_PATH=${PWD}/lib$LD_LIBRARY_PATH
cd bin
echo "DRsim"
./DRsim ${CMND}.mac ${SEED} ${WORKPLACE}/${filename}
cd ../Reco
echo "Reco"
./Reco ${SEED} ${WORKPLACE}/${filename}
cd ${WORKPLACE}
date
