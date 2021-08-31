#!/bin/bash
CMND=$1
BOX=$3
SEED=$(( $2 + $BOX ))
filename=$4
PARTICLE=$5
ENERGY=$6
num_EVENT=$7
WORKPLACE=${PWD}
mkdir box
date
hostname
#ls /pad
#cd /pad/yulee/dream/build
cd build
echo "pwd"
pwd
export LD_LIBRARY_PATH=${PWD}/lib:$LD_LIBRARY_PATH
cd DRsim
cp ${CMND}.mac buf.mac
echo "/gun/particle ${PARTICLE}" >> buf.mac
#echo "/gun/energy $( seq 20 0.01 50 | shuf | head -n1 ) GeV" >> buf.mac
echo "/gun/energy ${ENERGY} GeV" >> buf.mac
echo "/run/beamOn ${num_EVENT}" >> buf.mac
cat buf.mac
echo "DRsim"
./DRsim buf.mac ${SEED} ${WORKPLACE}/${filename}
cd ../Reco
echo "Reco"
./Reco ${SEED} ${WORKPLACE}/${filename}
cd ${WORKPLACE}
date
