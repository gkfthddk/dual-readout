#!/bin/bash
CMND=$1
BOX=$3
SEED=$(( $2 + $BOX ))
filename=$4
PARTICLE=$5
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
cp ${CMND}.mac buf.mac
echo "/gun/energy $( seq 20 .01 50.01 | shuf | head -n1 ) GeV" >> buf.mac
echo "/gun/particle ${PARTICLE}" >> buf.mac
cat buf.mac
echo "DRsim"
./DRsim buf.mac ${SEED} ${WORKPLACE}/${filename}
echo "Reco"
./Reco ${SEED} ${WORKPLACE}/${filename}
cd ${WORKPLACE}
date
