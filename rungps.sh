#!/bin/bash
CMND=$1
BOX=$3
SEED=$(( $2 + $BOX ))
filename=$4
PARTICLE=$5
MINENERGY=$6
MAXENERGY=$7
num_EVENT=$8
num_CPU=$9
WORKPLACE=${PWD}
date
begin=$( date )
mkdir box
hostname
ls /pad
echo "${WORKPLACE}/${filename}_${SEED}.root batch ${num_EVENT} from ${BOX}"

if grep -Fxq "${SEED}" build/${filename##*/}_skip.txt
then
    echo "skip seed ${SEED}"
    return 1 
    echo "exitting"
    exit 1
fi

echo "passed"
cd build
echo "pwd"
pwd
export LD_LIBRARY_PATH=${PWD}/lib:$LD_LIBRARY_PATH
cd DRsim
echo "/run/numberOfThreads ${num_CPU}" > buf.mac
cat ${CMND}.mac >> buf.mac
echo "/gps/particle ${PARTICLE}" >> buf.mac
#echo "/gun/energy $( seq 20 0.01 50 | shuf | head -n1 ) GeV" >> buf.mac
echo "/gps/ene/min ${MINENERGY} GeV" >> buf.mac
echo "/gps/ene/max ${MAXENERGY} GeV" >> buf.mac
echo "/run/beamOn ${num_EVENT}" >> buf.mac
cat buf.mac
echo "#macro end"
echo "DRsim"
./DRsim buf.mac ${SEED} ${WORKPLACE}/${filename}
date
cd ../Reco
echo "Reco"
./Reco ${SEED} ${WORKPLACE}/${filename}
cd ${WORKPLACE}
echo "${WORKPLACE}/${filename}_${SEED}.root from ${BOX}"
echo "batch ${num_EVENT}"
echo "begin" $begin
echo "end" $( date )
