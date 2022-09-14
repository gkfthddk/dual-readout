#!/bin/bash
CMND=$1
BOX=$3
SEED=$(( $2 + $BOX ))
filename=$4
WORKPLACE=${PWD}
num_EVENT=$5
P8=$6
MAX_ENERGY=$7
MOD=$8
date
begin=$( date )
mkdir box
hostname
ls
echo "${WORKPLACE}/${filename}_${SEED}.root batch ${num_EVENT} from ${BOX}"
cd build
echo "pwd" $PWD
export LD_LIBRARY_PATH=${PWD}/lib:$LD_LIBRARY_PATH
cd Gen
cp ../bin/P8* .
echo "Gen"
echo "Gen" >> /dev/stderr
date
#./${P8} ${CMND}.cmnd ${SEED} ${WORKPLACE}/${filename}
./${P8} ${CMND}.cmnd ${SEED} ${WORKPLACE}/${filename} ${MAX_ENERGY} ${MOD}
cd ../DRsim
cp ../bin/DRsim .
echo "./DRsim run_${CMND}.mac ${SEED} ${WORKPLACE}/${filename}"
echo "DRsim" >> /dev/stderr
date
./DRsim run_${CMND}.mac ${SEED} ${WORKPLACE}/${filename}
cd ../Reco
cp ../bin/Reco .
echo "./Reco ${SEED} ${WORKPLACE}/${filename}"
echo "Reco" >> /dev/stderr
date
./Reco ${SEED} ${WORKPLACE}/${filename}
cd ${WORKPLACE}
echo "ls box" $(ls box)
echo "${WORKPLACE}/${filename}_${SEED}.root from ${BOX}"
cd build/bin
./buf ${WORKPLACE}/${filename}_${SEED}.root ${WORKPLACE}/pass_${filename}_${SEED} 1
ls ${WORKPLACE}
cat ${WORKPLACE}/pass_${filename}_${SEED}
echo "batch ${num_EVENT}"
echo "begin" $begin
echo "end" $( date )
