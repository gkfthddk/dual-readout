cp bin/DRsim DRsim/
cp bin/Reco Reco/
mkdir DRsim/bin
mkdir Reco/bin
mkdir analysis/bin
cp -r Detector/DRcalo/compact build/DRsim/bin/
cp -r Detector/DRcalo/compact build/Reco/bin/
cp -r Detector/DRcalo/compact build/analysis/bin/
source bin/thisddDRcalo.sh
