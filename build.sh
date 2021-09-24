source setenv-cc7-gcc8.sh
source init_lcg.sh
mkdir build
cd build
mkdir box
cp ../DRsim/run_ele.mac DRsim/
cp ../DRsim/run_pi.mac DRsim/
cp ../DRsim/run_en.mac DRsim/
cmake ..
make -j8
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HEPMC_DIR/lib64:$FASTJET_DIR/lib:$PYTHIA_DIR/lib:$PWD/lib
source bin/thisddDRcalo.sh
cp rootIO/{librootIO_rdict.pcm,librootIO.rootmap} lib/
cp -r ../Detector/DRcalo/compact bin/
mkdir -p DRsim/bin
mkdir -p Reco/bin
mkdir -p analysis/bin
cp -r ../Detector/DRcalo/compact Reco/bin/
cp -r ../Detector/DRcalo/compact DRsim/bin/
cp -r ../Detector/DRcalo/compact analysis/bin/

cp bin/P8 Gen/
cp bin/DRsim DRsim/
cp bin/Reco Reco/
cp bin/process analysis/

cd ..
unset PYTHONPATH
unset PYTHONHOME
