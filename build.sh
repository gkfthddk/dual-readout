source setenv-cc7-gcc8.sh
source init_lcg.sh
mkdir build
cd build
mkdir box
cp ../DRsim/run_ele.mac DRsim/
cp ../DRsim/run_pi.mac DRsim/
cmake ..
make -j8
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HEPMC_DIR/lib64:$FASTJET_DIR/lib:$PYTHIA_DIR/lib:$PWD/lib
source bin/thisddDRcalo.sh
cp rootIO/{librootIO_rdict.pcm,librootIO.rootmap} lib/
cp -r ../Detector/DRcalo/compact bin/
mkdir -p DRsim/bin
mkdir -p Reco/bin
cp -r ../Detector/DRcalo/compact Reco/bin/
cp -r ../Detector/DRcalo/compact DRsim/bin/

cp bin/DRsim DRsim/
cp bin/Reco Reco/

./DRsim run_ele.mac 0 testel
