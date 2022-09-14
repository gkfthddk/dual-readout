source setenv-cc7-gcc8.sh
source init_lcg.sh
#source sh setenv-cc7-gcc8.sh
#source sh init_lcg.sh
cd build
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HEPMC_DIR/lib64:$FASTJET_DIR/lib:$PYTHIA_DIR/lib:$PWD/lib
source bin/thisddDRcalo.sh
#source sh bin/thisddDRcalo.sh
cd ..
unset PYTHONPATH
unset PYTHONHOME
