#!/bin/sh

export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.14.2/Linux-x86_64/bin:$PATH
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/8/x86_64-centos7/setup.sh

source /cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/Boost/1.72.0/x86_64-centos7-gcc8-opt/Boost-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/ROOT/v6.20.06/x86_64-centos7-gcc8-opt/ROOT-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/Geant4/10.06.p02/x86_64-centos7-gcc8-opt/Geant4-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/DD4hep/01-12-01/x86_64-centos7-gcc8-opt/DD4hep-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/DD4hep/01-12-01/x86_64-centos7-gcc8-opt/bin/thisdd4hep.sh

export BOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/Boost/1.72.0/x86_64-centos7-gcc8-opt
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/hepmc3/3.2.2/x86_64-centos7-gcc8-opt
export FASTJET_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/fastjet/3.3.2/x86_64-centos7-gcc8-opt
export FASTJET_ROOT_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/fastjet/3.3.2/x86_64-centos7-gcc8-opt
export PYTHIA_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/MCGenerators/pythia8/244/x86_64-centos7-gcc8-opt
export PYTHIA8_ROOT_DIR=/cvmfs/sft.cern.ch/lcg/views/LCG_97a_FCC_1/x86_64-centos7-gcc8-opt/
export XERCESC_ROOT_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/XercesC/3.1.3/x86_64-centos7-gcc8-opt
export CMAKETOOLS_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/cmaketools/1.8/x86_64-centos7-gcc8-opt/share/CMakeTools
#/cvmfs/sft.cern.ch/lcg/views/LCG_97a_FCC_1/x86_64-centos7-gcc8-opt/share/CMakeTools/CMakeToolsConfig.cmake

export PYTHIA8=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/MCGenerators/pythia8/244/x86_64-centos7-gcc8-opt
export PYTHIA8DATA=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/MCGenerators/pythia8/244/x86_64-centos7-gcc8-opt/share/Pythia8/xmldoc
export ROOT_INCLUDE_PATH=$HEPMC_DIR/include:$ROOT_INCLUDE_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HEPMC_DIR/lib64:$FASTJET_DIR/lib:$PYTHIA_DIR/lib


source $HEPMC_DIR/hepmc3-env.sh
source $FASTJET_DIR/fastjet-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/cmaketools/1.8/x86_64-centos7-gcc8-opt/cmaketools-env.sh

source init_lcg.sh
