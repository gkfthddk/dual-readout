# dual-readout
Repository for GEANT4 simulation &amp; analysis of the dual-readout calorimeter.

    git clone git@github.com:gkfthddk/dual-readout.git
    git pull origin 0022
    git checkout 0022

## How-to
### Compile
After fetching the repository, do (modify setuos.sh or use setenv*.sh at your system)

    source setuos.sh
    mkdir build
    cd build
    cmake3 ..
    make -j4
    cd ..
    source aftermake.sh # copy librootIO to lib/ and set ddhep environment

### Running GEANT4
particle gun macro script (modify rungun.sh and add *.mac at build/bin at your setting)
DRsim and Reco to test_0.root

    source rungun.sh run_ele 0 0 test
    
Condor with macro script, output files stored at box/, log at condor/

    mkdir box
    mkdir condor
    condor_submit runel.co


#### 1. GEANT4 standalone particle gun
In build/DRsim,

    cd build/bin
    ./DRsim <run_macro> <filenumber> <filename>

generates, `<filename>_<filenumber>.root`

#### 2. Using HepMC input
This requires the ROOT file generated from `Gen`. Assuming the name of the file `<filename>_<filenumber>.root`,

    ./DRsim run_hepmc.mac <filenumber> <filename>

### Reconstruction
This requires the ROOT file generated from `DRsim`. Assuming the name of the file `<filename>_<filenumber>.root`, in build/Reco,

    cd build/bin
    ./Reco <filenumber> <filename>
