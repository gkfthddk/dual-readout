# dual-readout
Repository for GEANT4 simulation &amp; analysis of the dual-readout calorimeter.

## How-to
### Compile
After fetching the repository, do

    source setenv-cc7-gcc8.sh
    

### Install
For a case that needs to install the package (e.g. `condor` requires file transfer), one can install the package via

    source build
    

### Running GEANT4
#### GEANT4 standalone particle gun
In build/DRsim,

    ./DRsim <run_macro> <filenumber> <filename>

generates, `<filename>_<filenumber>.root`

### Reconstruction
This requires the ROOT file generated from `DRsim`. Assuming the name of the file `<filename>_<filenumber>.root`, in build/Reco,

    ./Reco <filenumber> <filename>

### Analysis
This requires the ROOT file generated from `Reco`. Assuming the name of the file `<filename>_<filenumber>.root`, in build/analysis,

    ./<your_analysis_program> <filenumber> <filename>

source runen.sh run_en 0 0 box/ele e- 20 100

condor_submit run_ele.co
