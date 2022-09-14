# dual-readout
Repository for GEANT4 simulation &amp; analysis of the dual-readout calorimeter.

    git checkout -t origin/0023_img

## How-to
### Compile
After fetching the repository, do

    source setenv-cc7-gcc8.sh
    

### Install
For a case that needs to install the package (e.g. `condor` requires file transfer), one can install the package via

    /bin/bash build.sh

### Running GEANT4
In build/Gen, 

    ./P8generic <run_command> <filenumber> <filename>
       
or
    
    ./P8ptcgun <run_command> <filenumber> <filename> <MAX_ENERGY> <MOD>

generates `<filename>_<filenumber>.root`

### Running GEANT4
In build/DRsim, 

    ./DRsim <run_macro> <filenumber> <filename>

generates DRsim tree in `<filename>_<filenumber>.root`

### Reconstruction
This requires the ROOT file generated from `DRsim`. Assuming the name of the file `<filename>_<filenumber>.root`, in build/Reco,

    ./Reco <filenumber> <filename>
    
generates, `<filename>_<filenumber>_Reco.root`

### Analysis
This requires the ROOT file generated from `Reco`. Assuming the name of the file `<filename>_0.root ~ <filename>_9.root`, in build/analysis, out.root will be created.

    ./process <filename>_%d.root 0 10 <out> 0 

when jet case

    ./process <filename>_%d.root 0 10 <out> 1
  

source runen.sh run_en 0 0 box/ele e- 20 100

condor_submit run_ele.co

version 0.0.2.3

### Saving image
This requires the ROOT file generated from `process`. Edit `tonpz.py` variables `name_file,name_cls,save_path` with your file name and path.
The script is based on python3, you would run at clean environment (open another terminal session do not source simulation environment).

     python3 tonpz.py
