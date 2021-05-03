# dual-readout
Repository for GEANT4 simulation &amp; analysis of the dual-readout calorimeter.

## How-to
### Compile
After fetching the repository, do

    source setenv-cc7-gcc8.sh
    cmake3 .
    make -j4

There are compiled executables in bin/
Copy executable file in bin to each directory.

    cp bin/DRsim DRsim/
    cp bin/Reco Reco/
    cp bin/mergy analysis/
    cp bin/process analysis/
    
When you change code `make` again and update executable from `/bin` directory.

### Running GEANT4
#### 1. GEANT4 standalone particle gun
In DRsim,

    #./DRsim <run_macro> <filenumber> <filename>
    ./DRsim run_ele.mac 0 ../eltest

generates, `<filename>_<filenumber>.root`
I you mind `lib/librootIO_rdict.pcm file does not exist`, copy `rootIO/librootIO_rdict.pcm` to `lib`

#### 1.1 Using HepMC input (need to run Gen first)
This requires the ROOT file generated from `Gen`. Assuming the name of the file `<filename>_<filenumber>.root`,

    #./DRsim run_hepmc.mac <filenumber> <filename>  

Modify `#.mac` for other paticles or conditions

### Reconstruction
This requires the ROOT file generated from `DRsim`. Assuming the name of the file `<filename>_<filenumber>.root`, in build/Reco,

    #./Reco <filenumber> <filename>
    ./Reco 0 ../eltest

### generation macros

    #source rungun.sh <run_macro> <filenumber> <filename>
    source rungun.sh run_ele 0 box/ele

This generates box/ele_0.root with GEANT4 simulation

Run multiple jobs with with condor.
    
    mkdir log
    condor_submit runel.co

This will make box/ele_#.root, # from 0 to 9.
Number of jobs can be changed by editing `runel.co`.
Check status of jobs by
    
    condor_q
    
`Idle` - jobs not yet running, `Run` - jobs running, `Done` - completed jobs, `Hold` - jobs holded having some problem.
Check `log/` if jobs have problem.

### Analysis
This requires the ROOT file generated from `Reco`. Assuming the name of the file `<filename>_<filenumber>.root`, in build/analysis,

    ./<your_analysis_program> <filenumber> <filename>

Process image with files from condor jobs with name `box/ele_#.root`(# from 0 to 9) to `eltest.root`

    ./process ../box/ele_%.root 10 ../eltest

eltest.root have "event" tree contains processed data by each event.

for python contents can be inspected by accessing tree

    import ROOT
    import numpy as np
    
    infile=ROOT.TFile("eltest.root","read")
    event=infile.Get("event")
    event.Print()
    
    event.GetEntry(0)
    print(event.E_DRcorr,np.array(image_ecor_s).shape,np.array(fiber_ecor_s).shape))
    
images have 168 * 168 shape, but imaging algorithm is still in progress.

### Precaution
Since GEANT4 takes very large amount of time per an event, P8ptcgun, DRsim and Reco are assumed to run a few events only per ROOT file. The executables can be run on parallel using `torque` or `condor`, and can be merged before analysis step using `hadd` from ROOT.
