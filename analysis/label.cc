#include "RootInterface.h"
#include "RecoInterface.h"
#include "DRsimInterface.h"
#include "fastjetInterface.h"
#include "functions.h"

#include "GeoSvc.h"
#include "GridDRcalo.h"
#include "fastjet/PseudoJet.hh"

#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVectorF.h"
#include "TMatrixFSym.h"

#include <string>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;
//./process ../../box/pi0_0 ../../test_0.root 0
//./process ../../elenshower ../../elenshower.root 0
int main(int , char* argv[]){
  TString innameform = argv[1];
  TString outname = argv[2];
  int isjet = atoi(argv[3]);
  TString inname;
  TString reconame;

  int check_hep=0;
  int check_dr=0;
  int check_reco=0;
  int count=0;
  int pass=0;
  int code=0;
  TList* keys;
  TFile* box;
  int readcount=0;

  std::ofstream writeFile;
  writeFile.open(outname.Data());


    inname=innameform;
    reconame=inname.Copy();
    reconame.ReplaceAll(".root","_Reco.root");
    check_hep=0;
    check_dr=0;
    check_reco=0;
    pass=0;
    code=0;
    if(!fs::exists(inname.Data())){
      //printf("nofile %s;\n",inname.Data());
      pass=1;
    }
    if(pass==0){
      //printf("filechecked2 %s;\n",inname.Data());
      box=TFile::Open(inname.Data(),"read");
      keys = box->GetListOfKeys();
      //printf("filechecked3 %s;\n",inname.Data());
      if(keys->GetEntries()>0){
        for(int i =0;i<keys->GetEntries();i++){
          if(strcmp(keys->At(i)->GetName(),"hepmc3_tree")==0){
            code=1;
            check_hep=1;
          }
          if(strcmp(keys->At(i)->GetName(),"DRsim")==0){
            code=2;
            check_dr=1;
          }
          if(strcmp(keys->At(i)->GetName(),"Reco")==0){
            //gDirectory
            //gDirectory->Delete("Reco;1");
          }
        }
        if(isjet==1){
          if(check_hep==0){
            //printf("%s no hep\n",inname.Data());
            pass=1;
          }
        }
        if(check_dr==0){
          //printf("%s no dr\n",inname.Data());
          pass=1;
        }
      }
      else{
        //printf("%s no key\n",inname.Data());
        pass=1;
      }
    keys->Clear();
    box->Close();
      if(pass==0){
        if(!fs::exists(reconame.Data())){
          //printf("nofile %s;\n",reconame.Data());
          pass=1;
        }
        if(pass==0){
          box=TFile::Open(reconame.Data(),"read");
          keys = box->GetListOfKeys();
          if(keys->GetEntries()>0){
            for(int i =0;i<keys->GetEntries();i++){
              if(strcmp(keys->At(i)->GetName(),"Reco")==0){
                check_reco=1;
              }
            }
            if(check_reco==0){
              printf("%s no reco\n",reconame.Data());
              pass=1;
            }
          }
          else{
            printf("%s no key\n",reconame.Data());
            pass=1;
          }
        }
      }
    }
    //printf("load %s...;\n",inname.Data());
    //mychain.Add(inname+"/DRsim");
    if(pass==0){
      printf("%s\n",inname.Data());
      writeFile.write(std::string(inname.Data()).c_str(), std::string(inname.Data()).size());
      writeFile.write("\n", 1);
    }
writeFile.close();

return 0;
}
