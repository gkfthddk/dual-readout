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
  int num_file = std::stoi(argv[2]);
  TString outname = argv[3];  
  int isjet = atoi(argv[4]);
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

  std::ifstream readFile;
  std::ofstream writeFile;
  readFile.open(innameform.Data());
  writeFile.open(outname.Data());
  while (!readFile.eof()){
    readcount+=1;
    std::string str;
    getline(readFile, str);

  std::string delimiter = " ";

  size_t pos = 0;
  pass=0;
  std::string token;
  while ((pos = str.find(delimiter)) != std::string::npos) {
        token = str.substr(0, pos);
            //std::cout <<"token "<< token << std::endl;
                str.erase(0, pos + delimiter.length());
  //std::cout <<"str "<< str << std::endl;
  }


    inname=TString(str);
    reconame=inname.Copy();
    reconame.ReplaceAll(".root","_Reco.root");
    //printf("filechecking %s;\n",inname.Data());
    check_hep=0;
    check_dr=0;
    check_reco=0;
    pass=0;
    code=0;
    if(!fs::exists(inname.Data())){
      printf("nofile %s;\n",inname.Data());
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
            code+=1;
            check_hep=1;
          }
          if(strcmp(keys->At(i)->GetName(),"DRsim")==0){
            code+=10;
            check_dr=1;
          }
          if(strcmp(keys->At(i)->GetName(),"Reco")==0){
            code+=100;
            check_reco=1;
          }
        }
        if(isjet==1){
          if(check_hep==0){
            printf("%s no hep\n",inname.Data());
            pass=1;
          }
        }
        if(check_dr==0){
          printf("%s no dr\n",inname.Data());
          pass=1;
        }
      }
      else{
        printf("%s no key\n",inname.Data());
        pass=1;
      }
      if(pass==0){
        keys->Clear();
        box->Close();
        if(!fs::exists(reconame.Data())){
          printf("nofile %s;\n",reconame.Data());
          pass=1;
        }
        if(pass==0){
          box=TFile::Open(reconame.Data(),"read");
          keys = box->GetListOfKeys();
          if(keys->GetEntries()>0){
            for(int i =0;i<keys->GetEntries();i++){
              if(strcmp(keys->At(i)->GetName(),"Reco")==0)check_reco=1;
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
    //if(pass==0){
      writeFile.write(std::string(inname.Data()).c_str(), std::string(inname.Data()).size());
      writeFile.write(";", 1);
      writeFile.write(std::to_string(code).c_str(), std::to_string(code).size());
      writeFile.write("\n", 1);
    //}
      
  }
readFile.close();
writeFile.close();
printf("--done\n");
printf("%d files %d events filled in %s\n",readcount,count,outname.Data());

return 0;
}
