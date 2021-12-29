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
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TPaveStats.h"
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


  int count=0;
  int pass=0;
   TFile* box;
  TList* keys;
  int readcount=0;

  std::ifstream readFile;
  std::ofstream writeFile;
  readFile.open(innameform.Data());
  writeFile.open(outname.Data());
  while (!readFile.eof()){
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
    //printf("filechecking %s;\n",inname.Data());
    
    if(!fs::exists(inname.Data())){
      printf("nofile1 %s;\n",inname.Data());
      continue;}
    //printf("filechecked2 %s;\n",inname.Data());
    box=TFile::Open(inname.Data(),"read");
    keys = box->GetListOfKeys();
    
    //printf("filechecked3 %s;\n",inname.Data());
    if(keys->GetEntries()>1){
      if(isjet==1){
        if(strcmp(keys->At(0)->GetName(),"hepmc3_tree")==0 && strcmp(keys->At(keys->GetEntries()-2)->GetName(),"DRsim")==0 && strcmp(keys->At(keys->GetEntries()-1)->GetName(),"Reco")==0){
    RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(inname.Data()));
    RootInterface<RecoInterface::RecoEventData>* recoInterface = new RootInterface<RecoInterface::RecoEventData>(std::string(inname.Data()));
    HepMC3::ReaderRootTree reader(std::string(inname.Data()));
    //HepMC3::GenEvent genEvt;
    //reader.read_event(genEvt);
    drInterface->set("DRsim","DRsimEventData");
    recoInterface->set("Reco","RecoEventData");
    unsigned int entries = recoInterface->entries();
    if(reader.m_tree->GetEntries()==drInterface->entries() && drInterface->entries()==recoInterface->entries()){
      readcount+=1;
      pass=1;
      //std::cout<<inname.Data()<<" "<<reader.m_tree->GetEntries()<<" "<<drInterface->entries()<<" "<<recoInterface->entries()<<std::endl;
    }
    else{
      std::cout<<inname.Data()<<" "<<reader.m_tree->GetEntries()<<" "<<drInterface->entries()<<" "<<recoInterface->entries()<<std::endl;
      }
    delete drInterface;
    delete recoInterface;
        }
        else{
          printf("%s no tree: ",inname.Data());
          for(int i=0;i<keys->GetEntries();i++){
            printf("%s ",keys->At(i)->GetName());

          }
          printf("\n");
        }
      }
      if(isjet==0){
        if(strcmp(keys->At(0)->GetName(),"DRsim")!=0 || strcmp(keys->At(1)->GetName(),"Reco")!=0){
          printf("noevent2\n");
        }
      }
    }
    else{
      printf("nokey1\n");
    }
    box->Close();
    //printf("load %s...;\n",inname.Data());
    //mychain.Add(inname+"/DRsim");
    if(pass!=1){
      writeFile.write(std::string(inname.Data()).c_str(), std::string(inname.Data()).size());
      writeFile.write("\n", 1);
    }
      
  }
  readFile.close();
  writeFile.close();

printf("--done\n");
printf("%d files %d events filled in %s\n",readcount,count,outname.Data());

return 0;
}
