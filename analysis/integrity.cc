#include "RootInterface.h"
#include "RecoInterface.h"
#include "DRsimInterface.h"

#include "fastjet/PseudoJet.hh"

#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <exception>
#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
  TString innameform = argv[1];
  TString outname = argv[2];
  TString inname;
  TString reconame;
  TList* keys;
  TFile* box;
    inname=innameform;
    reconame=inname.Copy();
    reconame.ReplaceAll(".root","_Reco.root");
    int check_hep=0;
    int check_dr=0;
    int check_reco=0;
    int pass=0;
    int code=0;
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
        }
      }
      else{
        //printf("%s no key\n",inname.Data());
        pass=1;
      }
      keys->Clear();
      box->Close();
      if(fs::exists(reconame.Data())){
          //printf("nofile %s;\n",reconame.Data());
          box=TFile::Open(reconame.Data(),"read");
          keys = box->GetListOfKeys();
          if(keys->GetEntries()>0){
            for(int i =0;i<keys->GetEntries();i++){
              if(strcmp(keys->At(i)->GetName(),"Reco")==0){
                check_reco=1;
              }
            }
          }
      }
    }
  HepMC3::ReaderRootTree hepreader(inname.Data());
  unsigned int hep_entries = hepreader.m_tree->GetEntries();
  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(inname.Data()));
  RootInterface<RecoInterface::RecoEventData>* recoInterface = new RootInterface<RecoInterface::RecoEventData>(std::string(reconame.Data()));
  drInterface->set("DRsim","DRsimEventData");
  recoInterface->set("Reco","RecoEventData");
  unsigned int dr_entries = drInterface->entries();
  unsigned int reco_entries = recoInterface->entries();
  recoInterface->close();
  drInterface->close();
  hepreader.close();
  std::ofstream writeFile;
  writeFile.open(outname.Data());
  char buffer[100];
  int buf_size=0;
  buf_size=sprintf(buffer,"\ncheck hep %d dr %d reco %d\nentries hep %d dr %d reco %d",check_hep,check_dr,check_reco,hep_entries,dr_entries,reco_entries);
      writeFile.write(std::string(inname.Data()).c_str(), std::string(inname.Data()).size());
      writeFile.write(buffer, buf_size);
writeFile.close();
}
