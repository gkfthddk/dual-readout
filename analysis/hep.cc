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

int toroot(std::shared_ptr<HepMC3::GenParticle> cand){
  auto b = cand->status();
  //printf("root %d\n",cand->pid());
  if(cand->parents().size()==0){
    printf("root 1 %d %d\n",cand->pid(),cand->status());
    return 0;
  }
  if(cand->parents().at(0)->status()==-1){
    printf("root 2 %d %d\n",cand->pid(),cand->status());
    return 0;
  }
  printf("root 3 %d %d\n",cand->pid(),cand->status());
  for(int i =0;i < cand->parents().size();i++){
    toroot(cand->parents().at(i));
  }
  return 1;
}

int main(int argc, char* argv[]) {
  TString innameform = argv[1];
  TString eventnameform = argv[2];
  TString outname = argv[3];
  int isjet = atoi(argv[4]);
  int getroot=0;
  if (argc > 5) getroot = std::stoi(argv[6]);
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
  HepMC3::GenEvent genEvt;
  int pid_count[100] {0};
  int pids[100] {0};
  int pid_idx=-1;
  for(int i = 0; i<hep_entries; i++){
    hepreader.read_event(genEvt);
    for (auto ptc : genEvt.particles()) {
      int abspid = std::abs(ptc->pid());
      int stat = ptc->status();
      //printf("stat %d pid %d child %d\n",stat,abspid,ptc->children().size());
      if(stat==11){
        for (auto cptc : ptc->children()) {
          //printf("child stat %d pid %d\n",cptc->status(),cptc->pid());
          abspid = std::abs(cptc->pid());
          stat = cptc->status();
          pid_idx=-1;
          for(int j = 0; j<100; j++){
            //printf("!3!");
            if(pids[j]==0){
              pids[j]=abspid;
            }
            if(pids[j]==abspid){
              pid_idx=j;
              break;
            }
          }
          pid_count[pid_idx]+=1;
          //if(i<100)printf("ptc %d %d %d\n",i,pid_idx,pids[pid_idx]);

        //printf("i4i");
        if(getroot==1){
          if(i==9){
            printf("ent %d pid %d e %g stat %d\n",i,ptc->pid(),ptc->momentum().e(),ptc->status());
            //if(abspid<=16 and abspid>=11)continue;
            //qq[abspid]+=1;
            toroot(ptc);
          }
        }
        }
      }
    }
  }
  //printf("endhep\n");
  int main_pid=0;
  int count=0;
  int event_pid_count[100] {0};
  int event_pids[100] {0};
  for(int i = 0; i<100; i++){
    if(pids[i]!=0){
      event_pids[i]=pids[i];
      if(pid_count[i]>count){
        count=pid_count[i];
        main_pid=pids[i];
      }
    }
  }
  recoInterface->close();
  drInterface->close();
  hepreader.close();
  Int_t num_entry=-1;
  float E_C_img=0.;
  float E_Scorr_img=0.;
  std::vector<Int_t> * prt_pid=0;

  TFile infile(eventnameform.Data(), "read");
  TList* keys2;
  keys2 = infile.GetListOfKeys();
  //if(keys2->GetEntries()>1){
  //    for(int i=0;i<keys2->GetEntries();i++){
  //      printf("%d %s \n",i,keys2->At(i)->GetName());
  //    }
  //}
  //printf("~1~\n");
  TTree *eventtree = (TTree *) infile.Get("event");
  eventtree->SetBranchAddress("prt_pid",&prt_pid);
  eventtree->SetBranchAddress("E_C_img",&E_C_img);
  eventtree->SetBranchAddress("E_Scorr_img",&E_Scorr_img);

  int event_entries=eventtree->GetEntries();
  float ecimg[int(event_entries)]={0};
  float esimg[int(event_entries)]={0};
  int same_ec=0;
  int same_es=0;
  //printf("2 %d\n",event_entries);
  int ccount=0;
  for(int i = 0; i<event_entries; i++){
    eventtree->GetEntry(i);
    //printf("prt %d\n",prt_pid->size());
    for(int j = 0; j<prt_pid->size(); j++){
      int abspid=abs(prt_pid->at(j));
    //printf("prt %d %d\n",prt_pid->size(),abspid);
      pid_idx=-1;
      for(int k = 0; k<100; k++){
        if(event_pids[k]==0){
          event_pids[k]=abspid;
        }
        if(event_pids[k]==abspid){
          pid_idx=k;
          break;
        }
      }
      event_pid_count[pid_idx]+=1;
      //if(i<100)printf("%d %d %d %d\n",i,j,abspid,event_pid_count[pid_idx]);

    }
    /*ecimg[i]=E_C_img;
    esimg[i]=E_Scorr_img;
    for(int k = 0; k<i; k++){
      if(ecimg[k]==E_C_img){
        same_ec+=1;
      }

      if(esimg[k]==E_Scorr_img){
        same_es+=1;
      }
    }*/
  }
  //printf(". 3 \n");
  infile.Close();
  printf("same ec %d same es %d\n",same_ec,same_es);
  printf("\ncheck hep %d dr %d reco %d\nentries hep %d dr %d reco %d event %d\n",check_hep,check_dr,check_reco,hep_entries,dr_entries,reco_entries,event_entries);

  //printf("4\n");
  std::ofstream writeFile;
  writeFile.open(outname.Data());


  char buffer[100];
  int buf_size=0;
  buf_size=sprintf(buffer,"\ncheck hep %d dr %d reco %d\nentries hep %d dr %d reco %d event %d same_ec %d same_es %d\n",check_hep,check_dr,check_reco,hep_entries,dr_entries,reco_entries,event_entries,same_ec,same_es);
      writeFile.write(std::string(inname.Data()).c_str(), std::string(inname.Data()).size());
      writeFile.write(buffer, buf_size);
  for(int i = 0; i<100; i++){
    if(event_pids[i]!=0 || pids[i]!=0){
        buf_size=sprintf(buffer,"pid %d hep %d event %d %d\n",event_pids[i],pids[i],event_pid_count[i],i);
        //printf("pid %d hep %d event %d %d\n",event_pids[i],pids[i],event_pid_count[i],i);
        writeFile.write(buffer, buf_size);
    }
  }
writeFile.close();
}
