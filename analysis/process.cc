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
namespace fs = std::filesystem;
//./process ../../box/pi0_0 ../../test_0.root 0
//./process ../../elenshower ../../elenshower.root 0
//./process ../../box/enshowers/pi_en_%d.root 1200 ../../pienshower.root 0
int main(int , char* argv[]){
  TString innameform = argv[1];
  int num_file = std::stoi(argv[2]);
  TString outname = argv[3];  
  int isjet = atoi(argv[4]);
  TString inname;
  TFile* box;
  TList* keys;
  int readcount=0;
  for(int file_num=0; file_num<num_file;file_num++){
    inname.Form(innameform.Data(),file_num);
    if(!fs::exists(inname.Data())){
      printf("nofile %s;\n",inname.Data());
      continue;}
    box=TFile::Open(inname.Data(),"read");
    keys = box->GetListOfKeys();
    printf("file loaded %s..%d.;\n",inname.Data(),keys->GetEntries());
    if(keys->GetEntries()>1){
      /*if(isjet==1){
        if(strcmp(keys->At(1)->GetName(),"DRsim")!=0 || strcmp(keys->At(2)->GetName(),"Reco")!=0){
          printf("noevent1\n");
          continue;
        }
      }
      if(isjet==0){
        if(strcmp(keys->At(0)->GetName(),"DRsim")!=0 || strcmp(keys->At(1)->GetName(),"Reco")!=0){
          printf("noevent2\n");
          continue;
        }
      }*/
    }
    else{
      printf("nokey1 %d\n",keys->GetEntries());
      continue;
    }
    keys->Clear();
    box->Close();
  }
  new GeoSvc({"./bin/compact/DRcalo.xml"});

  auto m_geoSvc = GeoSvc::GetInstance();
  std::string m_readoutName = "DRcaloSiPMreadout";

  auto lcdd = m_geoSvc->lcdd();
  auto allReadouts = lcdd->readouts();
  if (allReadouts.find(m_readoutName) == allReadouts.end()) {
    throw std::runtime_error("Readout " + m_readoutName + " not found! Please check tool configuration.");
  } else {
    std::cout << "Reading EDM from the collection " << m_readoutName << std::endl;
  }

      auto segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDRcalo*>(m_geoSvc->lcdd()->readout(m_readoutName).segmentation().segmentation());
  
  
  
  
  
  float pi=TMath::Pi();
  float x=0.;
  float y=0.;
  float z=0.;
TFile outfile(outname.Data(), "recreate");
  std::vector<Int_t> tower_idx_eta;
  std::vector<Int_t> tower_idx_phi;
  std::vector<Float_t> tower_eta;
  std::vector<Float_t> tower_phi;
  std::vector<Int_t> tower_diff_eta;
  std::vector<Int_t> tower_diff_phi;
  std::vector<Int_t> tower_numx;
  std::vector<Int_t> tower_numy;
  float ptd=0.;
  float major_axis=0.;
  float minor_axis=0.;
  float deleta2=0.;
  float delphi2=0.;
  float width_Gen=0.;
  int cmult=0.;
  int mult=0.;
  int nmult=0.;
  int chad_mult=0.;
  int nhad_mult=0.;
  int electron_mult=0.;
  int muon_mult=0.;
  int photon_mult=0.;
  float mass_Gen=0.;
  float E_Gen=0.;
  float pt_Gen=0.;
  float E_C=0.;
  float E_S=0.;
  float E_Scorr=0.;
  float E_C_img=0.;
  float E_Scorr_img=0.;
  int n_C=0.;
  int n_S=0.;
  float E_DR=0.;
  float E_DRcorr=0.;
  float Eleak_nu=0;
  float Pleak =0;
  float phicut[46] = {-0.499591, -0.47728, -0.455022, -0.432663, -0.410903, -0.388499, -0.366236, -0.34418, -0.321856, -0.299631, -0.2775, -0.25534, -0.233038, -0.210916, -0.188759, -0.166583, -0.144445, -0.122228, -0.098836, -0.077592, -0.055648, -0.033304, -0.011214, 0.011112, 0.033254, 0.055531, 0.07785, 0.099895, 0.122083, 0.144345, 0.166667, 0.188781, 0.211015, 0.233232, 0.255433, 0.277524, 0.299768, 0.322163, 0.344244, 0.366127, 0.388528, 0.410884, 0.433001, 0.455096, 0.477211, 0.499555};
  float etacut[47] = {-0.509844, -0.487655, -0.465542, -0.443474, -0.421363, -0.399169, -0.376835, -0.35499, -0.332782, -0.310622, -0.288396, -0.266306, -0.244095, -0.221847, -0.198638, -0.17761, -0.155555, -0.133243, -0.110986, -0.088896, -0.066692, -0.044368, -0.022242, 1.8e-05, 0.022211, 0.04439, 0.066663, 0.088883, 0.111086, 0.133137, 0.155502, 0.17757, 0.199889, 0.222027, 0.244179, 0.266126, 0.288726, 0.310622, 0.332651, 0.354897, 0.377076, 0.399316, 0.421363, 0.443408, 0.465533, 0.48769, 0.509791};
  Float_t center_eta_img=0.;
  Float_t center_phi_img=0.;
  Float_t center_eta_gen=0.;
  Float_t center_phi_gen=0.;
  Float_t cen_Ecorr_eta=0.;
  Float_t cen_Ecorr_phi=0.;
  int num_phicut=46;
  int num_etacut=47;
  std::vector<Float_t> tower_e_s;
  std::vector<Float_t> tower_e_c;
  std::vector<Int_t> tower_n_s;
  std::vector<Int_t> tower_n_c;
  std::vector<Float_t> tower_ecor_s;
  std::vector<Float_t> tower_ecor_dr;
  std::vector<Float_t> fiber_e;
  std::vector<Float_t> fiber_ecor;
  std::vector<Float_t> fiber_ecor_s;
  std::vector<Float_t> fiber_ecor_c;
  std::vector<Int_t> fiber_n;
  std::vector<Int_t> fiber_itower;
  std::vector<Int_t> fiber_ix;
  std::vector<Int_t> fiber_iy;
  std::vector<Float_t> fiber_t;
  std::vector<Float_t> fiber_x;
  std::vector<Float_t> fiber_y;
  std::vector<Float_t> fiber_z;
  std::vector<Float_t> fiber_depth;
  std::vector<Float_t> fiber_r;
  std::vector<Float_t> fiber_theta;
  std::vector<Float_t> fiber_phi;
  std::vector<Float_t> fiber_eta;
  std::vector<Bool_t> fiber_iscerenkov;
  std::vector<Float_t> ptc_E;
  std::vector<Float_t> ptc_pt;
  std::vector<Float_t> ptc_theta;
  std::vector<Float_t> ptc_phi;
  std::vector<Float_t> ptc_eta;

  std::vector<Int_t> ptc_pid;

  Int_t num_tower=0;
  TTree eventtree("event","event info");
  //eventtree.SetAutoSave(0);
  eventtree.Branch("tower_eta","vector<Float_t>",&tower_eta);
  eventtree.Branch("tower_phi","vector<Float_t>",&tower_phi);
  eventtree.Branch("tower_diff_eta","vector<Int_t>",&tower_diff_eta);
  eventtree.Branch("tower_diff_phi","vector<Int_t>",&tower_diff_phi);
  eventtree.Branch("tower_idx_eta","vector<Int_t>",&tower_idx_eta);
  eventtree.Branch("tower_idx_phi","vector<Int_t>",&tower_idx_phi);
  eventtree.Branch("tower_numx","vector<Int_t>",&tower_numx);
  eventtree.Branch("tower_numy","vector<Int_t>",&tower_numy);
  eventtree.Branch("ptd",&ptd,"ptd/F");
  eventtree.Branch("major_axis",&major_axis,"major_axis/F");
  eventtree.Branch("minor_axis",&minor_axis,"minor_axis/F");
  eventtree.Branch("center_eta_gen",&center_eta_gen,"center_eta_gen/F");
  eventtree.Branch("center_phi_gen",&center_phi_gen,"center_phi_gen/F");
  eventtree.Branch("center_eta_img",&center_eta_img,"center_eta_img/F");
  eventtree.Branch("center_phi_img",&center_phi_img,"center_phi_img/F");
  eventtree.Branch("deleta2",&deleta2,"deleta2/F");
  eventtree.Branch("delphi2",&delphi2,"delphi2/F");
  eventtree.Branch("width_Gen",&width_Gen,"width_Gen/F");
  eventtree.Branch("mult",&mult,"mult/I");
  eventtree.Branch("cmult",&cmult,"cmult/I");
  eventtree.Branch("nmult",&nmult,"nmult/I");
  eventtree.Branch("chad_mult",&chad_mult,"chad_mult/I");
  eventtree.Branch("nhad_mult",&nhad_mult,"nhad_mult/I");
  eventtree.Branch("electron_mult",&electron_mult,"electron_mult/I");
  eventtree.Branch("muon_mult",&muon_mult,"muon_mult/I");
  eventtree.Branch("photon_mult",&photon_mult,"photon_mult/I");
  eventtree.Branch("E_Gen",&E_Gen,"E_Gen/F");
  eventtree.Branch("mass_Gen",&mass_Gen,"mass_Gen/F");
  eventtree.Branch("pt_Gen",&pt_Gen,"pt_Gen/F");
  eventtree.Branch("E_C",&E_C,"E_C/F");
  eventtree.Branch("E_C_img",&E_C,"E_C_img/F");
  eventtree.Branch("E_S",&E_S,"E_S/F");
  eventtree.Branch("E_Scorr",&E_Scorr,"E_Scorr/F");
  eventtree.Branch("E_Scorr_img",&E_Scorr_img,"E_Scorr_img/F");
  eventtree.Branch("n_C",&n_C,"n_C/I");
  eventtree.Branch("n_S",&n_S,"n_S/I");
  eventtree.Branch("E_DR",&E_DR,"E_DR/F");
  eventtree.Branch("E_DRcorr",&E_DRcorr,"E_DRcorr/F");
  eventtree.Branch("Eleak_nu",&Eleak_nu,"Eleak_nu/F");
  eventtree.Branch("Pleak",&Pleak,"Pleak/F");
  eventtree.Branch("tower_e_s","vector<Float_t>",&tower_e_s);
  eventtree.Branch("tower_e_c","vector<Float_t>",&tower_e_c);
  eventtree.Branch("tower_n_s","vector<Int_t>",&tower_n_s);
  eventtree.Branch("tower_n_c","vector<Int_t>",&tower_n_c);
  eventtree.Branch("tower_ecor_s","vector<Float_t>",&tower_ecor_s);
  eventtree.Branch("tower_ecor_dr","vector<Float_t>",&tower_ecor_dr);
  eventtree.Branch("fiber_e","vector<Float_t>",&fiber_e);
  eventtree.Branch("fiber_ecor","vector<Float_t>",&fiber_ecor);
  eventtree.Branch("fiber_ecor_s","vector<Float_t>",&fiber_ecor_s);
  eventtree.Branch("fiber_ecor_c","vector<Float_t>",&fiber_ecor_c);
  eventtree.Branch("fiber_n","vector<Int_t>",&fiber_n);
  eventtree.Branch("fiber_itower","vector<Int_t>",&fiber_itower);
  eventtree.Branch("fiber_ix","vector<Int_t>",&fiber_ix);
  eventtree.Branch("fiber_iy","vector<Int_t>",&fiber_iy);
  eventtree.Branch("fiber_t","vector<Float_t>",&fiber_t);
  eventtree.Branch("fiber_x","vector<Float_t>",&fiber_x);
  eventtree.Branch("fiber_y","vector<Float_t>",&fiber_y);
  eventtree.Branch("fiber_z","vector<Float_t>",&fiber_z);
  eventtree.Branch("fiber_depth","vector<Float_t>",&fiber_depth);
  eventtree.Branch("fiber_r","vector<Float_t>",&fiber_r);
  eventtree.Branch("fiber_theta","vector<Float_t>",&fiber_theta);
  eventtree.Branch("fiber_phi","vector<Float_t>",&fiber_phi);
  eventtree.Branch("fiber_eta","vector<Float_t>",&fiber_eta);
  eventtree.Branch("fiber_iscerenkov","vector<Bool_t>",&fiber_iscerenkov);
  eventtree.Branch("cen_Ecorr_eta",&cen_Ecorr_eta,"cen_Ecorr_eta/F");
  eventtree.Branch("cen_Ecorr_phi",&cen_Ecorr_phi,"cen_Ecorr_phi/F");
  eventtree.Branch("ptc_E","vector<Float_t>",&ptc_E);
  eventtree.Branch("ptc_pt","vector<Float_t>",&ptc_pt);
  eventtree.Branch("ptc_theta","vector<Float_t>",&ptc_theta);
  eventtree.Branch("ptc_phi","vector<Float_t>",&ptc_phi);
  eventtree.Branch("ptc_eta","vector<Float_t>",&ptc_eta);
  eventtree.Branch("ptc_pid","vector<Int_t>",&ptc_pid);

  //Float_t voxel_ecor_s_[729000];//90*90*90
  //Int_t voxel_n_s_[729000];
  Float_t image_ecor_c_[28224];//168*168
  Int_t image_n_c_[28224];
  Float_t image_ecor_s_[28224];//168*168
  Int_t image_n_s_[28224];
  Float_t point_2048_[8192];//90*90
  #define BRANCH_A_(name, size, suffix) eventtree.Branch(#name, & name##_, #name"["#size"]/"#suffix);
  #define BRANCH_AF(name, size)  BRANCH_A_(name, size, F);
  #define BRANCH_AI(name, size)  BRANCH_A_(name, size, I);
  #define FILL_ZERO(array, size) std::fill(array, array + size, 0.0);
  #define FILL_ZEROI(array, size) std::fill(array, array + size, 0);
  //BRANCH_AF(voxel_ecor_s, 729000);
  //BRANCH_AI(voxel_n_s, 729000);
  BRANCH_AF(image_ecor_c, 28224);
  BRANCH_AI(image_n_c, 28224);//don't forget to add FILL_ZERO in loop
  BRANCH_AF(image_ecor_s, 28224);
  BRANCH_AI(image_n_s, 28224);//don't forget to add FILL_ZERO in loop
  BRANCH_AI(point_2048, 8192);//don't forget to add FILL_ZERO in loop
  //FILL_ZERO(voxel_ecor_s_,729000);
  //FILL_ZEROI(voxel_n_s_,729000);
  FILL_ZERO(image_ecor_c_,28224);
  FILL_ZEROI(image_n_c_,28224);
  FILL_ZERO(image_ecor_s_,28224);
  FILL_ZEROI(image_n_s_,28224);
  FILL_ZEROI(point_2048_,8192);
  int depthbin=89;//22+1
  float depthmin=34;
  float depthmax=3060;
  int phibin=56;
  int etabin=56;

  //float phimin=-0.008;
  //float phimax=0.0081;
  //float etamin=0.003;
  //float etamax=0.0191;
  float phimin=-0.0072;
  float phimax=0.0078;
  float etamin=0.004;
  float etamax=0.019;
  if(isjet==1){
    printf("jet %d\n",isjet);
    //int n_in_bin=0.25; //number of towers in a bin
    phibin=90;
    etabin=90;
    phimin=phicut[0];
    phimax=phicut[num_phicut-1];
    etamin=etacut[0];
    etamax=etacut[num_etacut-1];
    //etamin90=-0.0221*90/4.;
    //etamax90=0.0221*90/4.;
  }
  float depthsize=1.*(depthmax-depthmin)/depthbin;
  float phisize=1.*(phimax-phimin)/phibin;
  float etasize=1.*(etamax-etamin)/etabin;
  float eta=0.;
  float phi=0.;
  Float_t w2=0.;
  cen_Ecorr_eta=0.;
  cen_Ecorr_phi=0.;
  Float_t pt_square=0.;
  Float_t m00=0.;
  Float_t m01=0.;
  Float_t m11=0.;
  float depth=0.;
  int depthindex=-1;
  int phiindex=-1;
  int etaindex=-1;
  int voxindex=-1;
  int imgindex=-1;
  int buf_index=0;
  bool isopposite = false;

  TVector3* fiberxyz = new TVector3();
  TVector3* towerxyz = new TVector3();
  TLorentzVector* ptc_p = new TLorentzVector();
  int count=0;
  printf("file will be checked %d\n",num_file);
  for(int file_num=0; file_num<num_file;file_num++){
    inname.Form(innameform.Data(),file_num);
    if(!fs::exists(inname.Data())){
      printf("nofile %s;\n",inname.Data());
      continue;}
    box=TFile::Open(inname.Data(),"read");
    keys = box->GetListOfKeys();
    
    if(keys->GetEntries()>1){
      /*if(isjet==1){
        if(strcmp(keys->At(1)->GetName(),"DRsim")!=0 || strcmp(keys->At(2)->GetName(),"Reco")!=0){
          printf("noevent1\n");
          continue;
        }
      }
      if(isjet==0){
        if(strcmp(keys->At(0)->GetName(),"DRsim")!=0 || strcmp(keys->At(1)->GetName(),"Reco")!=0){
          printf("noevent2\n");
          continue;
        }
      }*/
    }
    else{
      printf("nokey1 %d\n",keys->GetEntries());
      continue;
    }
    keys->Clear();
    box->Close();
    printf("loaded %s...;\n",inname.Data());
    //mychain.Add(inname+"/DRsim");
    RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(inname.Data()));
    RootInterface<RecoInterface::RecoEventData>* recoInterface = new RootInterface<RecoInterface::RecoEventData>(std::string(inname.Data()));
    drInterface->set("DRsim","DRsimEventData");
    recoInterface->set("Reco","RecoEventData");
    readcount+=1;
    unsigned int entries = recoInterface->entries();
    
    while (recoInterface->numEvt() < entries) {
    //printf("%d%%--",int(count/entries));
    //if(count%int(entries/10.)==0)printf("%d%%--\n",int(100.*count/entries));
    DRsimInterface::DRsimEventData drEvt;
    RecoInterface::RecoEventData recoEvt;
    drInterface->read(drEvt);
    recoInterface->read(recoEvt);
    E_C=recoEvt.E_C;
    E_S=recoEvt.E_S;
    E_Scorr=recoEvt.E_Scorr;
    n_C=recoEvt.n_C;
    n_S=recoEvt.n_S;
    E_DR=recoEvt.E_DR;
    E_DRcorr=recoEvt.E_DRcorr;
    Eleak_nu=0;
    Pleak =0;
    tower_eta.clear();
    tower_phi.clear();
    tower_idx_eta.clear();
    tower_idx_phi.clear();
    tower_numx.clear();
    tower_numy.clear();
    tower_e_s.clear();
    tower_e_c.clear();
    tower_ecor_s.clear();
    tower_ecor_dr.clear();
    tower_n_s.clear();
    tower_n_c.clear();
     for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      if ( std::abs(leak.pdgId)==12 || std::abs(leak.pdgId)==14 || std::abs(leak.pdgId)==16 ) {
        Eleak_nu += leak4vec.P();
      } else {
        Pleak += leak4vec.P();
      }
      }

    int fibercheck=0;
    for ( int num_jet=0; num_jet<2; num_jet++){
      int fibercount=0;
      fiber_e.clear();
      fiber_ecor.clear();
      fiber_ecor_s.clear();
      fiber_ecor_c.clear();
      fiber_n.clear();
      fiber_itower.clear();
      fiber_ix.clear();
      fiber_iy.clear();
      fiber_t.clear();
      fiber_x.clear();
      fiber_y.clear();
      fiber_z.clear();
      fiber_depth.clear();
      fiber_r.clear();
      fiber_theta.clear();
      fiber_phi.clear();
      fiber_eta.clear();
      fiber_iscerenkov.clear();
      ptc_E.clear();
      ptc_pt.clear();
      ptc_theta.clear();
      ptc_phi.clear();
      ptc_eta.clear();
      ptc_pid.clear();
      tower_diff_eta.clear();
      tower_diff_phi.clear();
      //FILL_ZERO(voxel_ecor_s_,729000);
      //FILL_ZEROI(voxel_n_s_,729000);
      FILL_ZERO(image_ecor_c_,28224);
      FILL_ZEROI(image_n_c_,28224);
      FILL_ZERO(image_ecor_s_,28224);
      FILL_ZEROI(image_n_s_,28224);
      if(num_jet==1 && isjet==0) break;
      m00 = 0.;
      m01 = 0.;
      m11 = 0.;
      pt_square = 0.;
      pt_Gen = 0.;
      cmult=0;
      nmult=0;
      chad_mult=0;
      nhad_mult=0;
      electron_mult=0;
      muon_mult=0;
      photon_mult=0;
      mult=0;
      E_Gen=0.;
      deleta2=0.;
      delphi2=0.;
      width_Gen=0.;
      center_eta_gen=0.;
      center_phi_gen=0.;
      cen_Ecorr_eta=0.;
      cen_Ecorr_phi=0.;
      TLorentzVector ptc_sum;
      
      for (auto genptc : drEvt.GenPtcs){
        ptc_p->SetPxPyPzE(genptc.px,genptc.py,genptc.pz,genptc.E);
        phi=float(ptc_p->Phi());
        if(isjet==1 && num_jet==0 && abs(phi)>pi/2.)continue;
        if(isjet==1 && num_jet==1){
          if(abs(phi)>pi/2.){
            isopposite=true;
            if(phi<0){
              phi=pi+phi;
            }
            else{
              phi=phi-pi;
            }
          }
          else continue;
        }
        E_Gen+=genptc.E;
        center_eta_gen+=ptc_p->Eta()*ptc_p->Pt();
        center_phi_gen+=phi*ptc_p->Pt();
        pt_Gen+=ptc_p->Pt();
        ptc_sum=ptc_sum+*ptc_p;
      }
      //printf("mass %g\n",ptc_sum.M());
      mass_Gen=ptc_sum.M();
      center_eta_gen=center_eta_gen/pt_Gen;
      center_phi_gen=center_phi_gen/pt_Gen;
      for (auto genptc : drEvt.GenPtcs){
        ptc_p->SetPxPyPzE(genptc.px,genptc.py,genptc.pz,genptc.E);
        phi=float(ptc_p->Phi());
        if(isjet==1 && num_jet==0 && abs(phi)>pi/2.)continue;
        if(isjet==1 && num_jet==1){
          if(abs(phi)>pi/2.){
            isopposite=true;
            if(phi<0){
              phi=pi+phi;
            }
            else{
              phi=phi-pi;
            }
          }
          else continue;
        }
        deleta2+=std::pow(ptc_p->Eta()-center_eta_gen,2)*ptc_p->Pt();
        delphi2+=std::pow(phi-center_phi_gen,2)*ptc_p->Pt();
        ptc_theta.push_back(ptc_p->Theta());
        ptc_eta.push_back(ptc_p->Eta());
        ptc_phi.push_back(phi);
        ptc_pid.push_back(genptc.pdgId);
        mult++;
        if(ptc_p->Pt()>500.){
        if(abs(genptc.pdgId)==11 || abs(genptc.pdgId)==13 || abs(genptc.pdgId)==321 || abs(genptc.pdgId)==211 || abs(genptc.pdgId)==2212){
          cmult++;
          if(abs(genptc.pdgId)==11){
            electron_mult++;
          }
          else if(abs(genptc.pdgId)==13){
            muon_mult++;
          }
          else{
            chad_mult++;
          }
        }
        else if(abs(genptc.pdgId)==12 || abs(genptc.pdgId)==14 || abs(genptc.pdgId)==16){
          
        }
        else if(genptc.pdgId!=0){
          nmult++;
          if(genptc.pdgId==22){
            photon_mult++;
          }
          else{
            nhad_mult++;
          }
        }
        }
        ptc_E.push_back(genptc.E);
        ptc_pt.push_back(ptc_p->Pt());
        x = ptc_p->Eta();
        y = phi;
        w2 = std::pow(ptc_p->Pt(), 2);
        m00 += w2 * std::pow(x, 2);
        m01 -= w2 * x * y;
        m11 += w2 * std::pow(y, 2);
        pt_square+=w2;

      }
      
      deleta2=deleta2/pt_Gen;
      delphi2=delphi2/pt_Gen;
      width_Gen=delphi2+deleta2;
      TMatrixFSym covariance_matrix(2);
      covariance_matrix(0, 0) = m00;
      covariance_matrix(0, 1) = m01;
      covariance_matrix(1, 1) = m11;
      TVectorF eigen_values;
      covariance_matrix.EigenVectors(eigen_values);
      major_axis = std::sqrt(eigen_values[0] / pt_square);
      minor_axis = std::sqrt(eigen_values[1] / pt_square);
      ptd=std::sqrt(pt_square)/pt_Gen;
      num_tower = -1;
      TVector3* towercenterxyz = new TVector3();//should be energy center
      
      auto towercenterpos = segmentation->towerposition(0,0);
      towercenterxyz->SetXYZ(towercenterpos.x(),towercenterpos.y(),towercenterpos.z());
      auto towercenterphi=float(towercenterxyz->Phi());
      auto towercentereta=float(towercenterxyz->Eta());
      center_eta_img=0.;
      center_phi_img=0.;
      E_Scorr_img=0.;
      E_C_img=0.;
      for (auto tower : recoEvt.towers) {

        tower_e_s.push_back(float(tower.E_S));
        tower_e_c.push_back(float(tower.E_C));
        tower_ecor_s.push_back(tower.E_Scorr);
        tower_ecor_dr.push_back(tower.E_DRcorr);
        tower_n_s.push_back(float(tower.n_S));
        tower_n_c.push_back(float(tower.n_C));
        tower_numx.push_back(tower.numx);
        tower_numy.push_back(tower.numy);

        //towerData.iTheta = segmentation->numEta(hit->GetSiPMnum());
        //towerData.iPhi = segmentation->numPhi(hit->GetSiPMnum());
        int noEta = tower.iTheta;
        int noPhi = tower.iPhi;
        /*DRparamBase* paramBase = nullptr;
        if ( fParamEndcap->unsignedTowerNo(noEta) >= fParamBarrel->GetTotTowerNum() ) paramBase = segmentation.fParamEndcap;
        else paramBase = segmentation.fParamBarrel;
        // This should not be called while building detector geometry
        if (!paramBase->IsFinalized()) throw std::runtime_error("GridDRcalo::position should not be called while building detector geometry!");
        paramBase->SetDeltaThetaByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
        paramBase->SetThetaOfCenterByTowerNo(noEta, fParamBarrel->GetTotTowerNum());
        paramBase->SetIsRHSByTowerNo(noEta);
        paramBase->init();
        auto towerpos = paramBase->GetTowerPos(noPhi);*/
        auto towerpos = segmentation->towerposition(noEta,noPhi);
        towerxyz->SetXYZ(towerpos.x(),towerpos.y(),towerpos.z());
        auto towerphi=float(towerxyz->Phi());
        auto towereta=float(towerxyz->Eta());
        tower_phi.push_back(towerphi);
        tower_eta.push_back(towereta);
        int diff_phi=int(1 +TMath::Nint((towerphi-towercenterphi)/0.022));
        int diff_eta=int(1 +TMath::Nint((towereta-towercentereta)/0.022));
        tower_diff_eta.push_back(diff_eta);
        tower_diff_phi.push_back(diff_phi);
        num_tower+=1;
        int checktower=0;
        
        
        for (auto fiber : tower.fibers) {
          
            auto pos = segmentation->position(fiber.fiberNum);
            isopposite=false;
            x=float(pos.x());
            y=float(pos.y());
            z=float(pos.z());
            fiberxyz->SetXYZ(x,y,z);
            phi=float(fiberxyz->Phi());
            eta=float(fiberxyz->Eta());
            if(isjet==1 && num_jet==0 && abs(phi)>pi/2.)break;
            if(isjet==1 && num_jet==1){
              if(abs(phi)>pi/2.){
                isopposite=true;
                if(phi<0){
                  phi=pi+phi;
                }
                else{
                  phi=phi-pi;
                }
              }
              else break;
            }
            tower_idx_eta.push_back(segmentation->numEta(fiber.fiberNum));
            tower_idx_phi.push_back(segmentation->numPhi(fiber.fiberNum));
            checktower=1;
            fibercheck=1;
            fibercount+=1;
            depthindex=-1;
            phiindex=-1;
            etaindex=-1;
            voxindex=-1;
            imgindex=-1;
            fiber_x.push_back(x);
            fiber_y.push_back(y);
            fiber_z.push_back(z);
            fiber_iscerenkov.push_back(segmentation->IsCerenkov(fiber.fiberNum));
            fiber_e.push_back(float(fiber.E));
            fiber_ecor.push_back(float(fiber.Ecorr));
            if(segmentation->IsCerenkov(fiber.fiberNum)){
              fiber_ecor_s.push_back(float(0.));
              fiber_ecor_c.push_back(float(fiber.Ecorr));
            }
            else{
              fiber_ecor_s.push_back(float(fiber.Ecorr));
              fiber_ecor_c.push_back(float(0.));
            }

            fiber_n.push_back(fiber.n);
            fiber_t.push_back(fiber.t);
            fiber_itower.push_back(num_tower);
            fiber_ix.push_back(segmentation->x(fiber.fiberNum));
            fiber_iy.push_back(segmentation->y(fiber.fiberNum));
            fiber_depth.push_back(float(fiber.depth));
            fiber_r.push_back(float(TMath::Sqrt(x*x+y*y+z*z)));
            depth=fiber.depth;
            fiber_phi.push_back(phi);
            fiber_theta.push_back(float(fiberxyz->Theta()));
            fiber_eta.push_back(eta);
            //std::sort(phis.begin(),phis.end());
            
            if(isjet==1){// need to revise#############################
                if(phimax>phi && phimin<=phi){
                if(etamax>eta && etamin<=eta){
                
                for(buf_index=0;buf_index<num_phicut-1;buf_index++){
                  if(28>int(segmentation->x(fiber.fiberNum))){
                    if(phi<(phicut[buf_index]+phicut[buf_index+1])/2.){
                      phiindex=Int_t(buf_index*2);
                    }
                    else{
                      phiindex=Int_t(buf_index*2+1);
                    }
                    break;
                  }
                }
                for(buf_index=0;buf_index<num_etacut-1;buf_index++){
                  if(etacut[buf_index]<eta && eta<etacut[buf_index+1]){
                    if(eta<(etacut[buf_index]+etacut[buf_index+1])/2.){
                      if(buf_index>0)etaindex=Int_t(buf_index*2-1);
                    }
                    else{
                      if(buf_index<num_etacut-2)etaindex=Int_t(buf_index*2);
                    }
                    break;
                  }
                }
                /*for(buf_index=0;buf_index<num_phicut-1;buf_index++){
                  if(phicut[buf_index]<phi && phi<phicut[buf_index+1]){
                    if(phi<(phicut[buf_index]+phicut[buf_index+1])/2.){
                      phiindex=Int_t(buf_index*2);
                    }
                    else{
                      phiindex=Int_t(buf_index*2+1);
                    }
                    break;
                  }
                }
                for(buf_index=0;buf_index<num_etacut-1;buf_index++){
                  if(etacut[buf_index]<eta && eta<etacut[buf_index+1]){
                    if(eta<(etacut[buf_index]+etacut[buf_index+1])/2.){
                      if(buf_index>0)etaindex=Int_t(buf_index*2-1);
                    }
                    else{
                      if(buf_index<num_etacut-2)etaindex=Int_t(buf_index*2);
                    }
                    break;
                  }
                }*/
                }}
            }
            else{
              if(!segmentation->IsCerenkov(fiber.fiberNum)){
              cen_Ecorr_eta+=eta*fiber.Ecorr;
              cen_Ecorr_phi+=phi*fiber.Ecorr;
              }
                
                if(diff_eta<=2 && diff_eta>=0 && diff_phi<=2 && diff_phi>=0){
                  phiindex = Int_t(55-segmentation->x(fiber.fiberNum));
                  etaindex = Int_t(55-segmentation->y(fiber.fiberNum));
                  if(diff_eta<1){
                    phiindex = Int_t(segmentation->x(fiber.fiberNum));
                    etaindex = Int_t(segmentation->y(fiber.fiberNum));
                  }
                }
                else{
                  phiindex=-1;
                  etaindex=-1;
                }
            }
          
            
            if(phiindex!=-1 && etaindex!=-1){
                imgindex=3*phibin*etabin*diff_eta+3*phibin*etaindex+phibin*diff_phi+phiindex;//[eta,phi]
                //imgindex=phibin*etabin*diff_phi+phibin*etaindex+*phibin+phiindex;//[eta,phi]
                //imgindex=phibin*etaindex+phiindex;//[eta,phi]
                if(segmentation->IsCerenkov(fiber.fiberNum)){
                  image_ecor_c_[imgindex]+=fiber.Ecorr;
                  E_C_img+=fiber.Ecorr;
                  image_n_c_[imgindex]+=fiber.n;
                }
                else{
                  image_ecor_s_[imgindex]+=fiber.Ecorr;
                  E_Scorr_img+=fiber.Ecorr;
                  center_eta_img+=fiber.Ecorr*eta;
                  center_phi_img+=fiber.Ecorr*phi;
                  image_n_s_[imgindex]+=fiber.n;
                  if(depthmax>depth && depthmin<=depth){
                  depthindex=Int_t(1.*(depth-depthmin)/depthsize);
                  }
                  if(depth==0. || depthindex!=-1){
                    depthindex+=1;
                    voxindex=phibin*etabin*depthindex+imgindex;//[depth,eta,phi]
                    //voxel_ecor_s_[voxindex]+=fiber.Ecorr;
                    //voxel_n_s_[voxindex]+=fiber.n;
                  }
                }
            } 
        }

        
       }
      center_eta_img=center_eta_img/E_Scorr_img;
      center_phi_img=center_phi_img/E_Scorr_img;
      cen_Ecorr_eta=cen_Ecorr_eta/recoEvt.E_Scorr;
      cen_Ecorr_phi=cen_Ecorr_phi/recoEvt.E_Scorr;
       if(fibercheck==1){
         eventtree.Fill();
       }
       else{
         printf("no fiber\n");
       }
     }
     count+=1;
   }
    
    delete drInterface;
    delete recoInterface;

  }

//TFile outfile(outname.Data(), "recreate");
outfile.Write("");
//outfile.Write("",TObject::kOverwrite);
outfile.Close();
printf("--done\n");
printf("%d files %d events filled in %s\n",readcount,int(eventtree.GetEntries()),outname.Data());

return 0;
}
