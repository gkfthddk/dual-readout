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


//./process ../../el2entest2_%d.root 0 1 ../../elentest2 0 0
//./process /hdfs/user/yulee/DRsim/simdata/cleanhad/${pid}_${energy}GeV_100_${i}_%d.root 0 1 $here/pack/${pid}_${energy}GeV_${i}.root 1 $here/hadbox/${pid}_${energy}GeV_100_${i}_%d_Reco.root
int main(int argc, char* argv[]){
  TString innameform = argv[1];
  int num_file_begin = std::stoi(argv[2]);
  int num_file_end = std::stoi(argv[3]);
  TString filename = argv[4];
  TString inname;
  TFile* box;
  TList* keys;
  int readcount=0;
  printf("file will be checked %d\n",num_file_end-num_file_begin);
  int check_hep=0;
  int nokey=0;
  for(int file_num=num_file_begin; file_num<num_file_end;file_num++){
    inname.Form(innameform.Data(),file_num);
    if(!fs::exists(inname.Data())){
      printf("nofile %s;\n",inname.Data());
      continue;}
    //printf("file loaded %s..%d.;\n",inname.Data(),keys->GetEntries());
    nokey=0;
    check_hep=0;
    box=TFile::Open(inname.Data(),"read");
    keys = box->GetListOfKeys();
    keys = box->GetListOfKeys();
    if(keys->GetEntries()>1){
      for(int i=0;i<keys->GetEntries();i++){
        if(strcmp(keys->At(i)->GetName(),"hepmc3_tree")==0)check_hep=1;
      }
    }
    else{
      nokey=1;
    }
    keys->Clear();
    box->Close();
    readcount+=1;
  }
  readcount=1;
  if(readcount==0){
    printf("no file found..\n");
    return(1);
  }
  readcount=0;
  float pi=TMath::Pi();
  float x=0.;
  float y=0.;
  float z=0.;
  TString outname=filename;
  TFile *outfile = new TFile(outname.Data(), "recreate");
  printf("outfile %s\n",outname.Data());
  float ptd=0.;
  float major_axis=0.;
  float minor_axis=0.;
  float deltheta2=0.;
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
  float p_Gen=0.;
  float E_C=0.;
  float E_S=0.;
  float Edep=0.;
  float E_Scorr=0.;
  int n_C=0.;
  int n_S=0.;
  int idx_phi=0.;
  int idx_theta=0.;
  float E_DR=0.;
  float E_DRcorr=0.;
  float Eleak_nu=0;
  float Pleak =0;
  int same_ec=0;
  int same_es=0;
  Float_t center_theta_gen=0.;
  Float_t center_phi_gen=0.;
  Float_t cen_Ecorr_theta=0.;
  Float_t cen_Ecorr_phi=0.;
  Int_t num_jet=0;
  Int_t num_entry=-1;
  Int_t not_phi=0;
  Int_t not_theta=0;
  std::vector<Float_t> prt_E;
  std::vector<Float_t> prt_pid;
  std::vector<Float_t> prt_phi;
  std::vector<Float_t> prt_theta;
  std::vector<Float_t> prt_px;
  std::vector<Float_t> prt_py;
  std::vector<Float_t> prt_pz;
  std::vector<Float_t> ptc_E;
  std::vector<Float_t> ptc_pt;
  std::vector<Float_t> ptc_p;
  std::vector<Float_t> ptc_theta;
  std::vector<Float_t> ptc_phi;
  std::vector<Float_t> ptc_eta;
  std::vector<Float_t> ptc_px;
  std::vector<Float_t> ptc_py;
  std::vector<Float_t> ptc_pz;
  std::vector<Float_t> ptc_vx;
  std::vector<Float_t> ptc_vy;
  std::vector<Float_t> ptc_vz;
  std::vector<Float_t> ptc_vt;
  std::vector<Int_t> ptc_pid;

  Int_t num_tower=0;
  TTree * eventtree = new TTree("event","event info");
  //eventtree->SetAutoSave(0);
  eventtree->Branch("ptd",&ptd,"ptd/F");
  eventtree->Branch("major_axis",&major_axis,"major_axis/F");
  eventtree->Branch("minor_axis",&minor_axis,"minor_axis/F");
  eventtree->Branch("center_theta_gen",&center_theta_gen,"center_theta_gen/F");
  eventtree->Branch("center_phi_gen",&center_phi_gen,"center_phi_gen/F");
  eventtree->Branch("deltheta2",&deltheta2,"deltheta2/F");
  eventtree->Branch("delphi2",&delphi2,"delphi2/F");
  eventtree->Branch("width_Gen",&width_Gen,"width_Gen/F");
  eventtree->Branch("mult",&mult,"mult/I");
  eventtree->Branch("cmult",&cmult,"cmult/I");
  eventtree->Branch("nmult",&nmult,"nmult/I");
  eventtree->Branch("chad_mult",&chad_mult,"chad_mult/I");
  eventtree->Branch("nhad_mult",&nhad_mult,"nhad_mult/I");
  eventtree->Branch("electron_mult",&electron_mult,"electron_mult/I");
  eventtree->Branch("muon_mult",&muon_mult,"muon_mult/I");
  eventtree->Branch("photon_mult",&photon_mult,"photon_mult/I");
  eventtree->Branch("num_jet",&num_jet,"num_jet/I");
  eventtree->Branch("num_entry",&num_entry,"num_entry/I");
  eventtree->Branch("not_theta",&not_theta,"not_theta/I");
  eventtree->Branch("not_phi",&not_phi,"not_phi/I");
  eventtree->Branch("E_Gen",&E_Gen,"E_Gen/F");
  eventtree->Branch("mass_Gen",&mass_Gen,"mass_Gen/F");
  eventtree->Branch("pt_Gen",&pt_Gen,"pt_Gen/F");
  eventtree->Branch("p_Gen",&p_Gen,"p_Gen/F");
  eventtree->Branch("Edep",&Edep,"Edep/F");
  eventtree->Branch("E_C",&E_C,"E_C/F");
  eventtree->Branch("E_S",&E_S,"E_S/F");
  eventtree->Branch("E_Scorr",&E_Scorr,"E_Scorr/F");
  eventtree->Branch("n_C",&n_C,"n_C/I");
  eventtree->Branch("n_S",&n_S,"n_S/I");
  eventtree->Branch("E_DR",&E_DR,"E_DR/F");
  eventtree->Branch("E_DRcorr",&E_DRcorr,"E_DRcorr/F");
  eventtree->Branch("Eleak_nu",&Eleak_nu,"Eleak_nu/F");
  eventtree->Branch("Pleak",&Pleak,"Pleak/F");
  eventtree->Branch("cen_Ecorr_theta",&cen_Ecorr_theta,"cen_Ecorr_theta/F");
  eventtree->Branch("cen_Ecorr_phi",&cen_Ecorr_phi,"cen_Ecorr_phi/F");
  eventtree->Branch("prt_E","vector<Float_t>",&prt_E);
  eventtree->Branch("prt_theta","vector<Float_t>",&prt_theta);
  eventtree->Branch("prt_phi","vector<Float_t>",&prt_phi);
  eventtree->Branch("prt_px","vector<Float_t>",&prt_px);
  eventtree->Branch("prt_py","vector<Float_t>",&prt_py);
  eventtree->Branch("prt_pz","vector<Float_t>",&prt_pz);
  eventtree->Branch("prt_pid","vector<Float_t>",&prt_pid);
  eventtree->Branch("ptc_E","vector<Float_t>",&ptc_E);
  eventtree->Branch("ptc_pt","vector<Float_t>",&ptc_pt);
  eventtree->Branch("ptc_p","vector<Float_t>",&ptc_p);
  eventtree->Branch("ptc_px","vector<Float_t>",&ptc_px);
  eventtree->Branch("ptc_py","vector<Float_t>",&ptc_py);
  eventtree->Branch("ptc_pz","vector<Float_t>",&ptc_pz);
  eventtree->Branch("ptc_vx","vector<Float_t>",&ptc_vx);
  eventtree->Branch("ptc_vy","vector<Float_t>",&ptc_vy);
  eventtree->Branch("ptc_vz","vector<Float_t>",&ptc_vz);
  eventtree->Branch("ptc_vt","vector<Float_t>",&ptc_vt);
  eventtree->Branch("ptc_theta","vector<Float_t>",&ptc_theta);
  eventtree->Branch("ptc_phi","vector<Float_t>",&ptc_phi);
  eventtree->Branch("ptc_eta","vector<Float_t>",&ptc_eta);
  eventtree->Branch("ptc_pid","vector<Int_t>",&ptc_pid);

  //Float_t voxel_ecor_s_[729000];//90*90*90
  //Int_t voxel_n_s_[729000];
  Float_t image_gen_[28224];//shower 168*168, jet 162*162
  Float_t point_2048_[8192];//90*90
  #define BRANCH_A_(name, size, suffix) eventtree->Branch(#name, & name##_, #name"["#size"]/"#suffix);
  #define BRANCH_AF(name, size)  BRANCH_A_(name, size, F);
  #define BRANCH_AI(name, size)  BRANCH_A_(name, size, I);
  #define FILL_ZERO(array, size) std::fill(array, array + size, 0.0);
  #define FILL_ZEROI(array, size) std::fill(array, array + size, 0);
  //BRANCH_AF(voxel_ecor_s, 729000);
  //BRANCH_AI(voxel_n_s, 729000);
  BRANCH_AF(image_gen, 28224);
  BRANCH_AI(point_2048, 8192);//don't forget to add FILL_ZERO in loop
  //FILL_ZERO(voxel_ecor_s_,729000);
  //FILL_ZEROI(voxel_n_s_,729000);
  FILL_ZERO(image_gen_,28224);
  FILL_ZEROI(point_2048_,8192);
  int depthbin=89;//22+1
  //float depthmin=34;
  //float depthmax=3060;
  float thetamin=0.79;
  float thetamax=2.35;
  float thetawidth=thetamax-thetamin;
  float phimin=-0.88;
  float phimax=0.88;
  float phiwidth=phimax-phimin;
  float etabin=56;
  float phibin=56;

  float eta=0.;
  float theta=0.;
  float phi=0.;
  Float_t w2=0.;
  cen_Ecorr_theta=0.;
  cen_Ecorr_phi=0.;
  Float_t pt_square=0.;
  Float_t m00=0.;
  Float_t m01=0.;
  Float_t m11=0.;
  float depth=0.;
  int depthindex=-1;
  int phiindex=-1;
  int thetaindex=-1;
  int voxindex=-1;
  int imgindex=-1;
  int buf_index=0;
  int diff_phi=0;
  int diff_theta=0;
  bool isopposite = false;
  TLorentzVector* ptc_pxyz = new TLorentzVector();
  Int_t count=0;
  /*
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

      auto segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDRcalo*>(m_geoSvc->lcdd()->readout(m_readoutName).segmentation().segmentation());*/
  for(int file_num=num_file_begin; file_num<num_file_end;file_num++){
    inname.Form(innameform.Data(),file_num);
    if(!fs::exists(inname.Data())){
      printf("nofile %s;\n",inname.Data());
      continue;}
    nokey=0;
    check_hep=0;
    box=TFile::Open(inname.Data(),"read");
    keys = box->GetListOfKeys();
    if(keys->GetEntries()>1){
      for(int i=0;i<keys->GetEntries();i++){
        if(strcmp(keys->At(i)->GetName(),"hepmc3_tree")==0)check_hep=1;
      }
    }
    else{
      nokey=1;
    }
    keys->Clear();
    box->Close();
    printf("loaded %s...;\n",inname.Data());
    //mychain.Add(inname+"/DRsim");
    HepMC3::ReaderRootTree hepreader(inname.Data());
    readcount+=1;
    int start=0;
    int readidx=0;
    char readpass[32];
    char key[] = "pass";
    int keyskip=0;
    unsigned int hep_entries = hepreader.m_tree->GetEntries();
    for(int hep_i = 0; hep_i<hep_entries; hep_i++){
    num_entry+=1;
    HepMC3::GenEvent genEvt;
    hepreader.read_event(genEvt);

    prt_pid.clear();
    prt_E.clear();
    prt_phi.clear();
    prt_theta.clear();
    prt_px.clear();
    prt_py.clear();
    prt_pz.clear();

      ptc_E.clear();
      ptc_pt.clear();
      ptc_p.clear();
      ptc_px.clear();
      ptc_py.clear();
      ptc_pz.clear();
      ptc_vx.clear();
      ptc_vy.clear();
      ptc_vz.clear();
      ptc_theta.clear();
      ptc_phi.clear();
      ptc_eta.clear();
      ptc_pid.clear();
    //float Etot = 0.;
    for (auto ptc : genEvt.particles()) {
      int abspid = std::abs(ptc->pid());
      int stat = ptc->status();
      //printf("%d %d\n",stat,abspid);
      if ( stat==2){
        prt_pid.push_back(ptc->pid());
        prt_E.push_back(ptc->momentum().e());
        prt_phi.push_back(ptc->momentum().phi());
        prt_theta.push_back(ptc->momentum().theta());
        prt_px.push_back(ptc->momentum().px());
        prt_py.push_back(ptc->momentum().py());
        prt_pz.push_back(ptc->momentum().pz());
        //printf("ent %d pid %d e %g stat %d\n",entry,ptc->pid(),ptc->momentum().e(),ptc->status());
      }
      if ( stat==1){
        ptc_pid.push_back(ptc->pid());
        ptc_E.push_back(ptc->momentum().e());
        ptc_phi.push_back(ptc->momentum().phi());
        ptc_theta.push_back(ptc->momentum().theta());
        ptc_px.push_back(ptc->momentum().px());
        ptc_py.push_back(ptc->momentum().py());
        ptc_pz.push_back(ptc->momentum().pz());
        //printf("ent %d pid %d e %g stat %d\n",entry,ptc->pid(),ptc->momentum().e(),ptc->status());
      }
      //if ( ptc->status() != 1 ) continue;
      //if ( abspid == 12 || abspid == 14 || abspid == 16 ) continue;
    }
    count+=1;
    Eleak_nu=0;
    Pleak =0;


    for ( num_jet=0; num_jet<1; num_jet++){
      //FILL_ZERO(voxel_ecor_s_,729000);
      //FILL_ZEROI(voxel_n_s_,729000);
      FILL_ZERO(image_gen_,28224);
      m00 = 0.;
      m01 = 0.;
      m11 = 0.;
      pt_square = 0.;
      pt_Gen = 0.;
      p_Gen = 0.;
      cmult=0;
      nmult=0;
      chad_mult=0;
      nhad_mult=0;
      electron_mult=0;
      muon_mult=0;
      photon_mult=0;
      mult=0;
      E_Gen=0.;
      deltheta2=0.;
      delphi2=0.;
      width_Gen=0.;
      center_theta_gen=0.;
      center_phi_gen=0.;
      cen_Ecorr_theta=0.;
      cen_Ecorr_phi=0.;
      not_theta=0;
      not_phi=0;
      TLorentzVector ptc_sum;
      imgindex=-1;
      phiindex=-1;
      thetaindex=-1;
      /*
      for (auto genptc : drEvt.GenPtcs){
        ptc_pxyz->SetPxPyPzE(genptc.px,genptc.py,genptc.pz,genptc.E);
        phi=float(ptc_pxyz->Phi());
        isopposite=false;
        if(ptc_pxyz->Theta()>thetamin && ptc_pxyz->Theta()<thetamax){
          if(phi>phimin && phi<phimax){
            thetaindex=int(162*(ptc_pxyz->Theta()-thetamin)/thetawidth);
            phiindex=int(162*(phi-phimin)/phiwidth);
            imgindex=162*thetaindex+phiindex;
            image_gen_[imgindex]+=genptc.E;
          }
        }
        center_theta_gen+=ptc_pxyz->Theta()*ptc_pxyz->Pt();
        center_phi_gen+=phi*ptc_pxyz->Pt();
        pt_Gen+=ptc_pxyz->Pt();
        p_Gen+=ptc_pxyz->P();
        E_Gen+=genptc.E;
        ptc_sum=ptc_sum+*ptc_pxyz;
      }
      //printf("mass %g\n",ptc_sum.M());
      mass_Gen=ptc_sum.M();
      center_theta_gen=center_theta_gen/pt_Gen;
      center_phi_gen=center_phi_gen/pt_Gen;
      for (auto genptc : drEvt.GenPtcs){
        ptc_pxyz->SetPxPyPzE(genptc.px,genptc.py,genptc.pz,genptc.E);
        phi=float(ptc_pxyz->Phi());
        isopposite=false;
        deltheta2+=std::pow(ptc_pxyz->Theta()-center_theta_gen,2)*ptc_pxyz->Pt();
        delphi2+=std::pow(phi-center_phi_gen,2)*ptc_pxyz->Pt();
        ptc_theta.push_back(ptc_pxyz->Theta());
        ptc_eta.push_back(ptc_pxyz->Eta());
        ptc_phi.push_back(phi);
        ptc_pid.push_back(genptc.pdgId);
        ptc_px.push_back(genptc.px/1000.);
        ptc_py.push_back(genptc.py/1000.);
        ptc_pz.push_back(genptc.pz/1000.);
        ptc_vx.push_back(genptc.vx);
        ptc_vy.push_back(genptc.vy);
        ptc_vz.push_back(genptc.vz);
        ptc_vt.push_back(genptc.vt);
        ptc_E.push_back(genptc.E/1000.);
        ptc_pt.push_back(ptc_pxyz->Pt()/1000.);
        ptc_p.push_back(ptc_pxyz->P()/1000.);
        mult++;
        if(genptc.E>100.){
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

        x = ptc_pxyz->Theta();
        y = phi;
        w2 = std::pow(ptc_pxyz->Pt(), 2);
        m00 += w2 * std::pow(x, 2);
        m01 -= w2 * x * y;
        m11 += w2 * std::pow(y, 2);
        pt_square+=w2;

      }
      pt_Gen=pt_Gen/1000.;
      p_Gen=p_Gen/1000.;
      E_Gen=E_Gen/1000.;
      Edep=0;
      if(1){
      deltheta2=deltheta2/pt_Gen;
      delphi2=delphi2/pt_Gen;
      width_Gen=delphi2+deltheta2;
      TMatrixFSym covariance_matrix(2);
      covariance_matrix(0, 0) = m00;
      covariance_matrix(0, 1) = m01;
      covariance_matrix(1, 1) = m11;
      TVectorF eigen_values;
      covariance_matrix.EigenVectors(eigen_values);
      major_axis = std::sqrt(eigen_values[0] / pt_square);
      minor_axis = std::sqrt(eigen_values[1] / pt_square);
      ptd=std::sqrt(pt_square)/pt_Gen;
      }
      */
     }
      //printf("entry end..\n");
      eventtree->Fill();
   }
    hepreader.close();

  }

//TFile outfile(outname.Data(), "recreate");
printf("%d files %d events filled in %s\n",readcount,int(eventtree->GetEntries()),outname.Data());
  //printf("writing..\n");
outfile->Write("");
//outfile->Write("",TObject::kOverwrite);
  //printf("closing..\n");
outfile->Close();
printf("--done\n");

return 0;
}
