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

class process_class {
    public:
      process_class(){
        /*
  ifstream towerpos_file("/pad/yulee/dream/repo/0023_img/tools/towerposition.txt");

  if(!towerpos_file) {
    printf("towerposition.txt not found\n");
    return 1;
  }
float temp;
count=0;
while (towerpos_file >> temp) {
    dest.push_back(temp);
    count+=1;
  }
towerpos_file.close();

towerpos_file=open("/pad/yulee/dream/repo/0023_img/tools/towerposition.txt","r")
towerpos=towerpos_file.readlines()
phi_list=[]
theta_list=[]
for i in range(283):
    phi_list.append(float(towerpos[1+i*92].split(" ")[2]))
for i in range(92):
    theta_list.append(float(towerpos[1+i].split(" ")[3]))
    */
        float phi_list[283]={0.0,0.0222021,0.0444041,0.0666062,0.0888083,0.11101,0.133212,
        0.155414,0.177617,0.199819,0.222021,0.244223,0.266425,0.288627,0.310829,
        0.333031,0.355233,0.377435,0.399637,0.421839,0.444041,0.466243,0.488445,
        0.510648,0.53285,0.555052,0.577254,0.599456,0.621658,0.64386,0.666062,
        0.688264,0.710466,0.732668,0.75487,0.777072,0.799274,0.821477,0.843679,
        0.865881,0.888083,0.910285,0.932487,0.954689,0.976891,0.999093,1.0213,
        1.0435,1.0657,1.0879,1.1101,1.13231,1.15451,1.17671,1.19891,
        1.22111,1.24332,1.26552,1.28772,1.30992,1.33212,1.35433,1.37653,
        1.39873,1.42093,1.44313,1.46534,1.48754,1.50974,1.53194,1.55414,
        1.57635,1.59855,1.62075,1.64295,1.66516,1.68736,1.70956,1.73176,
        1.75396,1.77617,1.79837,1.82057,1.84277,1.86497,1.88718,1.90938,
        1.93158,1.95378,1.97598,1.99819,2.02039,2.04259,2.06479,2.08699,
        2.1092,2.1314,2.1536,2.1758,2.198,2.22021,2.24241,2.26461,
        2.28681,2.30902,2.33122,2.35342,2.37562,2.39782,2.42003,2.44223,
        2.46443,2.48663,2.50883,2.53104,2.55324,2.57544,2.59764,2.61984,
        2.64205,2.66425,2.68645,2.70865,2.73085,2.75306,2.77526,2.79746,
        2.81966,2.84186,2.86407,2.88627,2.90847,2.93067,2.95288,2.97508,
        2.99728,3.01948,3.04168,3.06389,3.08609,3.10829,3.13049,-3.13049,
        -3.10829,-3.08609,-3.06389,-3.04168,-3.01948,-2.99728,-2.97508,-2.95288,
        -2.93067,-2.90847,-2.88627,-2.86407,-2.84186,-2.81966,-2.79746,-2.77526,
        -2.75306,-2.73085,-2.70865,-2.68645,-2.66425,-2.64205,-2.61984,-2.59764,
        -2.57544,-2.55324,-2.53104,-2.50883,-2.48663,-2.46443,-2.44223,-2.42003,
        -2.39782,-2.37562,-2.35342,-2.33122,-2.30902,-2.28681,-2.26461,-2.24241,
        -2.22021,-2.198,-2.1758,-2.1536,-2.1314,-2.1092,-2.08699,-2.06479,
        -2.04259,-2.02039,-1.99819,-1.97598,-1.95378,-1.93158,-1.90938,-1.88718,
        -1.86497,-1.84277,-1.82057,-1.79837,-1.77617,-1.75396,-1.73176,-1.70956,
        -1.68736,-1.66516,-1.64295,-1.62075,-1.59855,-1.57635,-1.55414,-1.53194,
        -1.50974,-1.48754,-1.46534,-1.44313,-1.42093,-1.39873,-1.37653,-1.35433,
        -1.33212,-1.30992,-1.28772,-1.26552,-1.24332,-1.22111,-1.19891,-1.17671,
        -1.15451,-1.13231,-1.1101,-1.0879,-1.0657,-1.0435,-1.0213,-0.999093,
        -0.976891,-0.954689,-0.932487,-0.910285,-0.888083,-0.865881,-0.843679,-0.821477,
        -0.799274,-0.777072,-0.75487,-0.732668,-0.710466,-0.688264,-0.666062,-0.64386,
        -0.621658,-0.599456,-0.577254,-0.555052,-0.53285,-0.510648,-0.488445,-0.466243,
        -0.444041,-0.421839,-0.399637,-0.377435,-0.355233,-0.333031,-0.310829,-0.288627,
        -0.266425,-0.244223,-0.222021,-0.199819,-0.177617,-0.155414,-0.133212,-0.11101,
        -0.0888083,-0.0666062,-0.0444041,-0.0222021};
        float theta_list[92]={1.55969,1.53748,1.51529,1.49314,1.47102,1.44896,1.42697,
        1.40505,1.38321,1.36147,1.33984,1.31832,1.29692,1.27566,1.25454,
        1.23357,1.21276,1.19212,1.17164,1.15135,1.13124,1.11132,1.0916,
        1.07208,1.05278,1.03369,1.01482,0.996166,0.977741,0.959546,0.941581,
        0.923851,0.906356,0.889096,0.872081,0.855311,0.838781,0.822496,0.806456,
        0.790661,0.775111,0.759806,0.744746,0.729926,0.715351,0.701021,0.686931,
        0.673081,0.659466,0.646086,0.632941,0.620026,0.607226,0.594426,0.581626,
        0.568826,0.556026,0.543226,0.530426,0.517626,0.504826,0.492026,0.479226,
        0.466426,0.453626,0.440826,0.428026,0.415226,0.402426,0.389626,0.376826,
        0.364026,0.351226,0.338426,0.325626,0.312826,0.300026,0.287226,0.274426,
        0.261626,0.248826,0.236026,0.223226,0.210426,0.197626,0.184826,0.172026,
        0.159226,0.146426,0.133626,0.120826,0.108026};
        for(int i =0;i<283;i++){
          _phi_list[i]=phi_list[i];
          if(i<92)_theta_list[i]=theta_list[i];
        }
        _pi=TMath::Pi();
      }

      int findphi(float phi){
        float diff=0;
        float mindiff=10;
        int idx=0;

        if(phi>2*_pi){
            phi=phi-2*_pi;
        }
        for(int i=0;i<283;i++){
            diff=abs(phi-_phi_list[i]);
            if(diff<mindiff){
                mindiff=diff;
                idx=i;
            }
        }
        return idx;
      }

      int findtheta(float theta){
        float diff=0;
        float mindiff=10;
        int idx=0;

        for(int i=0;i<92;i++){
            diff=abs(theta-_theta_list[i]);
            if(diff<mindiff){
                mindiff=diff;
                idx=i;
            }
            if(i!=91){
                diff=abs(theta-(_pi-_theta_list[i]));
                if(diff<mindiff){
                  mindiff=diff;
                  idx=-i-1;
                }
            }
        }
        return idx;
      }

    private:
      float _phi_list[283];
      float _theta_list[92];
      float _pi;
};

//./process ../../el2entest2_%d.root 0 1 ../../elentest2 0 0
//./process /hdfs/user/yulee/DRsim/simdata/cleanhad/${pid}_${energy}GeV_100_${i}_%d.root 0 1 $here/pack/${pid}_${energy}GeV_${i}.root 1 $here/hadbox/${pid}_${energy}GeV_100_${i}_%d_Reco.root
int main(int argc, char* argv[]){
  TString innameform = argv[1];
  int num_file_begin = std::stoi(argv[2]);
  int num_file_end = std::stoi(argv[3]);
  TString filename = argv[4];
  int isjet = atoi(argv[5]);
  TString reconameform;
  if(argc>6)reconameform = argv[6];
  TString inname;
  TString reconame;
  TFile* box;
  TList* keys;
  int readcount=0;
  printf("process file will be checked %d\n",num_file_end-num_file_begin);
  int check_hep=0;
  int check_dr=0;
  int check_reco=0;
  int nokey=0;
  for(int file_num=num_file_begin; file_num<num_file_end;file_num++){
    inname.Form(innameform.Data(),file_num);
    reconame=inname.Copy();
    reconame.ReplaceAll(".root","_Reco.root");
    if(argc>6)reconame.Form(reconameform.Data(),file_num);
    if(!fs::exists(inname.Data())){
      printf("nofile %s;\n",inname.Data());
      continue;}
    //printf("file loaded %s..%d.;\n",inname.Data(),keys->GetEntries());
    nokey=0;
    check_hep=0;
    check_dr=0;
    check_reco=0;
    box=TFile::Open(inname.Data(),"read");
    keys = box->GetListOfKeys();
    keys = box->GetListOfKeys();
    if(keys->GetEntries()>1){
      for(int i=0;i<keys->GetEntries();i++){
        if(strcmp(keys->At(i)->GetName(),"hepmc3_tree")==0)check_hep=1;
        if(strcmp(keys->At(i)->GetName(),"DRsim")==0)check_dr=1;
        if(strcmp(keys->At(i)->GetName(),"Reco")==0)check_reco=1;
      }
    }
    else{
      nokey=1;
    }
    if(check_hep==0 && isjet==1) nokey=1;
    if(check_dr==0) nokey=1;
    keys->Clear();
    box->Close();
    if(nokey==1){
      printf("nokey hep %d dr %d reco %d %s %s;\n",check_hep,check_dr,check_reco,inname.Data(),reconame.Data());
      continue;
    }
    if(fs::exists(reconame.Data())){
        box=TFile::Open(reconame.Data(),"read");
        keys = box->GetListOfKeys();
        //printf("%d %s\n",keys->GetEntries(),reconame.Data());
        if(keys->GetEntries()>0){
          //printf("%d %s\n",keys->GetEntries(),keys->At(0)->GetName());
          if(strcmp(keys->At(0)->GetName(),"Reco")==0)check_reco=1;
        }
        keys->Clear();
        box->Close();
    }
    else{
      printf("nofile %s;\n",reconame.Data());
      continue;
    }
    if(check_reco==0) nokey=1;
    if(nokey==1){
      printf("nokey hep %d dr %d reco %d %s %s;\n",check_hep,check_dr,check_reco,inname.Data(),reconame.Data());
      continue;
    }
    readcount+=1;
  }
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
  std::vector<Int_t> tower_idx_eta;
  std::vector<Int_t> tower_idx_phi;
  std::vector<Int_t> tower_no_eta;
  std::vector<Int_t> tower_no_phi;
  std::vector<Float_t> tower_eta;
  std::vector<Float_t> tower_theta;
  std::vector<Float_t> tower_phi;
  std::vector<Int_t> tower_notheta;
  std::vector<Int_t> tower_nophi;
  std::vector<Int_t> tower_diff_theta;
  std::vector<Int_t> tower_diff_phi;
  std::vector<Int_t> tower_numx;
  std::vector<Int_t> tower_numy;
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
  float E_C_img=0.;
  float E_Scorr_img=0.;
  int n_C=0.;
  int n_S=0.;
  int idx_phi=0.;
  int idx_theta=0.;
  float E_DR=0.;
  float E_DRcorr=0.;
  float Eleak_nu=0;
  float Pleak =0;
  float ecimg[10000]={0};
  float esimg[10000]={0};
  int same_ec=0;
  int same_es=0;
  Float_t center_theta_img=0.;
  Float_t center_phi_img=0.;
  Float_t center_theta_gen=0.;
  Float_t center_phi_gen=0.;
  Float_t cen_Ecorr_theta=0.;
  Float_t cen_Ecorr_phi=0.;
  Int_t num_jet=0;
  Int_t num_entry=-1;
  Int_t not_phi=0;
  Int_t not_theta=0;
  std::vector<Float_t> tower_e_s;
  std::vector<Float_t> tower_e_c;
  std::vector<Int_t> tower_n_s;
  std::vector<Int_t> tower_n_c;
  std::vector<Float_t> tower_ecor_s;
  std::vector<Float_t> tower_ecor_dr;
  std::vector<Float_t> fiber_energy;
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
  eventtree->Branch("tower_eta","vector<Float_t>",&tower_eta);
  eventtree->Branch("tower_theta","vector<Float_t>",&tower_theta);
  eventtree->Branch("tower_phi","vector<Float_t>",&tower_phi);
  eventtree->Branch("tower_notheta","vector<Int_t>",&tower_notheta);
  eventtree->Branch("tower_nophi","vector<Int_t>",&tower_nophi);
  eventtree->Branch("tower_diff_theta","vector<Int_t>",&tower_diff_theta);
  eventtree->Branch("tower_diff_phi","vector<Int_t>",&tower_diff_phi);
  eventtree->Branch("tower_idx_eta","vector<Int_t>",&tower_idx_eta);
  eventtree->Branch("tower_idx_phi","vector<Int_t>",&tower_idx_phi);
  eventtree->Branch("tower_no_eta","vector<Int_t>",&tower_no_eta);
  eventtree->Branch("tower_no_phi","vector<Int_t>",&tower_no_phi);
  eventtree->Branch("tower_numx","vector<Int_t>",&tower_numx);
  eventtree->Branch("tower_numy","vector<Int_t>",&tower_numy);
  eventtree->Branch("ptd",&ptd,"ptd/F");
  eventtree->Branch("major_axis",&major_axis,"major_axis/F");
  eventtree->Branch("minor_axis",&minor_axis,"minor_axis/F");
  eventtree->Branch("center_theta_gen",&center_theta_gen,"center_theta_gen/F");
  eventtree->Branch("center_phi_gen",&center_phi_gen,"center_phi_gen/F");
  eventtree->Branch("center_theta_img",&center_theta_img,"center_theta_img/F");
  eventtree->Branch("center_phi_img",&center_phi_img,"center_phi_img/F");
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
  eventtree->Branch("E_C_img",&E_C_img,"E_C_img/F");
  eventtree->Branch("E_S",&E_S,"E_S/F");
  eventtree->Branch("E_Scorr",&E_Scorr,"E_Scorr/F");
  eventtree->Branch("E_Scorr_img",&E_Scorr_img,"E_Scorr_img/F");
  eventtree->Branch("n_C",&n_C,"n_C/I");
  eventtree->Branch("n_S",&n_S,"n_S/I");
  eventtree->Branch("E_DR",&E_DR,"E_DR/F");
  eventtree->Branch("E_DRcorr",&E_DRcorr,"E_DRcorr/F");
  eventtree->Branch("Eleak_nu",&Eleak_nu,"Eleak_nu/F");
  eventtree->Branch("Pleak",&Pleak,"Pleak/F");
  eventtree->Branch("tower_e_s","vector<Float_t>",&tower_e_s);
  eventtree->Branch("tower_e_c","vector<Float_t>",&tower_e_c);
  eventtree->Branch("tower_n_s","vector<Int_t>",&tower_n_s);
  eventtree->Branch("tower_n_c","vector<Int_t>",&tower_n_c);
  eventtree->Branch("tower_ecor_s","vector<Float_t>",&tower_ecor_s);
  eventtree->Branch("tower_ecor_dr","vector<Float_t>",&tower_ecor_dr);
  eventtree->Branch("fiber_energy","vector<Float_t>",&fiber_energy);
  eventtree->Branch("fiber_ecor","vector<Float_t>",&fiber_ecor);
  eventtree->Branch("fiber_ecor_s","vector<Float_t>",&fiber_ecor_s);
  eventtree->Branch("fiber_ecor_c","vector<Float_t>",&fiber_ecor_c);
  eventtree->Branch("fiber_n","vector<Int_t>",&fiber_n);
  eventtree->Branch("fiber_itower","vector<Int_t>",&fiber_itower);
  eventtree->Branch("fiber_ix","vector<Int_t>",&fiber_ix);
  eventtree->Branch("fiber_iy","vector<Int_t>",&fiber_iy);
  eventtree->Branch("fiber_t","vector<Float_t>",&fiber_t);
  eventtree->Branch("fiber_x","vector<Float_t>",&fiber_x);
  eventtree->Branch("fiber_y","vector<Float_t>",&fiber_y);
  eventtree->Branch("fiber_z","vector<Float_t>",&fiber_z);
  eventtree->Branch("fiber_depth","vector<Float_t>",&fiber_depth);
  eventtree->Branch("fiber_r","vector<Float_t>",&fiber_r);
  eventtree->Branch("fiber_theta","vector<Float_t>",&fiber_theta);
  eventtree->Branch("fiber_phi","vector<Float_t>",&fiber_phi);
  eventtree->Branch("fiber_eta","vector<Float_t>",&fiber_eta);
  eventtree->Branch("fiber_iscerenkov","vector<Bool_t>",&fiber_iscerenkov);
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
  Float_t point_ecor_s_[2000];//point 2000
  Float_t point_ecor_c_[2000];//point 2000
  Int_t point_phi_idx_[2000];//point 2000
  Int_t point_theta_idx_[2000];//point 2000
  Float_t point_depth_[2000];//point 2000
  Int_t point_n_[10000];//point 2000*5 phi eta s c depth
  Float_t image_gen_[28224];//shower 168*168, jet 162*162
  Float_t image_ecor_c_[28224];//168*168
  Int_t image_n_c_[28224];
  Float_t image_ecor_s_[28224];//168*168
  Int_t image_n_s_[28224];
  Float_t point_2048_[8192];//90*90
  #define BRANCH_A_(name, size, suffix) eventtree->Branch(#name, & name##_, #name"["#size"]/"#suffix);
  #define BRANCH_AF(name, size)  BRANCH_A_(name, size, F);
  #define BRANCH_AI(name, size)  BRANCH_A_(name, size, I);
  #define FILL_ZERO(array, size) std::fill(array, array + size, 0.0);
  #define FILL_ZEROI(array, size) std::fill(array, array + size, 0);
  //BRANCH_AF(voxel_ecor_s, 729000);
  //BRANCH_AI(voxel_n_s, 729000);
  BRANCH_AF(point_ecor_s,2000);
  BRANCH_AF(point_ecor_c,2000);
  BRANCH_AI(point_phi_idx,2000);
  BRANCH_AI(point_theta_idx,2000);
  BRANCH_AF(point_depth,2000);
  BRANCH_AI(point_n,8000);
  BRANCH_AF(image_gen, 28224);
  BRANCH_AF(image_ecor_c, 28224);
  BRANCH_AI(image_n_c, 28224);//don't forget to add FILL_ZERO in loop
  BRANCH_AF(image_ecor_s, 28224);
  BRANCH_AI(image_n_s, 28224);//don't forget to add FILL_ZERO in loop
  BRANCH_AI(point_2048, 8192);//don't forget to add FILL_ZERO in loop
  //FILL_ZERO(voxel_ecor_s_,729000);
  //FILL_ZEROI(voxel_n_s_,729000);
  FILL_ZERO(image_gen_,28224);
  FILL_ZERO(image_ecor_c_,28224);
  FILL_ZEROI(image_n_c_,28224);
  FILL_ZERO(image_ecor_s_,28224);
  FILL_ZEROI(image_n_s_,28224);
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
  float f_threshold=0.001;
  float f_count=0;

int width=40;
int ymax=56;
int xmax=56;
int pix_per_tower=56;
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
  int fibercount=0;
  bool isopposite = false;
  int skipimg=0;
    int fibercheck=0;
  Float_t buf_fiber_e[5];
  Float_t buf_fiber_phi[5];
  Int_t buf_fiber_same=0;
  Int_t outrange=0;
  Int_t repeatfiber=0;
  Int_t nofiber=0;

  TVector3* fiberxyz = new TVector3();
  TVector3* towerxyz = new TVector3();
  TLorentzVector* ptc_pxyz = new TLorentzVector();
  FILE *fp1;
  Int_t count=0;
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
  
  for(int file_num=num_file_begin; file_num<num_file_end;file_num++){
    inname.Form(innameform.Data(),file_num);
    reconame=inname.Copy();
    reconame.ReplaceAll(".root","_Reco.root");
    if(argc>6){
      reconame.Form(reconameform.Data(),file_num);
      //printf("load reco file\n");
    }
    if(!fs::exists(inname.Data())){
      printf("nofile %s;\n",inname.Data());
      continue;}
    nokey=0;
    check_hep=0;
    check_dr=0;
    check_reco=0;
    box=TFile::Open(inname.Data(),"read");
    keys = box->GetListOfKeys();
    if(keys->GetEntries()>1){
      for(int i=0;i<keys->GetEntries();i++){
        if(strcmp(keys->At(i)->GetName(),"hepmc3_tree")==0)check_hep=1;
        if(strcmp(keys->At(i)->GetName(),"DRsim")==0)check_dr=1;
        if(strcmp(keys->At(i)->GetName(),"Reco")==0)check_reco=1;
      }
    }
    else{
      nokey=1;
    }
    if(check_hep==0 && isjet==1) nokey=1;
    if(check_dr==0) nokey=1;
    keys->Clear();
    box->Close();
    if(fs::exists(reconame.Data())){
        box=TFile::Open(reconame.Data(),"read");
        keys = box->GetListOfKeys();
        if(keys->GetEntries()>0){
          if(strcmp(keys->At(0)->GetName(),"Reco")==0)check_reco=1;
        }
        keys->Clear();
        box->Close();
    }
    else{
      printf("nofile %s;\n",reconame.Data());
      continue;
    }
    if(check_reco==0) nokey=1;
    if(nokey==1){
      printf("nokey hep %d dr %d reco %d %s;\n",check_hep,check_dr,check_reco,inname.Data());
      continue;
    }
    printf("loaded %s...;\n",inname.Data());
    //mychain.Add(inname+"/DRsim");
    HepMC3::ReaderRootTree hepreader(inname.Data());
    RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(inname.Data()));
    //if(!fs::exists(reconame.Data()))reconame=inname.Copy();
    //if(argc>6)reconame.Form(reconameform.Data(),file_num);
    RootInterface<RecoInterface::RecoEventData>* recoInterface = new RootInterface<RecoInterface::RecoEventData>(std::string(reconame.Data()));
    //RootInterface<RecoInterface::RecoEventData>* recoInterface = new RootInterface<RecoInterface::RecoEventData>("../../tools/box/gg_70GeV_799_Reco.root");
    //printf("Loading DR interface\n");
    drInterface->set("DRsim","DRsimEventData");
    //printf("Loading Reco interface\n");
    recoInterface->set("Reco","RecoEventData");
    readcount+=1;
    unsigned int entries = recoInterface->entries();
     
    int start=0;
    int readidx=0;
    char readpass[32];
    char key[] = "pass";
    int keyskip=0;
    while (recoInterface->numEvt() < entries) {
    //if(count%int(entries/10.)==0)printf("%d%%--\n",int(100.*count/entries));
    num_entry+=1;
    fibercheck=0;
    HepMC3::GenEvent genEvt;
    hepreader.read_event(genEvt);

    prt_pid.clear();
    prt_E.clear();
    prt_phi.clear();
    prt_theta.clear();
    prt_px.clear();
    prt_py.clear();
    prt_pz.clear();

    //float Etot = 0.;
    for (auto mptc : genEvt.particles()) {
      int abspid = std::abs(mptc->pid());
      int stat = mptc->status();
      if ( stat==11){
        for (auto ptc : mptc->children()) {
          abspid = std::abs(ptc->pid());
          stat = ptc->status();
          prt_pid.push_back(ptc->pid());
          prt_E.push_back(ptc->momentum().e());
          prt_phi.push_back(ptc->momentum().phi());
          prt_theta.push_back(ptc->momentum().theta());
          prt_px.push_back(ptc->momentum().px());
          prt_py.push_back(ptc->momentum().py());
          prt_pz.push_back(ptc->momentum().pz());
          //printf("ent %d pid %d e %g stat %d\n",entry,ptc->pid(),ptc->momentum().e(),ptc->status());
        }
      }
      //if ( ptc->status() != 1 ) continue;
      //if ( abspid == 12 || abspid == 14 || abspid == 16 ) continue;
    }
    DRsimInterface::DRsimEventData drEvt;
    RecoInterface::RecoEventData recoEvt;
    drInterface->read(drEvt);
    recoInterface->read(recoEvt);
    count+=1;
    E_C=recoEvt.E_C;
    E_S=recoEvt.E_S;
    E_Scorr=recoEvt.E_Scorr;
    //n_C=recoEvt.n_C;
    //n_S=recoEvt.n_S;
    if(E_DR==recoEvt.E_DR)printf("same energy %d\n",count);
    E_DR=recoEvt.E_DR;
    E_DRcorr=recoEvt.E_DRcorr;
    Eleak_nu=0;
    Pleak =0;
    tower_eta.clear();
    tower_theta.clear();
    tower_phi.clear();
    tower_notheta.clear();
    tower_nophi.clear();
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

    fibercheck=0;
    for ( num_jet=0; num_jet<2; num_jet++){
      if(num_jet==1 && isjet==0) break;
      fibercount=0;
      fiber_energy.clear();
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
      tower_idx_eta.clear();
      tower_idx_phi.clear();
      tower_no_eta.clear();
      tower_no_phi.clear();
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
      tower_diff_theta.clear();
      tower_diff_phi.clear();
      //FILL_ZERO(voxel_ecor_s_,729000);
      //FILL_ZEROI(voxel_n_s_,729000);
      FILL_ZERO(point_ecor_s_,2000);
      FILL_ZERO(point_ecor_c_,2000);
      FILL_ZERO(point_phi_idx_,2000);
      FILL_ZERO(point_theta_idx_,2000);
      FILL_ZERO(point_depth_,2000);
      FILL_ZERO(point_n_,10000);
      FILL_ZERO(image_gen_,28224);
      FILL_ZERO(image_ecor_c_,28224);
      FILL_ZEROI(image_n_c_,28224);
      FILL_ZERO(image_ecor_s_,28224);
      FILL_ZEROI(image_n_s_,28224);
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
      for (auto genptc : drEvt.GenPtcs){
        ptc_pxyz->SetPxPyPzE(genptc.px,genptc.py,genptc.pz,genptc.E);
        phi=float(ptc_pxyz->Phi());
        isopposite=false;
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
        if(ptc_pxyz->Theta()>thetamin && ptc_pxyz->Theta()<thetamax){
          if(phi>phimin && phi<phimax){
            thetaindex=int(162*(ptc_pxyz->Theta()-thetamin)/thetawidth);
            phiindex=int(162*(phi-phimin)/phiwidth);
            if(isjet==1) imgindex=162*thetaindex+phiindex;
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
      process_class process;
      int cenphi=process.findphi(center_phi_gen+num_jet*pi);
      int centheta=process.findtheta(center_theta_gen);
      for (auto genptc : drEvt.GenPtcs){
        ptc_pxyz->SetPxPyPzE(genptc.px,genptc.py,genptc.pz,genptc.E);
        phi=float(ptc_pxyz->Phi());
        isopposite=false;
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
      num_tower = -1;
      TVector3* towercenterxyz = new TVector3();//should be energy center

      auto towercenterpos = segmentation->towerposition(0,0);
      towercenterxyz->SetXYZ(towercenterpos.x(),towercenterpos.y(),towercenterpos.z());
      auto towercenterphi=float(towercenterxyz->Phi());
      auto towercentertheta=float(towercenterxyz->Theta());
      center_theta_img=0.;
      center_phi_img=0.;
      E_Scorr_img=0.;
      E_C_img=0.;
      buf_fiber_same=0;

      for (auto tower : recoEvt.towers) {
        //towerData.iTheta = segmentation->numEta(hit->GetSiPMnum());
        //towerData.iPhi = segmentation->numPhi(hit->GetSiPMnum());
        int noTheta = int(tower.iTheta);// -1,0,1 ...
        if(abs(noTheta)>91)continue;
        int noPhi = int(tower.iPhi);//282,0,1 ...
        if(isjet==0){
          if(abs(noTheta)>92 || noPhi<0 || noPhi>282){
            outrange+=1;
            fibercheck=0;
            break;
          }
        }
        if(buf_fiber_same>=5){
          repeatfiber+=1;
          fibercheck=0;
          break;
        }
        auto towerpos = segmentation->towerposition(noTheta,noPhi);
        towerxyz->SetXYZ(towerpos.x(),towerpos.y(),towerpos.z());
        float towerphi=float(towerxyz->Phi());
        float towereta=float(towerxyz->Eta());
        float towertheta=float(towerxyz->Theta());
        isopposite=false;
        if(isjet==1 && num_jet==0 && abs(towerphi)>pi/2.)continue;
        if(isjet==1 && num_jet==1){
          if(abs(towerphi)>pi/2.){
            isopposite=true;
            if(towerphi<0){
              towerphi=pi+towerphi;
            }
            else{
              towerphi=towerphi-pi;
            }
          }
          else continue;
        }
        tower_e_s.push_back(float(tower.E_S));
        tower_e_c.push_back(float(tower.E_C));
        tower_ecor_s.push_back(tower.E_Scorr);
        tower_ecor_dr.push_back(tower.E_DRcorr);
        //tower_n_s.push_back(int(tower.n_S));
        //tower_n_c.push_back(int(tower.n_C));
        tower_numx.push_back(tower.numx);
        tower_numy.push_back(tower.numy);
        tower_phi.push_back(towerphi);
        tower_eta.push_back(towereta);
        tower_theta.push_back(towertheta);
        tower_nophi.push_back(noPhi);
        tower_notheta.push_back(noTheta);
        diff_phi=int(1 +TMath::Nint((towerphi-towercenterphi)/0.022));
        diff_theta=int(1 -TMath::Nint((towertheta-towercentertheta)/0.022));
        if(isjet==1){
          diff_phi=int(22 +TMath::Nint((towerphi-towercenterphi)/0.022));
          diff_theta=int(22 +TMath::Nint((towertheta-towercentertheta)/0.022));
        }
        tower_diff_theta.push_back(diff_theta);
        tower_diff_phi.push_back(diff_phi);
        num_tower+=1;
        int checktower=0;
        int fiberintowercount=0;
        skipimg=0;
        if(isopposite==false){
          if(noTheta+40<0 || noTheta+40>80) skipimg=1;
          if(noPhi-282+40-1<0 && noPhi+40>80)skipimg=1;//phi 0 방향일때만
        }
        else{
          if(noTheta+40+1<0 || noTheta+40+1>80) skipimg=1;
          if(noPhi-142+40<0 || noPhi-142+40>80)skipimg=1;//phi pi 방향일때만
        }
        for (auto fiber : tower.fibers) {
            if(fiber.Ecorr==0)continue;
            if(noTheta!=(segmentation->numEta(fiber.fiberNum))){
              not_theta=1;
              //printf("no theta different\n");//tower 4 개로 나눠서 픽셀화할예정
            }
            if(noPhi!=(segmentation->numPhi(fiber.fiberNum))){
              not_phi=1;
              //printf("no phi different\n");
             }
            idx_theta=segmentation->numEta(fiber.fiberNum);//==notheta
            idx_phi=segmentation->numPhi(fiber.fiberNum);//==nophi
        if(abs(idx_theta)>91){
              printf("fiber fibernum %d notheta %d nophi %d\n",fiber.fiberNum,idx_theta,idx_phi);
              continue;
            }
            auto pos = segmentation->position(fiber.fiberNum);
            //isopposite=false;
            x=float(pos.x());
            y=float(pos.y());
            z=float(pos.z());
            fiberxyz->SetXYZ(x,y,z);
            phi=float(fiberxyz->Phi());
            /*if(num_tower==0&&num_jet==0&&fibercount<5&&fiber.Ecorr>0){
              if(buf_fiber_e[fibercount]==fiber.Ecorr and buf_fiber_phi[fibercount]==phi)buf_fiber_same+=1;
              buf_fiber_e[fibercount]=fiber.Ecorr;
              buf_fiber_phi[fibercount]=phi;
            }*/
            fiberintowercount+=1;
            eta=float(fiberxyz->Eta());
            theta=float(fiberxyz->Theta());
            if(isjet==1 && num_jet==0 && abs(phi)>pi/2.)continue; // ijset 할때 phi 반대쪽에서 nophi(phi index)도 연동해야 반대편이면 nophi가 반대쪽임
            if(isjet==1 && num_jet==1){
              if(abs(phi)>pi/2.){
                //isopposite=true;
                if(phi<0){
                  phi=pi+phi;
                }
                else{
                  phi=phi-pi;
                }
              }
              else continue;
            }
            diff_theta=int(22 +TMath::Nint((towertheta-towercentertheta)/0.022));
            diff_phi=int(22 +TMath::Nint((towertheta-towercentertheta)/0.022));
            tower_idx_eta.push_back(idx_theta);
            tower_idx_phi.push_back(idx_phi);
            tower_no_eta.push_back(noTheta);//tower 4 개로 나눠서 픽셀화할예정
            tower_no_phi.push_back(noPhi);
            /*if(isopposite==false){#not_phi,not_theta 있는경우엔 다시 사용, 원래는 서로 일치해야하나 파일이 깨지면
              if(idx_theta+40<0 || idx_theta+40>80) skipimg=1;
              if(idx_phi-282+40-1<0 && idx_phi+40>80)skipimg=1;//phi 0 방향일때만
            }
            else{
              if(idx_theta+40+1<0 || idx_theta+40+1>80) skipimg=1;
              if(idx_phi-142+40<0 || idx_phi-142+40>80)skipimg=1;//phi pi 방향일때만
            }*/
            checktower=1;
            fibercheck=1;
            depthindex=-1;
            phiindex=-1;
            thetaindex=-1;
            voxindex=-1;
            imgindex=-1;
            fiber_x.push_back(x);
            fiber_y.push_back(y);
            fiber_z.push_back(z);
            fiber_iscerenkov.push_back(segmentation->IsCerenkov(fiber.fiberNum));
            fiber_energy.push_back(float(fiber.E));
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
            int fillindex=-1;
//notheta=noTheta
          if(noTheta<centheta-width-0 or noTheta>centheta+width-0){
                  fillindex=-2;
          }
          if(cenphi>width and cenphi+width<=282){ // range between 0~282
              if(noPhi<cenphi-width or noPhi>cenphi+width){
                  fillindex=-2;
              }
          }
          else{
            if(cenphi<141){ // range ..0~, cenphi -141~141
              if(noPhi-282-(cenphi)+width-1<0 and noPhi-cenphi+width>width*2){
                  fillindex=-2;
              }
            }
            else{ // range ~282.. num_jet==None, cenphi 141~282
              if(noPhi-(cenphi)+width<0 and noPhi-cenphi+282>width-1){
                  fillindex=-2;
              }
            }
          }
          if(fillindex!=-2){
            if(noTheta>=0){
                thetaindex=pix_per_tower-1-int(1.*segmentation->y(fiber.fiberNum)/ymax*pix_per_tower);
                phiindex=pix_per_tower-1-int(1.*segmentation->x(fiber.fiberNum)/xmax*pix_per_tower);
            }
            else{
                thetaindex=int(1.*segmentation->y(fiber.fiberNum)/ymax*pix_per_tower);
                phiindex=int(1.*segmentation->x(fiber.fiberNum)/xmax*pix_per_tower);
            }
            /*
            if(num_jet==0){
                thetaindex+=(noTheta+width-centheta)*pix_per_tower;
                if(noPhi>141){
                    phiindex+=(noPhi-282+width-1-cenphi)*pix_per_tower;
                }
                else{
                    phiindex+=(noPhi+width-cenphi)*pix_per_tower;
                }
            }
            if(num_jet==1){
                thetaindex+=(noTheta+width+1-centheta)*pix_per_tower;
                phiindex+=(noPhi-142+width-cenphi)*pix_per_tower;
            }
            */
            thetaindex+=(noTheta+width-centheta)*pix_per_tower;
            if(cenphi>width && cenphi+width<=282){ // range between 0~282
                phiindex+=(noPhi+width-cenphi)*pix_per_tower;
            }
            else{// range around including 282 and 0
                if(cenphi<141){ // range ..0~, cenphi 0~140
                    if(noPhi>141){ // left
                        phiindex+=(noPhi-282+width-1-cenphi)*pix_per_tower;
                    }
                    else{ // right
                        phiindex+=(noPhi+width-cenphi)*pix_per_tower;
                    }
                }
                else{ // range ~282.., cenphi 141~282
                    if(noPhi>141){ // left
                        phiindex+=(noPhi+width-cenphi)*pix_per_tower;
                    }
                    else{ // right
                        phiindex+=(noPhi+width-cenphi+282+1)*pix_per_tower;
                    }
                }
            }
            if(phiindex!=-1 && thetaindex!=-1 && fiber.Ecorr>=f_threshold){
              for(int k=0;k<fibercount+1;k++){
                if(k==2000)break;
                if(point_phi_idx_[k]==phiindex && point_theta_idx_[k]==thetaindex){
                  fillindex=k;
                  break;
                }
              }
              if(fillindex==-1){
                  fillindex=fibercount;
                  fibercount+=1;
              }
              if(fillindex<2000){
                    if(num_entry<2&&fibercount<10)printf("%d %d\n",fillindex,fibercount);
                    point_phi_idx_[fillindex]=phiindex;
                    point_theta_idx_[fillindex]=thetaindex;
                    if(segmentation->IsCerenkov(fiber.fiberNum)){
                      point_ecor_c_[fillindex]+=fiber.Ecorr; //cerenkov energy
										}
                    else{
                      point_depth_[fillindex]=(fiber.depth*fiber.Ecorr+point_ecor_s_[fillindex]*point_depth_[fillindex])/(point_ecor_s_[fillindex]+fiber.Ecorr);
                      point_ecor_s_[fillindex]+=fiber.Ecorr; //scintilation energy
                    }
              }
            }
          }
            phiindex=-1;
            thetaindex=-1;
            //std::sort(phis.begin(),phis.end());
            if(isjet==1 &&skipimg==0){
                if(idx_theta>=0){
                  if((segmentation->y(fiber.fiberNum))<28) thetaindex=1;
                  else thetaindex=0;
                  if((segmentation->x(fiber.fiberNum))<28) phiindex=1;
                  else phiindex=0;
                }
                else{
                  if((segmentation->y(fiber.fiberNum))<28) thetaindex=0;
                  else thetaindex=1;
                  if((segmentation->x(fiber.fiberNum))<28) phiindex=0;
                  else phiindex=1;
                }
                if(num_jet==0){
                  thetaindex=thetaindex+(idx_theta+40)*2;
                  if(idx_phi>142) phiindex=phiindex+(idx_phi-282+40-1)*2;
                  else phiindex=phiindex+(idx_phi+40)*2;
                }
                if(num_jet==1){//반대편 그림
                  thetaindex=thetaindex+(idx_theta+40+1)*2;
                  phiindex=phiindex+(idx_phi-142+40)*2;
                }
                if(!segmentation->IsCerenkov(fiber.fiberNum)){
                  cen_Ecorr_theta+=theta*fiber.Ecorr;
                  cen_Ecorr_phi+=phi*fiber.Ecorr;
                }
            }
            else{
              if(!segmentation->IsCerenkov(fiber.fiberNum)){
              cen_Ecorr_theta+=theta*fiber.Ecorr;
              cen_Ecorr_phi+=phi*fiber.Ecorr;
              }
                if(diff_theta<=2 && diff_theta>=0 && diff_phi<=2 && diff_phi>=0){
                  if(diff_theta>0){
                    phiindex = Int_t(55-segmentation->x(fiber.fiberNum));
                    thetaindex = Int_t(55-segmentation->y(fiber.fiberNum));
                  }
                  else{
                    phiindex = Int_t(segmentation->x(fiber.fiberNum));
                    thetaindex = Int_t(segmentation->y(fiber.fiberNum));
                  }
                }
                else{
                  phiindex=-1;
                  thetaindex=-1;
                }
            }

            if(phiindex!=-1 && thetaindex!=-1){
                if(isjet==1) imgindex=162*thetaindex+phiindex;//[theta,phi]
                else imgindex=3*phibin*etabin*diff_theta+3*phibin*thetaindex+phibin*diff_phi+phiindex;//[theta,phi]
                //imgindex=phibin*etabin*diff_phi+phibin*thetaindex+*phibin+phiindex;//[eta,phi]
                //imgindex=phibin*thetaindex+phiindex;//[eta,phi]
                if(imgindex>=28224 || imgindex<0 ){
                  printf("imgindex %d  theta %d phi %d towtheta %d towphi%d \n",imgindex,thetaindex,phiindex,idx_theta,idx_phi);
                  continue;
                }
                if(segmentation->IsCerenkov(fiber.fiberNum)){
                  image_ecor_c_[imgindex]+=fiber.Ecorr;
                  E_C_img+=fiber.Ecorr;
                  image_n_c_[imgindex]+=fiber.n;
                }
                else{
                  image_ecor_s_[imgindex]+=fiber.Ecorr;
                  E_Scorr_img+=fiber.Ecorr;
                  center_theta_img+=fiber.Ecorr*theta;
                  center_phi_img+=fiber.Ecorr*phi;
                  image_n_s_[imgindex]+=fiber.n;
                }
            }
        }

       }
    /*for(int k = 0; k<count; k++){
      if(ecimg[k]==E_C_img){
        same_ec+=1;
      }

      if(esimg[k]==E_Scorr_img){
        same_es+=1;
      }
    }*/
      ecimg[count]=E_C_img;
      esimg[count]=E_Scorr_img;
      center_theta_img=center_theta_img/E_Scorr_img;
      center_phi_img=center_phi_img/E_Scorr_img;
      cen_Ecorr_theta=cen_Ecorr_theta/recoEvt.E_Scorr;
      cen_Ecorr_phi=cen_Ecorr_phi/recoEvt.E_Scorr;
      }
       if(fibercheck==1){
         eventtree->Fill();
       }
       else{
         //printf("no fiber\n");
         nofiber+=1;
       }
     }
      //printf("entry end..\n");
   }
    drInterface->close();
    recoInterface->close();
    hepreader.close();

  }

//TFile outfile(outname.Data(), "recreate");
printf("%d files %d events filled in %s\n",readcount,int(eventtree->GetEntries()),outname.Data());
//printf("same_es %d same_ec %d\n",same_es,same_ec);
  //printf("writing..\n");
outfile->Write("");
//outfile->Write("",TObject::kOverwrite);
  //printf("closing..\n");
outfile->Close();
printf("--done\n");
if(nofiber+outrange+repeatfiber>0){
  printf("        unfilled faults %d  tower %d fiber %d\n",nofiber,outrange,repeatfiber);
}

return 0;
}
