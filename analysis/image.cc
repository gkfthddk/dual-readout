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

//./image ../../tools/pack/gg_70GeV_2.root ../../gg_2 1
int main(int argc, char* argv[]){
  printf("begin;\n");
  TString innameform = argv[1];
  TString filename = argv[2];  
  int isjet = atoi(argv[3]);
  TString reconameform;
  TString inname;
  TString reconame;
  TList* keys;
  int readcount=0;
  int check_event=0;
  
  float pi=TMath::Pi();
  TString outname=filename;
  TFile infile(innameform.Data(), "read");
  keys = infile.GetListOfKeys();
  if(keys->GetEntries()>1){
      for(int i=0;i<keys->GetEntries();i++){
        if(strcmp(keys->At(i)->GetName(),"event")==0)check_event=1;
      }
  }
  if(check_event==0){
    printf("no file found..\n");
    return(1);
  }
  TFile outfile(outname.Data(), "recreate");
  std::vector<Int_t> * tower_idx_eta=0;
  std::vector<Int_t> * tower_idx_phi=0;
  std::vector<Float_t> * tower_eta=0;
  std::vector<Float_t> * tower_theta=0;
  std::vector<Float_t> * tower_phi=0;
  std::vector<Int_t> * tower_notheta=0;
  std::vector<Int_t> * tower_nophi=0;
  float mass_Gen=0.;
  float center_phi_gen=0;
  float center_theta_gen=0;
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
  Int_t num_jet=0;
  Int_t num_entry=-1;
  Int_t not_phi=0;
  Int_t not_theta=0;
  std::vector<Float_t> * fiber_energy=0;
  std::vector<Float_t> * fiber_ecor=0;
  std::vector<Float_t> * fiber_ecor_s=0;
  std::vector<Float_t> * fiber_ecor_c=0;
  std::vector<Int_t> * fiber_n=0;
  std::vector<Int_t> * fiber_itower=0;
  std::vector<Int_t> * fiber_ix=0;
  std::vector<Int_t> * fiber_iy=0;
  std::vector<Float_t> * fiber_t=0;
  std::vector<Float_t> * fiber_x=0;
  std::vector<Float_t> * fiber_y=0;
  std::vector<Float_t> * fiber_z=0;
  std::vector<Float_t> * fiber_depth=0;
  std::vector<Float_t> * fiber_r=0;
  std::vector<Float_t> * fiber_theta=0;
  std::vector<Float_t> * fiber_phi=0;
  std::vector<Float_t> * fiber_eta=0;
  std::vector<Bool_t> * fiber_iscerenkov=0;

  TTree imagetree("event","event info");
  TTree *eventtree = (TTree *) infile.Get("event");
  eventtree->SetBranchAddress("tower_eta",&tower_eta);
  //eventtree->SetAutoSave(0);
  eventtree->SetBranchAddress("tower_notheta",&tower_notheta);
  eventtree->SetBranchAddress("tower_nophi",&tower_nophi);
  eventtree->SetBranchAddress("tower_idx_eta",&tower_idx_eta);
  eventtree->SetBranchAddress("tower_idx_phi",&tower_idx_phi);
  eventtree->SetBranchAddress("num_jet",&num_jet);
  eventtree->SetBranchAddress("num_entry",&num_entry);
  eventtree->SetBranchAddress("not_theta",&not_theta);
  eventtree->SetBranchAddress("not_phi",&not_phi);
  eventtree->SetBranchAddress("E_Gen",&E_Gen);
  eventtree->SetBranchAddress("mass_Gen",&mass_Gen);
  eventtree->SetBranchAddress("pt_Gen",&pt_Gen);
  eventtree->SetBranchAddress("p_Gen",&p_Gen);
  eventtree->SetBranchAddress("Edep",&Edep);
  eventtree->SetBranchAddress("center_phi_gen",&center_phi_gen);
  eventtree->SetBranchAddress("center_theta_gen",&center_theta_gen);
  eventtree->SetBranchAddress("E_C",&E_C);
  eventtree->SetBranchAddress("E_S",&E_S);
  eventtree->SetBranchAddress("E_Scorr",&E_Scorr);
  eventtree->SetBranchAddress("n_C",&n_C);
  eventtree->SetBranchAddress("n_S",&n_S);
  eventtree->SetBranchAddress("E_DR",&E_DR);
  eventtree->SetBranchAddress("E_DRcorr",&E_DRcorr);
  eventtree->SetBranchAddress("Eleak_nu",&Eleak_nu);
  eventtree->SetBranchAddress("Pleak",&Pleak);
  eventtree->SetBranchAddress("fiber_energy",&fiber_energy);
  eventtree->SetBranchAddress("fiber_ecor",&fiber_ecor);
  eventtree->SetBranchAddress("fiber_ecor_s",&fiber_ecor_s);
  eventtree->SetBranchAddress("fiber_ecor_c",&fiber_ecor_c);
  eventtree->SetBranchAddress("fiber_n",&fiber_n);
  eventtree->SetBranchAddress("fiber_itower",&fiber_itower);
  eventtree->SetBranchAddress("fiber_ix",&fiber_ix);
  eventtree->SetBranchAddress("fiber_iy",&fiber_iy);
  eventtree->SetBranchAddress("fiber_t",&fiber_t);
  eventtree->SetBranchAddress("fiber_x",&fiber_x);
  eventtree->SetBranchAddress("fiber_y",&fiber_y);
  eventtree->SetBranchAddress("fiber_z",&fiber_z);
  eventtree->SetBranchAddress("fiber_depth",&fiber_depth);
  eventtree->SetBranchAddress("fiber_r",&fiber_r);
  eventtree->SetBranchAddress("fiber_theta",&fiber_theta);
  eventtree->SetBranchAddress("fiber_phi",&fiber_phi);
  eventtree->SetBranchAddress("fiber_eta",&fiber_eta);
  eventtree->SetBranchAddress("fiber_iscerenkov",&fiber_iscerenkov);

  //Float_t voxel_ecor_s_[729000];//90*90*90
  //Int_t voxel_n_s_[729000];
  Float_t image_ecor_c_[28224];//168*168
  Int_t image_n_c_[28224];
  Float_t image_ecor_s_[28224];//168*168
  Int_t image_n_s_[28224];
  Float_t point_ecor_s_[2000];//point 2000
    Float_t point_ecor_c_[2000];//point 2000
      Int_t point_phi_idx_[2000];//point 2000
        Int_t point_theta_idx_[2000];//point 2000
          Float_t point_depth_[2000];//point 2000
  Float_t point_2048_[8192];//90*90
  Float_t cen_Ecorr_theta=0.;
  Float_t cen_Ecorr_phi=0.;
  imagetree.Branch("cen_Ecorr_theta",&cen_Ecorr_theta,"cen_Ecorr_theta/F");
  imagetree.Branch("cen_Ecorr_phi",&cen_Ecorr_phi,"cen_Ecorr_phi/F");
  #define BRANCH_A_(name, size, suffix) imagetree.Branch(#name, & name##_, #name"["#size"]/"#suffix);
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
  BRANCH_AF(point_ecor_s,2000);
    BRANCH_AF(point_ecor_c,2000);
      BRANCH_AI(point_phi_idx,2000);
        BRANCH_AI(point_theta_idx,2000);
          BRANCH_AF(point_depth,2000);
  BRANCH_AI(point_2048, 8192);//don't forget to add FILL_ZERO in loop
  //FILL_ZERO(voxel_ecor_s_,729000);
  //FILL_ZEROI(voxel_n_s_,729000);
  FILL_ZERO(image_ecor_c_,28224);
  FILL_ZEROI(image_n_c_,28224);
  FILL_ZERO(image_ecor_s_,28224);
  FILL_ZEROI(image_n_s_,28224);
  FILL_ZEROI(point_2048_,8192);
  int width=40;
  int binsize=2*2*40+2;
  int maxindex=binsize*binsize;
  int phiindex=-1;
  int thetaindex=-1;
  int imgindex=-1;
  int noPhi=0;
  int noTheta=0;
  int fibercount=0;
  bool isopposite = false;
  int skipimg=0;
    int fibercheck=0;
  Int_t outrange=0;
  Int_t repeatfiber=0;
  Int_t nofiber=0;
int ymax=56;
int xmax=56;
int pix_per_tower=56;
  float f_threshold=0.001;

  for (int ent = 0; ent<eventtree->GetEntries(); ent++) {  
    eventtree->GetEntry(ent);
    process_class process;
    int cenphi=process.findphi(center_phi_gen+num_jet*pi);
    int centheta=process.findtheta(center_theta_gen);
    if(ent==11)printf("phi %d; theta %d %f %f %d\n",cenphi,centheta,center_phi_gen,center_theta_gen,num_jet);
    //if(count%int(entries/10.)==0)printf("%d%%--\n",int(100.*count/entries));
    fibercount=0; 

    fibercheck=0;
      cen_Ecorr_theta=0.;
      cen_Ecorr_phi=0.;
      FILL_ZERO(image_ecor_c_,28224);
      FILL_ZEROI(image_n_c_,28224);
      FILL_ZERO(image_ecor_s_,28224);
      FILL_ZEROI(image_n_s_,28224);
      FILL_ZERO(point_ecor_s_,2000);
      FILL_ZERO(point_ecor_c_,2000);
      FILL_ZERO(point_phi_idx_,2000);
      FILL_ZERO(point_theta_idx_,2000);
      FILL_ZERO(point_depth_,2000);
      if(num_jet==0)isopposite=false;
      if(num_jet==1)isopposite=true;
      E_Scorr_img=0.;
      E_C_img=0.;
      for(int i = 0; i<int(fiber_n->size()); i++){
        skipimg=0;
        if(fiber_ecor->at(i)==0)continue;
        noTheta=tower_idx_eta->at(i);
        noPhi=tower_idx_phi->at(i);
        if(isopposite==false){
          if(noTheta+width<0 || noTheta+width>80) skipimg=1;
          if(noPhi-282+width-1<0 && noPhi+width>80)skipimg=1;//phi 0 방향일때만
        }
        else{
          if(noTheta+width+1<0 || noTheta+width+1>80) skipimg=1;
          if(noPhi-142+width<0 || noPhi-142+width>80)skipimg=1;//phi pi 방향일때만
        }
            phiindex=-1;
            thetaindex=-1;
        if(noTheta<centheta-width-0 or noTheta>centheta+width-0){
                continue;
        }
        if(cenphi>width and cenphi+width<=282){ // range between 0~282
            if(noPhi<cenphi-width or noPhi>cenphi+width){
                continue;
            }
        }
        else{
          if(cenphi<141){ // range ..0~, cenphi -141~141
            if(noPhi-282-(cenphi)+width-1<0 and noPhi-cenphi+width>width*2){
              continue;
            }
          }
          else{ // range ~282.. num_jet==None, cenphi 141~282
            if(noPhi-(cenphi)+width<0 and noPhi-cenphi+282>width-1){
              continue;
            }
          }
        }
        if(noTheta>=0){
            thetaindex=pix_per_tower-1-int(1.*fiber_iy->at(i)/ymax*pix_per_tower);
            phiindex=pix_per_tower-1-int(1.*fiber_ix->at(i)/xmax*pix_per_tower);
        }
        else{
            thetaindex=int(1.*fiber_iy->at(i)/ymax*pix_per_tower);
            phiindex=int(1.*fiber_ix->at(i)/xmax*pix_per_tower);
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
            if(phiindex!=-1 && thetaindex!=-1 && fiber_ecor->at(i)>f_threshold){
							int fillindex=-1;
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
                  if(ent==11&&fibercount<15)printf("%d %d; iy %d ix %d; phi %d theta %d; %d %d\n",ent,i,fiber_iy->at(i),fiber_ix->at(i),phiindex,thetaindex,fillindex,fibercount);
                  //if(ent==11)printf(" %d %d %f %f;\n",point_phi_idx_[0],point_theta_idx_[0],point_ecor_s_[0],point_ecor_c_[0]);
                  point_phi_idx_[fillindex]=phiindex;
                  point_theta_idx_[fillindex]=thetaindex;
                  if(fiber_iscerenkov->at(i)){
                    point_ecor_c_[fillindex]+=fiber_ecor->at(i); //cerenkov energy
							    }
                  else{
                    point_depth_[fillindex]=(fiber_depth->at(i)*fiber_ecor->at(i)+point_ecor_s_[fillindex]*point_depth_[fillindex])/(point_ecor_s_[fillindex]+fiber_ecor->at(i));
                    point_ecor_s_[fillindex]+=fiber_ecor->at(i);
                  }
              }
            }
            fibercheck=1;
            phiindex=-1;
            thetaindex=-1;
            imgindex=-1;
            
            if(skipimg==0){
                if(noTheta>=0){
                  if(fiber_iy->at(i)<28) thetaindex=1;
                  else thetaindex=0;
                  if((fiber_ix->at(i))<28) phiindex=1;
                  else phiindex=0;
                }
                else{
                  if((fiber_iy->at(i))<28) thetaindex=0;
                  else thetaindex=1;
                  if(fiber_ix->at(i)<28) phiindex=0;
                  else phiindex=1;
                }
                if(num_jet==0){
                  thetaindex=thetaindex+(noTheta+width)*2;
                  if(noPhi>142) phiindex=phiindex+(noPhi-282+width-1)*2;
                  else phiindex=phiindex+(noPhi+width)*2;
                }
                if(num_jet==1){//반대편 그림
                  thetaindex=thetaindex+(noTheta+width+1)*2;
                  phiindex=phiindex+(noPhi-142+width)*2;
                }
                //if(!segmentation->IsCerenkov(fiber.fiberNum)){
                //  cen_Ecorr_theta+=theta*fiber.Ecorr;
                //  cen_Ecorr_phi+=phi*fiber.Ecorr;
                //}
            }

            if(phiindex!=-1 && thetaindex!=-1){
                imgindex=binsize*thetaindex+phiindex;//[theta,phi]
                if(imgindex>=maxindex || imgindex<0 ){
                  printf("imgindex %d  theta %d phi %d towtheta %d towphi%d \n",imgindex,thetaindex,phiindex,idx_theta,idx_phi);
                  continue;
                }
                if(fiber_iscerenkov->at(i)){
                  image_ecor_c_[imgindex]+=fiber_ecor->at(i);
                  E_C_img+=fiber_ecor->at(i);
                  image_n_c_[imgindex]+=fiber_n->at(i);
                }
                else{
                  image_ecor_s_[imgindex]+=fiber_ecor->at(i);
                  E_Scorr_img+=fiber_ecor->at(i);
                  image_n_s_[imgindex]+=fiber_n->at(i);
                }
            }

        
      //center_theta_img=center_theta_img/E_Scorr_img;
      //center_phi_img=center_phi_img/E_Scorr_img;
      //cen_Ecorr_theta=cen_Ecorr_theta/recoEvt.E_Scorr;
      //cen_Ecorr_phi=cen_Ecorr_phi/recoEvt.E_Scorr;
      }
       if(fibercheck==1){
         imagetree.Fill();
       }
       else{
         //printf("no fiber\n");
         nofiber+=1;
       }
       //if(ent==15)break;

   }
    

//TFile outfile(outname.Data(), "recreate");
  printf("writing..\n");
outfile.Write("");
//outfile.Write("",TObject::kOverwrite);
  printf("closing..\n");
outfile.Close();
printf("--done\n");
printf("%d files %d events filled in %s\n",readcount,int(imagetree.GetEntries()),outname.Data());
printf("        unfilled faults %d  tower %d fiber %d\n",nofiber,outrange,repeatfiber);

return 0;
}
