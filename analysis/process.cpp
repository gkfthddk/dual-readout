#include "functions.h"//FIXME need for E_DR

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SparseVectorCollection.h"
#include "edm4hep/CalorimeterHit_TowerCollection.h"
//#include "DD4hep/DD4hepUnits.h"//FIXME

#include "podio/ROOTReader.h"
#include "podio/EventStore.h"


#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
//#include "Math/GenVector/Vector3D.h"
#include "TVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"

#include <iostream>
#include <cmath>
#include <string>

using namespace ROOT::Math;

class pixel {
private:
    int _numX[92] = {56, 56, 56, 56, 56, 56, 56, 55, 55, 55, 55, 55, 55, 55, 54, 54, 54, 54, 53, 53, 53, 53, 52, 52, 52, 52, 51, 51, 51, 50, 50, 50, 49, 49, 49, 48, 48, 48, 47, 47, 47, 46, 46, 46, 45, 45, 45, 44, 44, 44, 44, 43, 43, 43, 43, 43, 42, 42, 42, 42, 42, 42, 41, 41, 41, 41, 41, 41, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39};
    int _numY[92] = {56, 56, 56, 56, 56, 56, 56, 56, 56, 55, 55, 55, 55, 55, 55, 54, 54, 54, 54, 54, 53, 53, 53, 53, 52, 52, 52, 52, 51, 51, 51, 50, 50, 50, 49, 49, 49, 48, 48, 48, 47, 47, 47, 46, 46, 46, 46, 45, 45, 45, 44, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 21, 20, 19, 18, 17, 16, 15, 14, 14, 13, 12, 11, 10, 9, 8, 7};
    Int_t phi_idx;
    Int_t theta_idx;
    Int_t time_idx;
    Float_t minE_s=100.;
    Float_t minE_c=100.;
    Int_t minIdx_s=-1;
    Int_t minIdx_c=-1;
    Int_t count=0;
    int cen_phi=0;
    int cen_theta=0;
    int width = 560*2;
    int wing = 0;
    int edge_theta_left=-1;
    int edge_theta_right=1;
    int numx_edge(int side);
    int diff_theta_count(int diff_theta);
    int pixelize(int ix,int iy,int noPhi,int noTheta, float time);
    int check(short is_Cerenkov,float en);
public:
    int point_phi_idx[5000];
    int point_theta_idx[5000];
    int point_time_idx[5000];
    float point_e_s[5000];
    float point_e_c[5000];
    void Reset();
    int numX(int noTheta);
    int numY(int noTheta);
    void SetCenter(int _phi,int _theta);
    void SetWidth(int _width);
    void Add(int ix, int iy, int noPhi, int noTheta, short is_Cerenkov, float en, float time);

};

void pixel::Reset(){
  std::fill(point_e_s, point_e_s + 5000, 0.0);
  std::fill(point_e_c, point_e_c + 5000, 0.0);
  std::fill(point_phi_idx, point_phi_idx + 5000, 0);
  std::fill(point_theta_idx, point_theta_idx + 5000, 0);
  std::fill(point_time_idx, point_time_idx + 5000, 0);
  minE_s=100.;
  minE_c=100.;
  minIdx_s=-1;
  minIdx_c=-1;
  count=0;
}

int pixel::numX(int noTheta){
    if(noTheta<0) noTheta=-noTheta-1;
    return _numX[noTheta];
}
int pixel::numY(int noTheta){
    if(noTheta<0) noTheta=-noTheta-1;
    return _numY[noTheta];
}

int pixel::numx_edge(int side){
    int edge_count=0;
    for(int i=1; i<92; i++){
        int edge_theta=cen_theta+side*i;
        edge_count+=numX(edge_theta);
        if(edge_count>wing || edge_theta==side*91+int((-1+side)/2)) return edge_theta;
    }
    return cen_theta+side*91;
}

void pixel::SetCenter(int _phi,int _theta){
    cen_phi=_phi;
    cen_theta=_theta;
}
void pixel::SetWidth(int _width){
    width=_width;
    wing=int((width-numX(cen_theta))/2);
		edge_theta_left = numx_edge(-1);
		edge_theta_right = numx_edge(+1);
}

int pixel::diff_theta_count(int diff_theta){
    int ro_count=0;
    if(diff_theta==0) return ro_count;
    else if(diff_theta>0){
        for(int i=0; i<diff_theta; i++ ) ro_count+=numX(cen_theta+i);
        return ro_count;
    }
    else{
        for(int i=1; i<-diff_theta+1; i++) ro_count-=numX(cen_theta-i);
        return ro_count;
    }
}

int pixel::pixelize(int ix,int iy,int noPhi,int noTheta,float time){

		int diff_theta=noTheta-cen_theta;
		if(noTheta<edge_theta_left || noTheta>edge_theta_right)return -1;

		int width_phi = int(width/numY(noTheta));
		int diff_phi=noPhi-cen_phi;
    if(cen_phi-width_phi<0){//over left
			if(noPhi<cen_phi-width_phi+283 && diff_phi>width_phi) return -2;
	  	if(diff_phi>=-width_phi+283){
		  	diff_phi-=283;
		  }
		}
		else if(cen_phi+width_phi>282){//over right
			if(noPhi<cen_phi-width_phi && diff_phi>width_phi-283) return -3;
		  if(diff_phi<=width_phi-283){
		  	diff_phi+=283;
		  }
		}
		else{
			if(abs(diff_phi)>width_phi) return -4;
		}

		ix=ix-1;
		iy=iy-1;
		if(noTheta>=0){
			iy=numX(noTheta)-iy;
			ix=numY(noTheta)-ix;
		}

		phi_idx=diff_phi*numY(noTheta)+ix;
		theta_idx=diff_theta_count(diff_theta)+iy;

    time_idx=std::round(time*10);

    return 0;
}
int pixel::check(short is_Cerenkov,float en){
  int idxlimit=5000;
  minE_s=100.;
  minE_c=100.;
  for(int i=0; i<std::min(count,idxlimit); i++){
    if(count>idxlimit){
      if(minE_s>point_e_s[i] && point_e_c[i]==0.){
        minE_s=point_e_s[i];
        minIdx_s=i;
      }
      if(minE_c>point_e_c[i] && point_e_s[i]==0.){
        minE_c=point_e_c[i];
        minIdx_c=i;
      }
    }
    if(point_phi_idx[i]==phi_idx && point_theta_idx[i]==theta_idx && point_time_idx[i]==time_idx){
      count+=1;
      return i;
    }
  }
  if(count<idxlimit){
    count+=1;
    return count-1;
  }
  else{
    if(is_Cerenkov==1){
      if(en>minE_c){
        point_e_s[minIdx_c]=0.;
        point_e_c[minIdx_c]=0.;
        return minIdx_c;
      }
    }
    else{
      if(en>minE_s){
        point_e_s[minIdx_s]=0.;
        point_e_c[minIdx_s]=0.;
        return minIdx_s;
      }
    }
  }
  return -1;
}
void pixel::Add(int ix, int iy, int noPhi, int noTheta,short is_Cerenkov,float en, float time){
    int aa=pixelize(ix, iy, noPhi, noTheta, time);
    //std::cout<<aa<<" "<<phi_idx<<" "<<theta_idx<<std::endl;//FIXME
    if(aa!=0)return;
    //if(pixelize(ix, iy, noPhi, noTheta)!=0)return;
    int fillindex=-1;
    fillindex=int(check(is_Cerenkov,en));
    if(fillindex!=-1){
      point_phi_idx[fillindex]=phi_idx;
      point_theta_idx[fillindex]=theta_idx;
      point_time_idx[fillindex]=time_idx;
      if(is_Cerenkov==1){
        point_e_s[fillindex]=0.;
        point_e_c[fillindex]+=en;
      }
      else{
        point_e_s[fillindex]+=en;
        point_e_c[fillindex]=0.;
      }
    }
}

float getidx(float xmin,float xbinsize,float x){
  return int((x-xmin)/xbinsize);
}
void fillpoint(int phi_idx, int theta_idx, Int_t* point_phi, Int_t* point_theta,float en, Float_t* point_e_, int* count){

      int fillindex=-1;
      for(int k=0;k<*count+1;k++){
          if(k==2000)break;
          if(point_phi[k]==phi_idx && point_theta[k]==theta_idx){
            fillindex=k;
            break;
          }
      }
      if(fillindex==-1){
            fillindex= *count;
            *count+=1;
      }
      if(fillindex<2000){
              point_phi[fillindex]=phi_idx;
              point_theta[fillindex]=theta_idx;
              point_e_[fillindex]+=en; //cerenkov energy
      }
}

void filltower(int nophi, int notheta, Int_t* point_phi, Int_t* point_theta, int ix, int iy, Int_t* point_ix, Int_t* point_iy,float en, Float_t* point_e_, float* minE, int* minfill, int* count,Float_t* point_e_other){

      int fillindex=-1;
      for(int k=0;k<*count+1;k++){
          if(k==10000)break;
          if(point_phi[k]==nophi && point_theta[k]==notheta && point_ix[k]==ix && point_iy[k]==iy){
            if(point_e_other[k]!=0)continue;
            fillindex=k;
            break;
          }
          if(*count>=10000){
            if(point_e_[k]<*minE and point_e_other[k]<*minE){
              *minE=point_e_[k];
              *minfill=k;
            }
          }
      }
      if(fillindex==-1){
        if(*count<10000){
            fillindex= *count;
            *count+=1;
        }
        else if(*minE<en){
          fillindex=*minfill;
        }
      }
      if(fillindex<10000){
              point_phi[fillindex]=nophi;
              point_theta[fillindex]=notheta;
              point_ix[fillindex]=ix;
              point_iy[fillindex]=iy;
              point_e_[fillindex]+=en; //cerenkov energy
              point_e_other[fillindex]=0.;
              if(*count>=10000)*minE=100.;
      }
}

int main(int , char* argv[]) {
  TString filename = argv[1];

  auto pReader = std::make_unique<podio::ROOTReader>();
  pReader->openFile(static_cast<std::string>(filename));

  auto pStore = std::make_unique<podio::EventStore>();
  pStore->setReader(pReader.get());

  std::string filenameStd = static_cast<std::string>(filename);
  std::string extension = "_reco.root";
  auto where = filenameStd.find(extension);
  if (where != std::string::npos) {
    filenameStd.replace(where, extension.length(),"_img.root");
  }
  std::string extension2 = "box";
  auto where2 = filenameStd.find(extension2);
  if (where2 != std::string::npos) {
    filenameStd.replace(where2, extension2.length(),"imgs");
  }
  filename = static_cast<TString>(filenameStd);
  std::cout << "filename " << filename <<std::endl;
  TFile *outfile = new TFile(filename.Data(), "recreate");

  TTree * tree = new TTree("event","event info");

  //tree->SetAutoSave(0);
  std::vector<float> E_Ss,E_Cs;
  std::vector<int> leak_pdg,leak_genstat,leak_simstat;
  std::vector<float> leak_p,leak_px,leak_py,leak_pz;
  std::vector<float> leak_phi,leak_theta;

  std::vector<Float_t> ptc_p;
  Int_t num_tower=0;

  tree->Branch("ptc_p","vector<Float_t>",&ptc_p);
  tree->Branch("leak_p","vector<Float_t>",&leak_p);
  tree->Branch("leak_px","vector<Float_t>",&leak_px);
  tree->Branch("leak_py","vector<Float_t>",&leak_py);
  tree->Branch("leak_pz","vector<Float_t>",&leak_pz);
  tree->Branch("leak_phi","vector<Float_t>",&leak_phi);
  tree->Branch("leak_theta","vector<Float_t>",&leak_theta);
  tree->Branch("leak_genstat","vector<Int_t>",&leak_genstat);
  tree->Branch("leak_simstat","vector<Int_t>",&leak_simstat);
  tree->Branch("leak_pdg","vector<Int_t>",&leak_pdg);
  tree->Branch("num_tower",&num_tower,"ptd/I");

  #define BRANCH_I(name) tree->Branch(#name, & name, #name"/I");
  #define BRANCH_F(name) tree->Branch(#name, & name, #name"/F");
  #define BRANCH_A_(name, size, suffix) tree->Branch(#name, & name, #name"["#size"]/"#suffix);
  #define BRANCH_AF(name, size)  BRANCH_A_(name, size, F);
  #define BRANCH_AI(name, size)  BRANCH_A_(name, size, I);
  #define FILL_ZERO(array, size) std::fill(array, array + size, 0.0);
  #define FILL_ZEROI(array, size) std::fill(array, array + size, 0);

  float E_S=0.;
  float E_C=0.;
  float E_DR=0.;
  float E_DR291=0.;
  BRANCH_F(E_S);
  BRANCH_F(E_C);
  BRANCH_F(E_DR);
  BRANCH_F(E_DR291);

  int pdg=0;
  BRANCH_I(pdg);
  float genx=0;
  float geny=0;
  float genz=0;
  float genp=0;
  BRANCH_F(genx);
  BRANCH_F(geny);
  BRANCH_F(genz);
  BRANCH_F(genp);

  float genphi=0;
  float gentheta=0;
  BRANCH_F(genphi);
  BRANCH_F(gentheta);

  int count_reco3d_s=0;
  int count_reco3d_c=0;
  BRANCH_I(count_reco3d_s);
  BRANCH_I(count_reco3d_c);
  int count_reco2d_s=0;
  int count_reco2d_c=0;
  BRANCH_I(count_reco2d_s);
  BRANCH_I(count_reco2d_c);
  int count_sim2d_s=0;
  int count_sim2d_c=0;
  BRANCH_I(count_sim2d_s);
  BRANCH_I(count_sim2d_c);

  int cen_phi=0;
  int cen_theta=0;
  BRANCH_I(cen_phi);
  BRANCH_I(cen_theta);

  Float_t reco3d_phi[10000];//reco3d 2000
  Float_t reco3d_theta[10000];//reco3d 2000
  Float_t reco3d_x[10000];//reco3d 2000
  Float_t reco3d_y[10000];//reco3d 2000
  Float_t reco3d_z[10000];//reco3d 2000
  Int_t reco3d_phi_idx[2000];//reco3d 2000
  Int_t reco3d_theta_idx[2000];//reco3d 2000
  Float_t reco3d_e_s[10000];
  Float_t reco3d_e_c[10000];
  Float_t reco3d_e_s_idx[2000];
  Float_t reco3d_e_c_idx[2000];
  BRANCH_AF(reco3d_phi,10000);
  BRANCH_AF(reco3d_theta,10000);
  BRANCH_AF(reco3d_x,10000);
  BRANCH_AF(reco3d_y,10000);
  BRANCH_AF(reco3d_z,10000);
  BRANCH_AI(reco3d_phi_idx,2000);
  BRANCH_AI(reco3d_theta_idx,2000);
  BRANCH_AF(reco3d_e_s,10000);
  BRANCH_AF(reco3d_e_c,10000);
  BRANCH_AF(reco3d_e_s_idx,2000);
  BRANCH_AF(reco3d_e_c_idx,2000);
  FILL_ZERO(reco3d_phi,2000);
  FILL_ZERO(reco3d_theta,2000);
    FILL_ZERO(reco3d_x,10000);
    FILL_ZERO(reco3d_y,10000);
    FILL_ZERO(reco3d_z,10000);
  FILL_ZEROI(reco3d_phi_idx,2000);
  FILL_ZEROI(reco3d_theta_idx,2000);
  FILL_ZERO(reco3d_e_s,10000);
  FILL_ZERO(reco3d_e_c,10000);
  FILL_ZERO(reco3d_e_s_idx,2000);
  FILL_ZERO(reco3d_e_c_idx,2000);

  Int_t reco2d_phi_idx[2000];//reco2d 2000
  Int_t reco2d_theta_idx[2000];//reco2d 2000
  Float_t reco2d_e_s[2000];
  Float_t reco2d_e_c[2000];
  BRANCH_AI(reco2d_phi_idx,2000);
  BRANCH_AI(reco2d_theta_idx,2000);
  BRANCH_AF(reco2d_e_s,2000);
  BRANCH_AF(reco2d_e_c,2000);
  FILL_ZEROI(reco2d_phi_idx,2000);
  FILL_ZEROI(reco2d_theta_idx,2000);
  FILL_ZERO(reco2d_e_s,2000);
  FILL_ZERO(reco2d_e_c,2000);

  Int_t sim2d_phi_idx[2000];//sim2d 2000
  Int_t sim2d_theta_idx[2000];//sim2d 2000
  Float_t sim2d_e_s[2000];
  Float_t sim2d_e_c[2000];
  BRANCH_AI(sim2d_phi_idx,2000);
  BRANCH_AI(sim2d_theta_idx,2000);
  BRANCH_AF(sim2d_e_s,2000);
  BRANCH_AF(sim2d_e_c,2000);
  FILL_ZEROI(sim2d_phi_idx,2000);
  FILL_ZEROI(sim2d_theta_idx,2000);
  FILL_ZERO(sim2d_e_s,2000);
  FILL_ZERO(sim2d_e_c,2000);


  Int_t fiber_idx_phi[5000];
  Int_t fiber_idx_theta[5000];
  Int_t fiber_idx_time[5000];
  Float_t fiber_idx_e_s[5000];
  Float_t fiber_idx_e_c[5000];
  Int_t fiber_nophi[10000];
  Int_t fiber_notheta[10000];
  Int_t fiber_iy[10000];
  Int_t fiber_ix[10000];
  Float_t fiber_e_s[10000];
  Float_t fiber_e_c[10000];
  BRANCH_AI(fiber_idx_phi,5000);
  BRANCH_AI(fiber_idx_theta,5000);
  BRANCH_AI(fiber_idx_time,5000);
  BRANCH_AF(fiber_idx_e_s,5000);
  BRANCH_AF(fiber_idx_e_c,5000);
  BRANCH_AI(fiber_nophi,10000);
  BRANCH_AI(fiber_notheta,10000);
  BRANCH_AI(fiber_iy,10000);
  BRANCH_AI(fiber_ix,10000);
  BRANCH_AF(fiber_e_s,10000);
  BRANCH_AF(fiber_e_c,10000);
  FILL_ZEROI(fiber_nophi,10000);
  FILL_ZEROI(fiber_notheta,10000);
  FILL_ZEROI(fiber_iy,10000);
  FILL_ZEROI(fiber_ix,10000);
  FILL_ZERO(fiber_e_s,10000);
  FILL_ZERO(fiber_e_c,10000);
  Int_t tower_idx_phi[5000];
  Int_t tower_idx_theta[5000];
  Float_t tower_idx_e[5000];


  unsigned int entries = pReader->getEntries();
  for (unsigned int iEvt = 0; iEvt < entries; iEvt++) {
    //if (iEvt % 100 == 0) printf("Analyzing %dth event ...\n", iEvt);

    auto& genParticles = pStore->get<edm4hep::MCParticleCollection>("GenParticles");
    auto& edepHits = pStore->get<edm4hep::SimCalorimeterHitCollection>("SimCalorimeterHits");
    auto& edep3dHits = pStore->get<edm4hep::SimCalorimeterHitCollection>("Sim3dCalorimeterHits");
    auto& calo3dHits = pStore->get<edm4hep::CalorimeterHitCollection>("DRcalo3dHits");
    auto& calo2dTowers = pStore->get<edm4hep::CalorimeterHit_TowerCollection>("DRcalo2dTowers");
    auto& calo3dTowers = pStore->get<edm4hep::CalorimeterHit_TowerCollection>("DRcalo3dTowers");
    auto& Leakages = pStore->get<edm4hep::MCParticleCollection>("Leakages");

    auto& calo2dHits = pStore->get<edm4hep::CalorimeterHitCollection>("DRcalo2dHits");
    auto& rawTimeStructs = pStore->get<edm4hep::SparseVectorCollection>("RawTimeStructs");
    auto& rawWavlenStructs = pStore->get<edm4hep::SparseVectorCollection>("RawWavlenStructs");
    auto& digiWaveforms = pStore->get<edm4hep::SparseVectorCollection>("DigiWaveforms");
    auto& rawHits = pStore->get<edm4hep::RawCalorimeterHitCollection>("RawCalorimeterHits");
    auto& digiHits = pStore->get<edm4hep::RawCalorimeterHitCollection>("DigiCalorimeterHits");
    auto& procTimes = pStore->get<edm4hep::SparseVectorCollection>("DRpostprocTime");
/*std::cout << "entry " <<iEvt<<std::endl;
std::cout << "genParticles" << genParticles.size() << std::endl;
std::cout << "edepHits" << edepHits.size() << std::endl;
std::cout << "edep3dHits" << edep3dHits.size() << std::endl;
std::cout << "calo3dHits" << calo3dHits.size() << std::endl;

std::cout << "calo2dHits" << calo2dHits.size() << std::endl;
std::cout << "rawTimeStructs" << rawTimeStructs.size() << std::endl;
std::cout << "rawWavlenStructs" << rawWavlenStructs.size() << std::endl;
std::cout << "digiWaveforms" << digiWaveforms.size() << std::endl;
std::cout << "rawHits" << rawHits.size() << std::endl;
std::cout << "digiHits" << digiHits.size() << std::endl;
std::cout << "procTimes" << procTimes.size() << std::endl;*/
    float Edep = 0.;
    float width=0.033;
    int res=84;
    float phibin=width*2./res;
    float thetabin=width*2./res;
    genphi=0;
    gentheta=0;
    float en_threshold=0.0005;
    for (unsigned int igen = 0; igen < genParticles.size(); igen++){
      genParticles[igen].getCharge();
      genParticles[igen].getTime();
      genParticles[igen].getMass();
      auto vertex = genParticles[igen].getVertex();
      auto momentum = genParticles[igen].getMomentum();
      Polar3DVector v1;
      Polar3DVector v2;
      v1.SetXYZ(momentum.x,momentum.y,momentum.z);
      v2.SetXYZ(vertex.x,vertex.y,vertex.z);
      v1.SetR(1800);
      auto shift=v1+v2;
      genphi=shift.Phi();
      gentheta=shift.Theta();
      pdg=genParticles[igen].getPDG();
      genx=momentum.x;
      geny=momentum.y;
      genz=momentum.z;
      genp=sqrt(pow(genx,2)+pow(geny,2)+pow(genz,2));
    }
    for (unsigned int iEdep = 0; iEdep < edepHits.size(); iEdep++){
    }
    for (unsigned int iEdep = 0; iEdep < edep3dHits.size(); iEdep++){
      /*float en = caloHit.getEnergy();
      int nhits = rawHit.getAmplitude();
      (type==0) ? en_S += en : en_C += en;
      edep3dHits[iEdep].getEnergy
      edep3dHits[iEdep].getPosition
      edep3dHits[iEdep].getCellID
      edep3dHits[iEdep].getType*/
    }

    FILL_ZERO(reco3d_phi,10000);
    FILL_ZERO(reco3d_theta,10000);
    FILL_ZERO(reco3d_x,10000);
    FILL_ZERO(reco3d_y,10000);
    FILL_ZERO(reco3d_z,10000);
    FILL_ZEROI(reco3d_phi_idx,2000);
    FILL_ZEROI(reco3d_theta_idx,2000);
    FILL_ZERO(reco3d_e_s,10000);
    FILL_ZERO(reco3d_e_c,10000);
    FILL_ZERO(reco3d_e_s_idx,2000);
    FILL_ZERO(reco3d_e_c_idx,2000);

    FILL_ZEROI(reco2d_phi_idx,2000);
    FILL_ZEROI(reco2d_theta_idx,2000);
    FILL_ZERO(reco2d_e_s,2000);
    FILL_ZERO(reco2d_e_c,2000);

    FILL_ZEROI(sim2d_phi_idx,2000);
    FILL_ZEROI(sim2d_theta_idx,2000);
    FILL_ZERO(sim2d_e_s,2000);
    FILL_ZERO(sim2d_e_c,2000);
  FILL_ZEROI(fiber_idx_phi,5000);
  FILL_ZEROI(fiber_idx_theta,5000);
  FILL_ZEROI(fiber_idx_time,5000);
  FILL_ZERO(fiber_idx_e_s,5000);//change to fiber_e_s if class pixel works.
  FILL_ZERO(fiber_idx_e_c,5000);//change to fiber_e_s if class pixel works.
  FILL_ZEROI(fiber_nophi,10000);
  FILL_ZEROI(fiber_notheta,10000);
  FILL_ZEROI(fiber_iy,10000);
  FILL_ZEROI(fiber_ix,10000);
  FILL_ZERO(fiber_e_s,10000);
  FILL_ZERO(fiber_e_c,10000);
  FILL_ZEROI(tower_idx_phi,5000);
  FILL_ZEROI(tower_idx_theta,5000);
  FILL_ZERO(tower_idx_e,5000);
  leak_pdg.clear();
  leak_genstat.clear();
  leak_simstat.clear();
  leak_phi.clear();
  leak_theta.clear();
  leak_p.clear();
  leak_px.clear();
  leak_py.clear();
  leak_pz.clear();
  for (unsigned int idx = 0; idx < Leakages.size(); idx++) {
      const auto& leak = Leakages.at(idx);
      leak_pdg.push_back(leak.getPDG());
      auto momentum = leak.getMomentum();
      Polar3DVector v1;
      v1.SetXYZ(momentum.x,momentum.y,momentum.z);
      leak_p.push_back(sqrt(v1.Mag2()));
      leak_theta.push_back(v1.Theta());
      leak_phi.push_back(v1.Phi());
      leak_px.push_back(momentum.x);
      leak_py.push_back(momentum.y);
      leak_pz.push_back(momentum.z);
      leak_genstat.push_back(leak.getGeneratorStatus());
      leak_simstat.push_back(leak.getSimulatorStatus());
    }

    float thetamin=gentheta-width;
    float phimin=genphi-width;
    count_reco3d_s=0;
    count_reco3d_c=0;
    int fibercount=0;
    int fillindex=-1;
    int minindex_s=0;
    int minindex_c=0;
    float min_s=100;
    float min_c=100;
    float enS=0.;
    float enC=0.;
    int count_fiber=0;
    pixel draw;
    draw.Reset();
    cen_phi=0;
    cen_theta=0;
    float maxE=0.;
    int tower_count=0;
    for (unsigned int idx = 0; idx < calo2dTowers.size(); idx++) {
      const auto& caloTower = calo2dTowers.at(idx);
      if(caloTower.getType()!=0)continue;
      float en = caloTower.getEnergy();
      for (unsigned int i = 0; i<5000; i++){
        if(tower_idx_phi[i]==caloTower.getNoPhi() && tower_idx_theta[i]==caloTower.getNoTheta()){
          tower_idx_e[i]+=en;
          break;
        }
        if(tower_idx_e[i]==0. && tower_idx_phi[i]==0 && tower_idx_theta[i]==0){
          tower_idx_phi[i]=caloTower.getNoPhi();
          tower_idx_theta[i]=caloTower.getNoTheta();
          tower_idx_e[i]+=en;
          tower_count+=1;
          break;
        }
      }
    }
    for (unsigned int i = 0; i<5000; i++){
      if(tower_idx_e[i]>maxE){
        maxE=tower_idx_e[i];
        cen_phi=tower_idx_phi[i];
        cen_theta=tower_idx_theta[i];
      }
    }
    draw.SetCenter(cen_phi,cen_theta);
    draw.SetWidth(560*2);
    for (unsigned int idx = 0; idx < calo2dTowers.size(); idx++) {
      const auto& caloTower = calo2dTowers.at(idx);
      float en = caloTower.getEnergy();
      float time = caloTower.getTime();
      draw.Add( caloTower.getIx(),caloTower.getIy(),caloTower.getNoPhi(),caloTower.getNoTheta(),caloTower.getType(),en,time);
      if(en > en_threshold){
        if(caloTower.getType()==0) filltower(caloTower.getNoPhi(),caloTower.getNoTheta(), fiber_nophi,fiber_notheta, caloTower.getIx(),caloTower.getIy(), fiber_ix,fiber_iy, en, fiber_e_s,&min_s,&minindex_s,&count_fiber,fiber_e_c);
        else filltower(caloTower.getNoPhi(),caloTower.getNoTheta(), fiber_nophi,fiber_notheta, caloTower.getIx(),caloTower.getIy(), fiber_ix,fiber_iy, en, fiber_e_c,&min_c,&minindex_c,&count_fiber,fiber_e_s);
      /*
        if(fibercount<10000){
          fillindex=fibercount;
          if(caloTower.getType()==0){
            fiber_e_s[fillindex]=en;
            fiber_e_c[fillindex]=0.;
            if(en<min_s){
              min_s=en;
              minindex_s=fillindex;
            }
          }
          else{
            fiber_e_s[fillindex]=0.;
            fiber_e_c[fillindex]=en;
            if(en<min_c){
              min_c=en;
              minindex_c=fillindex;
            }
          }
        }
        else{
          if(caloTower.getType()==0){
            if(min_s<en){
              fillindex=minindex_s;
              fiber_e_s[fillindex]=en;
              fiber_e_c[fillindex]=0.;
              min_s=100;
              for( unsigned int k = 0; k < 10000; k++){
                if(fiber_e_s[k]<min_s){
                  min_s=fiber_e_s[k];
                  minindex_s=k;
                }
              }
            }
            else{
              fillindex=-1;
            }
          }
          else{
            if(min_c<en){
              fillindex=minindex_c;
              fiber_e_s[fillindex]=0.;
              fiber_e_c[fillindex]=en;
              min_c=100;
              for( unsigned int k = 0; k < 10000; k++){
                if(fiber_e_c[k]<min_c){
                  min_c=fiber_e_c[k];
                  minindex_c=k;
                }
              }
            }
            else{
              fillindex=-1;
            }
          }
        }
        if(fillindex!=-1){
          fiber_nophi[fillindex]=caloTower.getNoPhi();
          fiber_notheta[fillindex]=caloTower.getNoTheta();
          fiber_iy[fillindex]=caloTower.getIy();
          fiber_ix[fillindex]=caloTower.getIx();
        }
        fibercount+=1;
        */
      }
    }
    for (unsigned int idx = 0; idx<5000; idx++){
      fiber_idx_phi[idx]=draw.point_phi_idx[idx];
      fiber_idx_theta[idx]=draw.point_theta_idx[idx];
      fiber_idx_time[idx]=draw.point_time_idx[idx];
      fiber_idx_e_s[idx]=draw.point_e_s[idx];
      fiber_idx_e_c[idx]=draw.point_e_c[idx];
    }
    fibercount=0;
    fillindex=-1;
    minindex_s=0;
    minindex_c=0;
    min_s=100;
    min_c=100;
    enS=0.;
    enC=0.;
    for (unsigned int idx = 0; idx < calo3dHits.size(); idx++) {
      const auto& caloHit = calo3dHits.at(idx);
      float en = caloHit.getEnergy();
      if(en > en_threshold){
        auto position = caloHit.getPosition();
        Polar3DVector v1;
        v1.SetXYZ(position.x,position.y,position.z);
        int phi_idx = getidx(phimin, phibin,v1.Phi());
        int theta_idx = getidx(thetamin, thetabin,v1.Theta());
        if(caloHit.getType()==0) fillpoint(phi_idx, theta_idx, reco3d_phi_idx, reco3d_theta_idx, en, reco3d_e_s_idx, &count_reco3d_s);
        else fillpoint(phi_idx, theta_idx, reco3d_phi_idx, reco3d_theta_idx, en, reco3d_e_c_idx, &count_reco3d_c);
        if(fibercount<10000){
          fillindex=fibercount;
          if(caloHit.getType()==0){
            reco3d_e_s[fillindex]=en;
            reco3d_e_c[fillindex]=0.;
            if(en<min_s){
              min_s=en;
              minindex_s=fillindex;
            }
          }
          else{
            reco3d_e_s[fillindex]=0.;
            reco3d_e_c[fillindex]=en;
            if(en<min_c){
              min_c=en;
              minindex_c=fillindex;
            }
          }
        }
        else{
          if(caloHit.getType()==0){
            if(min_s<en){
              fillindex=minindex_s;
              reco3d_e_s[fillindex]=en;
              reco3d_e_c[fillindex]=0.;
              min_s=100;
              for( unsigned int k = 0; k < 10000; k++){
                if(reco3d_e_s[k]<min_s){
                  min_s=reco3d_e_s[k];
                  minindex_s=k;
                }
              }
            }
            else{
              fillindex=-1;
            }
          }
          else{
            if(min_c<en){
              fillindex=minindex_c;
              reco3d_e_s[fillindex]=0.;
              reco3d_e_c[fillindex]=en;
              min_c=100;
              for( unsigned int k = 0; k < 10000; k++){
                if(reco3d_e_c[k]<min_c){
                  min_c=reco3d_e_c[k];
                  minindex_c=k;
                }
              }
            }
            else{
              fillindex=-1;
            }
          }
        }
        if(fillindex!=-1){
          reco3d_phi[fillindex]=v1.Phi();
          reco3d_theta[fillindex]=v1.Theta();
          reco3d_x[fillindex]=position.x;
          reco3d_y[fillindex]=position.y;
          reco3d_z[fillindex]=position.z;
        }
        fibercount+=1;
      /*int fillindex=-1;
        if(caloHit.getType()==0) fibercount=count_reco3d_s;
        else fibercount=count_reco3d_c;

        for(int k=0;k<fibercount+1;k++){
          if(k==2000)break;
          if(point_phi_idx[k]==phi_idx && point_theta_idx[k]==theta_idx){
            fillindex=k;
            break;
          }
        }
        if(fillindex==-1){
            fillindex=fibercount;
            fibercount+=1;
        }
        if(fillindex<2000){
              point_phi_idx[fillindex]=phi_idx;
              point_theta_idx[fillindex]=theta_idx;
              if(caloHit.getType()==0){
                point_e_s[fillindex]+=en; //scintilation energy
                count_reco3d_s=fibercount;
              }
              else{
                point_e_c[fillindex]+=en; //cerenkov energy
                count_reco3d_c=fibercount;
              }
        }*/
      }
      if(caloHit.getType()==0)enS+=en;
      else enC+=en;
    }
    float enS_cut=0.;
    float enC_cut=0.;
    for(int i=0;i<2000;i++){
      enS_cut+=reco3d_e_s[i];
      enC_cut+=reco3d_e_c[i];
    }

    count_reco2d_s=0;
    count_reco2d_c=0;
    E_S=0.;
    E_C=0.;
    for (unsigned int idx = 0; idx < calo2dHits.size(); idx++) {
      const auto& caloHit = calo2dHits.at(idx);
      float en = caloHit.getEnergy();
      if(en > en_threshold){
        auto position = caloHit.getPosition();
        Polar3DVector v1;
        v1.SetXYZ(position.x,position.y,position.z);
        int phi_idx = getidx(phimin, phibin,v1.Phi());
        int theta_idx = getidx(thetamin, thetabin,v1.Theta());
        if(caloHit.getType()==0) fillpoint(phi_idx, theta_idx, reco2d_phi_idx, reco2d_theta_idx, en, reco2d_e_s, &count_reco2d_s);
        else fillpoint(phi_idx, theta_idx, reco2d_phi_idx, reco2d_theta_idx, en, reco2d_e_c, &count_reco2d_c);
      }
      if(caloHit.getType()==0)E_S+=en;
      else E_C+=en;
    }
    E_DR=functions::E_DR(E_C,E_S);
    E_DR291=functions::E_DR291(E_C,E_S);
    count_sim2d_s=0;
    count_sim2d_c=0;
    for (unsigned int idx = 0; idx < edepHits.size(); idx++) {
      const auto& caloHit = edepHits.at(idx);
      float en = caloHit.getEnergy();
      if(en > en_threshold){
        auto position = caloHit.getPosition();
        Polar3DVector v1;
        v1.SetXYZ(position.x,position.y,position.z);
        int phi_idx = getidx(phimin, phibin,v1.Phi());
        int theta_idx = getidx(thetamin, thetabin,v1.Theta());
        fillpoint(phi_idx, theta_idx, sim2d_phi_idx, sim2d_theta_idx, en, sim2d_e_s, &count_sim2d_s);
        //if(caloHit.getType()==0) fillpoint(phi_idx, theta_idx, sim2d_phi_idx, sim2d_theta_idx, en, sim2d_e_s, &count_sim2d_s);FIXME
        //else fillpoint(phi_idx, theta_idx, sim2d_phi_idx, sim2d_theta_idx, en, sim2d_e_c, &count_sim2d_c);
      }
    }
    //printf("fiber %d %d\n",count_reco3d_s,count_reco3d_c);
    //std::cout<<"enS "<<enS_cut<< " "<<enS<<std::endl;
    //std::cout<<"enC "<<enC_cut<< " "<<enC<<std::endl;
    tree->Fill();

    /*
    for (unsigned int idx = 0; idx < calo2dHits.size(); idx++) {
      const auto& caloHit = calo2dHits.at(idx);
      const auto& timeStruct = rawTimeStructs.at(idx);
      const auto& wavlenStruct = rawWavlenStructs.at(idx);
      const auto& waveform = digiWaveforms.at(idx);
      const auto& rawHit = rawHits.at(idx);
      const auto& digiHit = digiHits.at(idx);
      const auto& procTime = procTimes.at(idx);
      caloHit.getPosition();
      caloHit.getTime();
      int type = caloHit.getType();
      float en = caloHit.getEnergy();
      int nhits = rawHit.getAmplitude();
      
      (type==0) ? en_S += en : en_C += en;

      tNhit->Fill(nhits);
      tToa->Fill(caloHit.getTime());

      if ( digiHit.getAmplitude() > 0 )
        tInt->Fill(static_cast<double>(digiHit.getAmplitude()));

      for (unsigned int bin = 0; bin < timeStruct.centers_size(); bin++) {
        float timeBin = timeStruct.getCenters(bin);
        tT->Fill(timeBin,timeStruct.getContents(bin));
      }

      for (unsigned int bin = 0; bin < wavlenStruct.centers_size(); bin++) {
        float wavlenBin = wavlenStruct.getCenters(bin);
        tWav->Fill(wavlenBin,wavlenStruct.getContents(bin));
      }

      for (unsigned int bin = 0; bin < waveform.centers_size(); bin++) {
        float timeBin = waveform.getCenters(bin);
        tD->Fill(timeBin,waveform.getContents(bin));
      }

      for (unsigned int bin = 0; bin < procTime.centers_size(); bin++) {
        float timeBin = procTime.getCenters(bin);
        tProc->Fill(timeBin,procTime.getContents(bin));
      }
    }

    E_Ss.push_back(en_S);
    E_Cs.push_back(en_C);*/

    pStore->clear();
    pReader->endOfEvent();
  } // event loop



outfile->Write("");
  //printf("closing..\n");
outfile->Close();
}
