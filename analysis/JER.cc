#include "RootInterface.h"
#include "RecoInterface.h"
#include "DRsimInterface.h"
#include "fastjetInterface.h"
#include "functions.h"

#include "fastjet/PseudoJet.hh"

#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TGraph.h"

#include <iostream>
#include <string>
#include <cmath>

int main(int argc, char* argv[]) {
  TString filename = argv[1];
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);
  float cen = std::stof(argv[4]);


  HepMC3::ReaderRootTree reader(std::string(filename)+".root");

  // fastjetInterface fjTower_S;
  // fjTower_S.set(recoInterface->getTree(),"RecoTowerJets_S");
  // fastjetInterface fjFiber_S;
  // fjFiber_S.set(recoInterface->getTree(),"RecoFiberJets_S");
  // fastjetInterface fjFiber_C;
  // fjFiber_C.set(recoInterface->getTree(),"RecoFiberJets_C");
  // fastjetInterface fjGen;
  // fjGen.set(reader.m_tree,"GenJets");
  //unsigned int entries = reader->entries();
  //unsigned int entries = recoInterface->entries();
  std::vector<Float_t> gen_energy;
    float Etot = 0.;
  TFile outfile("test.root", "recreate");
  TTree eventtree("event","event info");
  eventtree.Branch("gen_energy","vector<Float_t>",&gen_energy);
  eventtree.Branch("Etot",&Etot,"Etot/F");
  for (int num=0; num<10;num++) {

    HepMC3::GenEvent genEvt;
    reader.read_event(genEvt);
    gen_energy.clear();

    Etot = 0.;
    for (auto ptc : genEvt.particles()) {
      int abspid = std::abs(ptc->pid());
      if( ptc->status() == 23){
        gen_energy.push_back(ptc->momentum().e());

      }
      if ( ptc->status() != 1 ) continue;
      if ( abspid == 12 || abspid == 14 || abspid == 16 ) continue;

      auto mom = ptc->momentum();

      Etot += mom.e();
    }
    eventtree.Fill();
  } // event loop

  reader.close();
outfile.Write("");
outfile.Close();
}
