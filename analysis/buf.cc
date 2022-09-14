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
  new GeoSvc({"compact/DRcalo.xml"});

  auto m_geoSvc = GeoSvc::GetInstance();
  std::string m_readoutName = "DRcaloSiPMreadout";

  auto lcdd = m_geoSvc->lcdd();
  auto allReadouts = lcdd->readouts();
  if (allReadouts.find(m_readoutName) == allReadouts.end()) {
    throw std::runtime_error("Readout " + m_readoutName + " not found! Please check tool configuration.");
  } else {
    std::cout << "Reading EDM from the collection " << m_readoutName << std::endl;
  }
  TVector3* towerxyz = new TVector3();

      auto segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDRcalo*>(m_geoSvc->lcdd()->readout(m_readoutName).segmentation().segmentation());
  std::ofstream writeFile;
  writeFile.open("towerposition.txt");
  char buffer[200];
  int buf_size=0;
  buf_size=sprintf(buffer,"nophi notheta towerphi towertheta towerpos.x towerpos.y towerpos.z\n");
  writeFile.write(buffer, buf_size);
  for(int i =0;i<283;i++){
    for(int j =-91;j<92;j++){
      printf("nophi %d notheta %d\n",i,j);

      auto towerpos = segmentation->towerposition(j,i);
      towerxyz->SetXYZ(towerpos.x(),towerpos.y(),towerpos.z());
      auto towerphi=float(towerxyz->Phi());
      auto towertheta=float(towerxyz->Theta());
      buf_size=sprintf(buffer,"%d %d %g %g %g %g %g\n",i,j,towerphi,towertheta,towerpos.x(),towerpos.y(),towerpos.z());
      writeFile.write(buffer, buf_size);
  }
  }
    //printf("load %s...;\n",inname.Data());
    //mychain.Add(inname+"/DRsim");
writeFile.close();

return 0;
}
