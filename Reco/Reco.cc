#include "DRsimInterface.h"
#include "RootInterface.h"
#include "fastjetInterface.h"
#include "RecoTower.h"

#include "GeoSvc.h"

#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
  std::string filenum = std::string(argv[1]);
  std::string filename = std::string(argv[2]);
  std::string filename2 = std::string(argv[2]);
  if (argc > 3) filename2 = std::string(argv[3]);

  RootInterface<RecoInterface::RecoEventData>* recoInterface = new RootInterface<RecoInterface::RecoEventData>(filename2+"_"+filenum+"_Reco.root");
  recoInterface->create("Reco","RecoEventData");

  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(filename+"_"+filenum+".root");
  drInterface->set("DRsim","DRsimEventData");

  new GeoSvc({"./bin/compact/DRcalo.xml"});

  RecoTower* recoTower = new RecoTower();
  recoTower->readCSV();

  unsigned int entries = drInterface->entries();
  while (drInterface->numEvt() < entries) {
    recoTower->getFiber()->clear();

    DRsimInterface::DRsimEventData evt;
    RecoInterface::RecoEventData* recoEvt = new RecoInterface::RecoEventData();
    drInterface->read(evt);

    for (auto towerItr = evt.towers.begin(); towerItr != evt.towers.end(); ++towerItr) {
      auto tower = *towerItr;
      recoTower->reconstruct(tower,*recoEvt);

      auto theTower = recoTower->getTower();
      recoEvt->E_C += theTower.E_C;
      recoEvt->E_S += theTower.E_S;
      recoEvt->E_Scorr += theTower.E_Scorr;
      //recoEvt->n_C += theTower.n_C;
      //recoEvt->n_S += theTower.n_S;
    } // tower loop
    recoEvt->E_DR = RecoTower::E_DR(recoEvt->E_C,recoEvt->E_S);
    recoEvt->E_DRcorr = RecoTower::E_DR(recoEvt->E_C,recoEvt->E_Scorr);

    recoInterface->fill(recoEvt);
    delete recoEvt;

  } // event loop

  drInterface->close();
  recoInterface->write();
  recoInterface->close();

  return 0;
}
