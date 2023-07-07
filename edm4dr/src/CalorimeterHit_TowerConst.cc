// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "edm4hep/CalorimeterHit_Tower.h"
#include "edm4hep/CalorimeterHit_TowerConst.h"
#include "edm4hep/CalorimeterHit_TowerObj.h"
#include "edm4hep/CalorimeterHit_TowerData.h"
#include "edm4hep/CalorimeterHit_TowerCollection.h"


#include <ostream>

namespace edm4hep {


ConstCalorimeterHit_Tower::ConstCalorimeterHit_Tower() : m_obj(new CalorimeterHit_TowerObj()) {
  m_obj->acquire();
}

ConstCalorimeterHit_Tower::ConstCalorimeterHit_Tower(float energy, float energyError, float time, short type, unsigned long long cellID, int iy, int ix, int numy, int numx, int noPhi, int noTheta) : m_obj(new CalorimeterHit_TowerObj()) {
  m_obj->acquire();
  m_obj->data.energy = energy;
  m_obj->data.energyError = energyError;
  m_obj->data.time = time;
  m_obj->data.type = type;
  m_obj->data.cellID = cellID;
  m_obj->data.iy = iy;
  m_obj->data.ix = ix;
  m_obj->data.numy = numy;
  m_obj->data.numx = numx;
  m_obj->data.noPhi = noPhi;
  m_obj->data.noTheta = noTheta;
}

ConstCalorimeterHit_Tower::ConstCalorimeterHit_Tower(const ConstCalorimeterHit_Tower& other) : m_obj(other.m_obj) {
  m_obj->acquire();
}

ConstCalorimeterHit_Tower& ConstCalorimeterHit_Tower::operator=(const ConstCalorimeterHit_Tower& other) {
  if (m_obj) m_obj->release();
  m_obj = other.m_obj;
  return *this;
}

ConstCalorimeterHit_Tower::ConstCalorimeterHit_Tower( CalorimeterHit_TowerObj* obj) : m_obj(obj) {
  if (m_obj) m_obj->acquire();
}

ConstCalorimeterHit_Tower ConstCalorimeterHit_Tower::clone() const {
  return {new CalorimeterHit_TowerObj(*m_obj)};
}

ConstCalorimeterHit_Tower::~ConstCalorimeterHit_Tower() {
  if (m_obj) m_obj->release();
}
const float& ConstCalorimeterHit_Tower::getEnergy() const { return m_obj->data.energy; }
const float& ConstCalorimeterHit_Tower::getEnergyError() const { return m_obj->data.energyError; }
const float& ConstCalorimeterHit_Tower::getTime() const { return m_obj->data.time; }
const short& ConstCalorimeterHit_Tower::getType() const { return m_obj->data.type; }
const unsigned long long& ConstCalorimeterHit_Tower::getCellID() const { return m_obj->data.cellID; }
const int& ConstCalorimeterHit_Tower::getIy() const { return m_obj->data.iy; }
const int& ConstCalorimeterHit_Tower::getIx() const { return m_obj->data.ix; }
const int& ConstCalorimeterHit_Tower::getNumy() const { return m_obj->data.numy; }
const int& ConstCalorimeterHit_Tower::getNumx() const { return m_obj->data.numx; }
const int& ConstCalorimeterHit_Tower::getNoPhi() const { return m_obj->data.noPhi; }
const int& ConstCalorimeterHit_Tower::getNoTheta() const { return m_obj->data.noTheta; }






bool ConstCalorimeterHit_Tower::isAvailable() const {
  if (m_obj) {
    return true;
  }
  return false;
}

const podio::ObjectID ConstCalorimeterHit_Tower::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{podio::ObjectID::invalid, podio::ObjectID::invalid};
}

bool ConstCalorimeterHit_Tower::operator==(const CalorimeterHit_Tower& other) const {
  return m_obj == other.m_obj;
}

} // namespace edm4hep

