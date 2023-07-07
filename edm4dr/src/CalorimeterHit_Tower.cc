// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "edm4hep/CalorimeterHit_Tower.h"
#include "edm4hep/CalorimeterHit_TowerConst.h"
#include "edm4hep/CalorimeterHit_TowerObj.h"
#include "edm4hep/CalorimeterHit_TowerData.h"
#include "edm4hep/CalorimeterHit_TowerCollection.h"


#include <ostream>

namespace edm4hep {


CalorimeterHit_Tower::CalorimeterHit_Tower() : m_obj(new CalorimeterHit_TowerObj()) {
  m_obj->acquire();
}

CalorimeterHit_Tower::CalorimeterHit_Tower(float energy, float energyError, float time, short type, unsigned long long cellID, int iy, int ix, int numy, int numx, int noPhi, int noTheta) : m_obj(new CalorimeterHit_TowerObj()) {
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

CalorimeterHit_Tower::CalorimeterHit_Tower(const CalorimeterHit_Tower& other) : m_obj(other.m_obj) {
  m_obj->acquire();
}

CalorimeterHit_Tower& CalorimeterHit_Tower::operator=(const CalorimeterHit_Tower& other) {
  if (m_obj) m_obj->release();
  m_obj = other.m_obj;
  return *this;
}

CalorimeterHit_Tower::CalorimeterHit_Tower( CalorimeterHit_TowerObj* obj) : m_obj(obj) {
  if (m_obj) m_obj->acquire();
}

CalorimeterHit_Tower CalorimeterHit_Tower::clone() const {
  return {new CalorimeterHit_TowerObj(*m_obj)};
}

CalorimeterHit_Tower::~CalorimeterHit_Tower() {
  if (m_obj) m_obj->release();
}
CalorimeterHit_Tower::operator ConstCalorimeterHit_Tower() const { return ConstCalorimeterHit_Tower(m_obj); }

const float& CalorimeterHit_Tower::getEnergy() const { return m_obj->data.energy; }
const float& CalorimeterHit_Tower::getEnergyError() const { return m_obj->data.energyError; }
const float& CalorimeterHit_Tower::getTime() const { return m_obj->data.time; }
const short& CalorimeterHit_Tower::getType() const { return m_obj->data.type; }
const unsigned long long& CalorimeterHit_Tower::getCellID() const { return m_obj->data.cellID; }
const int& CalorimeterHit_Tower::getIy() const { return m_obj->data.iy; }
const int& CalorimeterHit_Tower::getIx() const { return m_obj->data.ix; }
const int& CalorimeterHit_Tower::getNumy() const { return m_obj->data.numy; }
const int& CalorimeterHit_Tower::getNumx() const { return m_obj->data.numx; }
const int& CalorimeterHit_Tower::getNoPhi() const { return m_obj->data.noPhi; }
const int& CalorimeterHit_Tower::getNoTheta() const { return m_obj->data.noTheta; }


void CalorimeterHit_Tower::setEnergy(float value) { m_obj->data.energy = value; }
void CalorimeterHit_Tower::setEnergyError(float value) { m_obj->data.energyError = value; }
void CalorimeterHit_Tower::setTime(float value) { m_obj->data.time = value; }
void CalorimeterHit_Tower::setType(short value) { m_obj->data.type = value; }
void CalorimeterHit_Tower::setCellID(unsigned long long value) { m_obj->data.cellID = value; }
void CalorimeterHit_Tower::setIy(int value) { m_obj->data.iy = value; }
void CalorimeterHit_Tower::setIx(int value) { m_obj->data.ix = value; }
void CalorimeterHit_Tower::setNumy(int value) { m_obj->data.numy = value; }
void CalorimeterHit_Tower::setNumx(int value) { m_obj->data.numx = value; }
void CalorimeterHit_Tower::setNoPhi(int value) { m_obj->data.noPhi = value; }
void CalorimeterHit_Tower::setNoTheta(int value) { m_obj->data.noTheta = value; }







bool CalorimeterHit_Tower::isAvailable() const {
  if (m_obj) {
    return true;
  }
  return false;
}

const podio::ObjectID CalorimeterHit_Tower::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{podio::ObjectID::invalid, podio::ObjectID::invalid};
}

bool CalorimeterHit_Tower::operator==(const ConstCalorimeterHit_Tower& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const ConstCalorimeterHit_Tower& value) {
  o << " id: " << value.id() << '\n';
  o << " energy : " << value.getEnergy() << '\n';
  o << " energyError : " << value.getEnergyError() << '\n';
  o << " time : " << value.getTime() << '\n';
  o << " type : " << value.getType() << '\n';
  o << " cellID : " << value.getCellID() << '\n';
  o << " iy : " << value.getIy() << '\n';
  o << " ix : " << value.getIx() << '\n';
  o << " numy : " << value.getNumy() << '\n';
  o << " numx : " << value.getNumx() << '\n';
  o << " noPhi : " << value.getNoPhi() << '\n';
  o << " noTheta : " << value.getNoTheta() << '\n';



  return o;
}

} // namespace edm4hep

