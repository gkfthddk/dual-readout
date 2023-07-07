// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_CalorimeterHit_Tower_H
#define EDM4HEP_CalorimeterHit_Tower_H

#include "edm4hep/CalorimeterHit_TowerConst.h"
#include "edm4hep/CalorimeterHit_TowerObj.h"

#include "podio/ObjectID.h"
#include <ostream>



namespace edm4hep {


class CalorimeterHit_TowerCollection;
class CalorimeterHit_TowerCollectionIterator;
class ConstCalorimeterHit_Tower;

/** @class CalorimeterHit_Tower
 *  User class to store calorimeterhit
 *  @author: Yunjae Lee
 */
class CalorimeterHit_Tower {

  friend CalorimeterHit_TowerCollection;
  friend CalorimeterHit_TowerCollectionIterator;
  friend ConstCalorimeterHit_Tower;

public:

  /// default constructor
  CalorimeterHit_Tower();
  CalorimeterHit_Tower(float energy, float energyError, float time, short type, unsigned long long cellID, int iy, int ix, int numy, int numx, int noPhi, int noTheta);

  /// constructor from existing CalorimeterHit_TowerObj
  CalorimeterHit_Tower(CalorimeterHit_TowerObj* obj);

  /// copy constructor
  CalorimeterHit_Tower(const CalorimeterHit_Tower& other);

  /// copy-assignment operator
  CalorimeterHit_Tower& operator=(const CalorimeterHit_Tower& other);

  /// support cloning (deep-copy)
  CalorimeterHit_Tower clone() const;

  /// destructor
  ~CalorimeterHit_Tower();

  /// conversion to const object
  operator ConstCalorimeterHit_Tower() const;

public:

  /// Access the energy of the hit in [GeV].
  const float& getEnergy() const;

  /// Access the error of the hit energy in [GeV].
  const float& getEnergyError() const;

  /// Access the time of the hit in [ns].
  const float& getTime() const;

  /// Access the type of hit. Mapping of integer types to names via collection parameters "CalorimeterHitTypeNames" and "CalorimeterHitTypeValues".
  const short& getType() const;

  /// Access the detector specific (geometrical) cell id.
  const unsigned long long& getCellID() const;

  /// Access the readout index for theta-wise.
  const int& getIy() const;

  /// Access the readout index for phi-wise.
  const int& getIx() const;

  /// Access the number of readout for phi-wise.
  const int& getNumy() const;

  /// Access the number of readout theta-wise.
  const int& getNumx() const;

  /// Access the tower index for phi-wise.
  const int& getNoPhi() const;

  /// Access the tower index for theta-wise.
  const int& getNoTheta() const;



  /// Set the energy of the hit in [GeV].
  void setEnergy(float value);

  /// Set the error of the hit energy in [GeV].
  void setEnergyError(float value);

  /// Set the time of the hit in [ns].
  void setTime(float value);

  /// Set the type of hit. Mapping of integer types to names via collection parameters "CalorimeterHitTypeNames" and "CalorimeterHitTypeValues".
  void setType(short value);

  /// Set the detector specific (geometrical) cell id.
  void setCellID(unsigned long long value);

  /// Set the readout index for theta-wise.
  void setIy(int value);

  /// Set the readout index for phi-wise.
  void setIx(int value);

  /// Set the number of readout for phi-wise.
  void setNumy(int value);

  /// Set the number of readout theta-wise.
  void setNumx(int value);

  /// Set the tower index for phi-wise.
  void setNoPhi(int value);

  /// Set the tower index for theta-wise.
  void setNoTheta(int value);






  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from CalorimeterHit_TowerObj instance
  void unlink() { m_obj = nullptr; }

  bool operator==(const CalorimeterHit_Tower& other) const { return m_obj == other.m_obj; }
  bool operator==(const ConstCalorimeterHit_Tower& other) const;

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const CalorimeterHit_Tower& other) const { return m_obj < other.m_obj; }

  unsigned int id() const { return getObjectID().collectionID * 10000000 + getObjectID().index; }

  const podio::ObjectID getObjectID() const;

private:
  CalorimeterHit_TowerObj* m_obj;
};

std::ostream& operator<<(std::ostream& o, const ConstCalorimeterHit_Tower& value);

} // namespace edm4hep


#endif
