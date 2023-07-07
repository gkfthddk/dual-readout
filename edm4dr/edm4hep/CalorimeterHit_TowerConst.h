// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_ConstCalorimeterHit_Tower_H
#define EDM4HEP_ConstCalorimeterHit_Tower_H

#include "edm4hep/CalorimeterHit_TowerObj.h"

#include "podio/ObjectID.h"



namespace edm4hep {


class CalorimeterHit_Tower;
class CalorimeterHit_TowerCollection;
class CalorimeterHit_TowerCollectionIterator;

/** @class ConstCalorimeterHit_Tower
 *  User class to store calorimeterhit
 *  @author: Yunjae Lee
 */
class ConstCalorimeterHit_Tower {

  friend CalorimeterHit_Tower;
  friend CalorimeterHit_TowerCollection;
  friend CalorimeterHit_TowerCollectionIterator;

public:
  /// default constructor
  ConstCalorimeterHit_Tower();
  ConstCalorimeterHit_Tower(float energy, float energyError, float time, short type, unsigned long long cellID, int iy, int ix, int numy, int numx, int noPhi, int noTheta);

  /// constructor from existing CalorimeterHit_TowerObj
  ConstCalorimeterHit_Tower(CalorimeterHit_TowerObj* obj);

  /// copy constructor
  ConstCalorimeterHit_Tower(const ConstCalorimeterHit_Tower& other);

  /// copy-assignment operator
  ConstCalorimeterHit_Tower& operator=(const ConstCalorimeterHit_Tower& other);

  /// support cloning (deep-copy)
  ConstCalorimeterHit_Tower clone() const;

  /// destructor
  ~ConstCalorimeterHit_Tower();


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





  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from CalorimeterHit_TowerObj instance
  void unlink() { m_obj = nullptr; }

  bool operator==(const ConstCalorimeterHit_Tower& other) const { return m_obj == other.m_obj; }
  bool operator==(const CalorimeterHit_Tower& other) const;

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const ConstCalorimeterHit_Tower& other) const { return m_obj < other.m_obj; }

  unsigned int id() const { return getObjectID().collectionID * 10000000 + getObjectID().index; }

  const podio::ObjectID getObjectID() const;

private:
  CalorimeterHit_TowerObj* m_obj;
};

} // namespace edm4hep


#endif
