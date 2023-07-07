// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_CalorimeterHit_TowerOBJ_H
#define EDM4HEP_CalorimeterHit_TowerOBJ_H

// data model specific includes
#include "edm4hep/CalorimeterHit_TowerData.h"

#include "podio/ObjBase.h"



namespace edm4hep {

class CalorimeterHit_Tower;
class ConstCalorimeterHit_Tower;

class CalorimeterHit_TowerObj : public podio::ObjBase {
public:
  /// constructor
  CalorimeterHit_TowerObj();
  /// copy constructor (does a deep-copy of relation containers)
  CalorimeterHit_TowerObj(const CalorimeterHit_TowerObj&);
  /// constructor from ObjectID and CalorimeterHit_TowerData
  /// does not initialize the internal relation containers
  CalorimeterHit_TowerObj(const podio::ObjectID id, CalorimeterHit_TowerData data);
  virtual ~CalorimeterHit_TowerObj();

public:
  CalorimeterHit_TowerData data;
};

} // namespace edm4hep


#endif
