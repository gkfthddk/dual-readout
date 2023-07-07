// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "edm4hep/CalorimeterHit_TowerObj.h"
namespace edm4hep {

CalorimeterHit_TowerObj::CalorimeterHit_TowerObj() :
  ObjBase{{podio::ObjectID::untracked, podio::ObjectID::untracked}, 0},
  data()
{  }

CalorimeterHit_TowerObj::CalorimeterHit_TowerObj(const podio::ObjectID id, CalorimeterHit_TowerData data) :
  ObjBase{id, 0}, data(data)
{  }

CalorimeterHit_TowerObj::CalorimeterHit_TowerObj(const CalorimeterHit_TowerObj& other) :
  ObjBase{{podio::ObjectID::untracked, podio::ObjectID::untracked}, 0},
  data(other.data)
{
}

CalorimeterHit_TowerObj::~CalorimeterHit_TowerObj() {

}
} // namespace edm4hep

