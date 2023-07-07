// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "edm4hep/SparseVectorObj.h"
namespace edm4hep {

SparseVectorObj::SparseVectorObj() :
  ObjBase{{podio::ObjectID::untracked, podio::ObjectID::untracked}, 0},
  data(),
  m_centers(new std::vector<float>()),
  m_contents(new std::vector<float>())
{  }

SparseVectorObj::SparseVectorObj(const podio::ObjectID id, SparseVectorData data) :
  ObjBase{id, 0}, data(data)
{  }

SparseVectorObj::SparseVectorObj(const SparseVectorObj& other) :
  ObjBase{{podio::ObjectID::untracked, podio::ObjectID::untracked}, 0},
  data(other.data),
  m_centers(new std::vector<float>(*(other.m_centers))),
  m_contents(new std::vector<float>(*(other.m_contents)))
{
}

SparseVectorObj::~SparseVectorObj() {
  if (id.index == podio::ObjectID::untracked) {
    delete m_centers;
    delete m_contents;
  }

}
} // namespace edm4hep

