// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_SparseVectorOBJ_H
#define EDM4HEP_SparseVectorOBJ_H

// data model specific includes
#include "edm4hep/SparseVectorData.h"
#include "podio/RelationRange.h"
#include <vector>

#include "podio/ObjBase.h"
#include <vector>


namespace edm4hep {

class SparseVector;
class ConstSparseVector;

class SparseVectorObj : public podio::ObjBase {
public:
  /// constructor
  SparseVectorObj();
  /// copy constructor (does a deep-copy of relation containers)
  SparseVectorObj(const SparseVectorObj&);
  /// constructor from ObjectID and SparseVectorData
  /// does not initialize the internal relation containers
  SparseVectorObj(const podio::ObjectID id, SparseVectorData data);
  virtual ~SparseVectorObj();

public:
  SparseVectorData data;
  std::vector<float>* m_centers;
  std::vector<float>* m_contents;
};

} // namespace edm4hep


#endif
