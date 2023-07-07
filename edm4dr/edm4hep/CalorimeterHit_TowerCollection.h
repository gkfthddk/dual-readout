// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_CalorimeterHit_TowerCollection_H
#define EDM4HEP_CalorimeterHit_TowerCollection_H

// datamodel specific includes
#include "edm4hep/CalorimeterHit_TowerData.h"
#include "edm4hep/CalorimeterHit_Tower.h"
#include "edm4hep/CalorimeterHit_TowerObj.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/CollectionIDTable.h"

#include <string>
#include <vector>
#include <deque>
#include <array>
#include <algorithm>
#include <ostream>

namespace edm4hep {


using CalorimeterHit_TowerDataContainer = std::vector<CalorimeterHit_TowerData>;
using CalorimeterHit_TowerObjPointerContainer = std::deque<CalorimeterHit_TowerObj*>;

class CalorimeterHit_TowerCollectionIterator {
public:
  CalorimeterHit_TowerCollectionIterator(int index, const CalorimeterHit_TowerObjPointerContainer* collection) : m_index(index), m_object(nullptr), m_collection(collection) {}

  bool operator!=(const CalorimeterHit_TowerCollectionIterator& x) const {
    return m_index != x.m_index; // TODO: may not be complete
  }

  const CalorimeterHit_Tower operator*() const;
  const CalorimeterHit_Tower* operator->() const;
  const CalorimeterHit_TowerCollectionIterator& operator++() const;

private:
  mutable int m_index;
  mutable CalorimeterHit_Tower m_object;
  const CalorimeterHit_TowerObjPointerContainer* m_collection;
};

/**
A Collection is identified by an ID.
*/
class CalorimeterHit_TowerCollection : public podio::CollectionBase {

public:
  using const_iterator = const CalorimeterHit_TowerCollectionIterator;

  CalorimeterHit_TowerCollection();
//  CalorimeterHit_TowerCollection(const CalorimeterHit_TowerCollection& ) = delete; // deletion doesn't work w/ ROOT IO ! :-(
//  CalorimeterHit_TowerCollection(CalorimeterHit_TowerVector* data, int collectionID);
  ~CalorimeterHit_TowerCollection();

  void clear() override final;

  /// operator to allow pointer like calling of members a la LCIO
  CalorimeterHit_TowerCollection* operator->() { return (CalorimeterHit_TowerCollection*) this; }

  /// Append a new object to the collection, and return this object.
  CalorimeterHit_Tower create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given
  template<typename... Args>
  CalorimeterHit_Tower create(Args&&... args);

  /// number of elements in the collection
  size_t size() const override final;

  /// fully qualified type name of elements - with namespace
  std::string getValueTypeName() const override { return std::string("edm4hep::CalorimeterHit_Tower"); }

  /// Returns the const object of given index
  const CalorimeterHit_Tower operator[](unsigned int index) const;
  /// Returns the object of a given index
  CalorimeterHit_Tower operator[](unsigned int index);
  /// Returns the const object of given index
  const CalorimeterHit_Tower at(unsigned int index) const;
  /// Returns the object of given index
  CalorimeterHit_Tower at(unsigned int index);


  /// Append object to the collection
  void push_back(ConstCalorimeterHit_Tower object);

  void prepareForWrite() override final;
  void prepareAfterRead() override final;
  void setBuffer(void* address) override final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) override final;

  podio::CollRefCollection* referenceCollections() override final { return &m_refCollections; }

  podio::VectorMembersInfo* vectorMembers() override { return &m_vecmem_info; }

  void setID(unsigned ID) override final {
    m_collectionID = ID;
    std::for_each(m_entries.begin(),m_entries.end(),
                  [ID] (CalorimeterHit_TowerObj* obj) { obj->id = {obj->id.index, static_cast<int>(ID)}; }
    );
  };

  unsigned getID() const override final {
    return m_collectionID;
  }

  bool isValid() const override final {
    return m_isValid;
  }

  // support for the iterator protocol
  const const_iterator begin() const {
    return const_iterator(0, &m_entries);
  }
  const const_iterator end() const {
    return const_iterator(m_entries.size(), &m_entries);
  }

  /// returns the address of the pointer to the data buffer
  void* getBufferAddress() override final { return (void*)&m_data; }

  /// Returns the pointer to the data buffer
  std::vector<CalorimeterHit_TowerData>* _getBuffer() { return m_data; }

  template<size_t arraysize>
  const std::array<float, arraysize> energy() const;
  template<size_t arraysize>
  const std::array<float, arraysize> energyError() const;
  template<size_t arraysize>
  const std::array<float, arraysize> time() const;
  template<size_t arraysize>
  const std::array<short, arraysize> type() const;
  template<size_t arraysize>
  const std::array<unsigned long long, arraysize> cellID() const;
  template<size_t arraysize>
  const std::array<int, arraysize> iy() const;
  template<size_t arraysize>
  const std::array<int, arraysize> ix() const;
  template<size_t arraysize>
  const std::array<int, arraysize> numy() const;
  template<size_t arraysize>
  const std::array<int, arraysize> numx() const;
  template<size_t arraysize>
  const std::array<int, arraysize> noPhi() const;
  template<size_t arraysize>
  const std::array<int, arraysize> noTheta() const;

private:
  bool m_isValid;
  bool m_isReadFromFile{false};
  int m_collectionID;
  CalorimeterHit_TowerObjPointerContainer m_entries;

  // members to handle 1-to-N-relations

  // members to handle vector members
  // members to handle streaming
  podio::CollRefCollection m_refCollections;
  podio::VectorMembersInfo m_vecmem_info;
  CalorimeterHit_TowerDataContainer* m_data;
};

std::ostream& operator<<(std::ostream& o, const CalorimeterHit_TowerCollection& v);

template<typename... Args>
CalorimeterHit_Tower CalorimeterHit_TowerCollection::create(Args&&... args) {
  const int size = m_entries.size();
  auto obj = new CalorimeterHit_TowerObj({size, m_collectionID}, {args...});
  m_entries.push_back(obj);
  return CalorimeterHit_Tower(obj);
}

template<size_t arraysize>
const std::array<float, arraysize> CalorimeterHit_TowerCollection::energy() const {
  std::array<float, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.energy;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<float, arraysize> CalorimeterHit_TowerCollection::energyError() const {
  std::array<float, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.energyError;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<float, arraysize> CalorimeterHit_TowerCollection::time() const {
  std::array<float, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.time;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<short, arraysize> CalorimeterHit_TowerCollection::type() const {
  std::array<short, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.type;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<unsigned long long, arraysize> CalorimeterHit_TowerCollection::cellID() const {
  std::array<unsigned long long, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.cellID;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<int, arraysize> CalorimeterHit_TowerCollection::iy() const {
  std::array<int, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.iy;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<int, arraysize> CalorimeterHit_TowerCollection::ix() const {
  std::array<int, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.ix;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<int, arraysize> CalorimeterHit_TowerCollection::numy() const {
  std::array<int, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.numy;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<int, arraysize> CalorimeterHit_TowerCollection::numx() const {
  std::array<int, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.numx;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<int, arraysize> CalorimeterHit_TowerCollection::noPhi() const {
  std::array<int, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.noPhi;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<int, arraysize> CalorimeterHit_TowerCollection::noTheta() const {
  std::array<int, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.noTheta;
  }
  return tmp;
}


} // namespace edm4hep


#endif
