// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_SparseVectorCollection_H
#define EDM4HEP_SparseVectorCollection_H

// datamodel specific includes
#include "edm4hep/SparseVectorData.h"
#include "edm4hep/SparseVector.h"
#include "edm4hep/SparseVectorObj.h"

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


using SparseVectorDataContainer = std::vector<SparseVectorData>;
using SparseVectorObjPointerContainer = std::deque<SparseVectorObj*>;

class SparseVectorCollectionIterator {
public:
  SparseVectorCollectionIterator(int index, const SparseVectorObjPointerContainer* collection) : m_index(index), m_object(nullptr), m_collection(collection) {}

  bool operator!=(const SparseVectorCollectionIterator& x) const {
    return m_index != x.m_index; // TODO: may not be complete
  }

  const SparseVector operator*() const;
  const SparseVector* operator->() const;
  const SparseVectorCollectionIterator& operator++() const;

private:
  mutable int m_index;
  mutable SparseVector m_object;
  const SparseVectorObjPointerContainer* m_collection;
};

/**
A Collection is identified by an ID.
*/
class SparseVectorCollection : public podio::CollectionBase {

public:
  using const_iterator = const SparseVectorCollectionIterator;

  SparseVectorCollection();
//  SparseVectorCollection(const SparseVectorCollection& ) = delete; // deletion doesn't work w/ ROOT IO ! :-(
//  SparseVectorCollection(SparseVectorVector* data, int collectionID);
  ~SparseVectorCollection();

  void clear() override final;

  /// operator to allow pointer like calling of members a la LCIO
  SparseVectorCollection* operator->() { return (SparseVectorCollection*) this; }

  /// Append a new object to the collection, and return this object.
  SparseVector create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given
  template<typename... Args>
  SparseVector create(Args&&... args);

  /// number of elements in the collection
  size_t size() const override final;

  /// fully qualified type name of elements - with namespace
  std::string getValueTypeName() const override { return std::string("edm4hep::SparseVector"); }

  /// Returns the const object of given index
  const SparseVector operator[](unsigned int index) const;
  /// Returns the object of a given index
  SparseVector operator[](unsigned int index);
  /// Returns the const object of given index
  const SparseVector at(unsigned int index) const;
  /// Returns the object of given index
  SparseVector at(unsigned int index);


  /// Append object to the collection
  void push_back(ConstSparseVector object);

  void prepareForWrite() override final;
  void prepareAfterRead() override final;
  void setBuffer(void* address) override final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) override final;

  podio::CollRefCollection* referenceCollections() override final { return &m_refCollections; }

  podio::VectorMembersInfo* vectorMembers() override { return &m_vecmem_info; }

  void setID(unsigned ID) override final {
    m_collectionID = ID;
    std::for_each(m_entries.begin(),m_entries.end(),
                  [ID] (SparseVectorObj* obj) { obj->id = {obj->id.index, static_cast<int>(ID)}; }
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
  std::vector<SparseVectorData>* _getBuffer() { return m_data; }

  template<size_t arraysize>
  const std::array<float, arraysize> sampling() const;
  template<size_t arraysize>
  const std::array<edm4hep::ObjectID, arraysize> assocObj() const;

private:
  bool m_isValid;
  bool m_isReadFromFile{false};
  int m_collectionID;
  SparseVectorObjPointerContainer m_entries;

  // members to handle 1-to-N-relations

  // members to handle vector members
  std::vector<float>* m_vec_centers; /// combined vector of all objects in collection
  std::vector<std::vector<float>*> m_vecs_centers; /// pointers to individual member vectors
  std::vector<float>* m_vec_contents; /// combined vector of all objects in collection
  std::vector<std::vector<float>*> m_vecs_contents; /// pointers to individual member vectors
  // members to handle streaming
  podio::CollRefCollection m_refCollections;
  podio::VectorMembersInfo m_vecmem_info;
  SparseVectorDataContainer* m_data;
};

std::ostream& operator<<(std::ostream& o, const SparseVectorCollection& v);

template<typename... Args>
SparseVector SparseVectorCollection::create(Args&&... args) {
  const int size = m_entries.size();
  auto obj = new SparseVectorObj({size, m_collectionID}, {args...});
  m_entries.push_back(obj);
  return SparseVector(obj);
}

template<size_t arraysize>
const std::array<float, arraysize> SparseVectorCollection::sampling() const {
  std::array<float, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.sampling;
  }
  return tmp;
}

template<size_t arraysize>
const std::array<edm4hep::ObjectID, arraysize> SparseVectorCollection::assocObj() const {
  std::array<edm4hep::ObjectID, arraysize> tmp;
  const auto valid_size = std::min(arraysize, m_entries.size());
  for (unsigned i = 0; i < valid_size; ++i) {
    tmp[i] = m_entries[i]->data.assocObj;
  }
  return tmp;
}


} // namespace edm4hep


#endif
