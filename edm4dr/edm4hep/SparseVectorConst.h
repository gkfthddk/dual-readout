// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_ConstSparseVector_H
#define EDM4HEP_ConstSparseVector_H

#include "edm4hep/SparseVectorObj.h"

#include "edm4hep/ObjectID.h"
#include <vector>
#include "podio/ObjectID.h"



namespace edm4hep {


class SparseVector;
class SparseVectorCollection;
class SparseVectorCollectionIterator;

/** @class ConstSparseVector
 *  User class to store sparse vector
 *  @author: Sanghyun Ko
 */
class ConstSparseVector {

  friend SparseVector;
  friend SparseVectorCollection;
  friend SparseVectorCollectionIterator;

public:
  /// default constructor
  ConstSparseVector();
  ConstSparseVector(float sampling, edm4hep::ObjectID assocObj);

  /// constructor from existing SparseVectorObj
  ConstSparseVector(SparseVectorObj* obj);

  /// copy constructor
  ConstSparseVector(const ConstSparseVector& other);

  /// copy-assignment operator
  ConstSparseVector& operator=(const ConstSparseVector& other);

  /// support cloning (deep-copy)
  ConstSparseVector clone() const;

  /// destructor
  ~ConstSparseVector();


public:

  /// Access the size of a bin
  const float& getSampling() const;

  /// Access the associated object ID
  const edm4hep::ObjectID& getAssocObj() const;



  unsigned int centers_size() const;
  float getCenters(unsigned int) const;
  std::vector<float>::const_iterator centers_begin() const;
  std::vector<float>::const_iterator centers_end() const;
  podio::RelationRange<float> getCenters() const;
  unsigned int contents_size() const;
  float getContents(unsigned int) const;
  std::vector<float>::const_iterator contents_begin() const;
  std::vector<float>::const_iterator contents_end() const;
  podio::RelationRange<float> getContents() const;


  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from SparseVectorObj instance
  void unlink() { m_obj = nullptr; }

  bool operator==(const ConstSparseVector& other) const { return m_obj == other.m_obj; }
  bool operator==(const SparseVector& other) const;

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const ConstSparseVector& other) const { return m_obj < other.m_obj; }

  unsigned int id() const { return getObjectID().collectionID * 10000000 + getObjectID().index; }

  const podio::ObjectID getObjectID() const;

private:
  SparseVectorObj* m_obj;
};

} // namespace edm4hep


#endif
