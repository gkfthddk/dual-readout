// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_SparseVector_H
#define EDM4HEP_SparseVector_H

#include "edm4hep/SparseVectorConst.h"
#include "edm4hep/SparseVectorObj.h"

#include "edm4hep/ObjectID.h"
#include <vector>
#include "podio/ObjectID.h"
#include <ostream>



namespace edm4hep {


class SparseVectorCollection;
class SparseVectorCollectionIterator;
class ConstSparseVector;

/** @class SparseVector
 *  User class to store sparse vector
 *  @author: Sanghyun Ko
 */
class SparseVector {

  friend SparseVectorCollection;
  friend SparseVectorCollectionIterator;
  friend ConstSparseVector;

public:

  /// default constructor
  SparseVector();
  SparseVector(float sampling, edm4hep::ObjectID assocObj);

  /// constructor from existing SparseVectorObj
  SparseVector(SparseVectorObj* obj);

  /// copy constructor
  SparseVector(const SparseVector& other);

  /// copy-assignment operator
  SparseVector& operator=(const SparseVector& other);

  /// support cloning (deep-copy)
  SparseVector clone() const;

  /// destructor
  ~SparseVector();

  /// conversion to const object
  operator ConstSparseVector() const;

public:

  /// Access the size of a bin
  const float& getSampling() const;

  /// Access the associated object ID
  const edm4hep::ObjectID& getAssocObj() const;



  /// Set the size of a bin
  void setSampling(float value);

  /// Set the associated object ID
  void setAssocObj(edm4hep::ObjectID value);
  /// Get reference to associated object ID
  edm4hep::ObjectID& assocObj();



  void addToCenters(float);
  unsigned int centers_size() const;
  float getCenters(unsigned int) const;
  std::vector<float>::const_iterator centers_begin() const;
  std::vector<float>::const_iterator centers_end() const;
  podio::RelationRange<float> getCenters() const;
  void addToContents(float);
  unsigned int contents_size() const;
  float getContents(unsigned int) const;
  std::vector<float>::const_iterator contents_begin() const;
  std::vector<float>::const_iterator contents_end() const;
  podio::RelationRange<float> getContents() const;

 unsigned int centers_size_begin() const {return m_obj->data.centers_begin;}
 unsigned int centers_size_end() const {return m_obj->data.centers_end;}
 


  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from SparseVectorObj instance
  void unlink() { m_obj = nullptr; }

  bool operator==(const SparseVector& other) const { return m_obj == other.m_obj; }
  bool operator==(const ConstSparseVector& other) const;

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const SparseVector& other) const { return m_obj < other.m_obj; }

  unsigned int id() const { return getObjectID().collectionID * 10000000 + getObjectID().index; }

  const podio::ObjectID getObjectID() const;

private:
  SparseVectorObj* m_obj;
};

std::ostream& operator<<(std::ostream& o, const ConstSparseVector& value);

} // namespace edm4hep


#endif
