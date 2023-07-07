// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "edm4hep/SparseVector.h"
#include "edm4hep/SparseVectorConst.h"
#include "edm4hep/SparseVectorObj.h"
#include "edm4hep/SparseVectorData.h"
#include "edm4hep/SparseVectorCollection.h"


#include <ostream>

namespace edm4hep {


ConstSparseVector::ConstSparseVector() : m_obj(new SparseVectorObj()) {
  m_obj->acquire();
}

ConstSparseVector::ConstSparseVector(float sampling, edm4hep::ObjectID assocObj) : m_obj(new SparseVectorObj()) {
  m_obj->acquire();
  m_obj->data.sampling = sampling;
  m_obj->data.assocObj = assocObj;
}

ConstSparseVector::ConstSparseVector(const ConstSparseVector& other) : m_obj(other.m_obj) {
  m_obj->acquire();
}

ConstSparseVector& ConstSparseVector::operator=(const ConstSparseVector& other) {
  if (m_obj) m_obj->release();
  m_obj = other.m_obj;
  return *this;
}

ConstSparseVector::ConstSparseVector( SparseVectorObj* obj) : m_obj(obj) {
  if (m_obj) m_obj->acquire();
}

ConstSparseVector ConstSparseVector::clone() const {
  return {new SparseVectorObj(*m_obj)};
}

ConstSparseVector::~ConstSparseVector() {
  if (m_obj) m_obj->release();
}
const float& ConstSparseVector::getSampling() const { return m_obj->data.sampling; }
const edm4hep::ObjectID& ConstSparseVector::getAssocObj() const { return m_obj->data.assocObj; }



std::vector<float>::const_iterator ConstSparseVector::centers_begin() const {
  auto ret_value = m_obj->m_centers->begin();
  std::advance(ret_value, m_obj->data.centers_begin);
  return ret_value;
}

std::vector<float>::const_iterator ConstSparseVector::centers_end() const {
  auto ret_value = m_obj->m_centers->begin();
  std::advance(ret_value, m_obj->data.centers_end);
  return ret_value;
}

unsigned int ConstSparseVector::centers_size() const {
  return m_obj->data.centers_end - m_obj->data.centers_begin;
}

float ConstSparseVector::getCenters(unsigned int index) const {
  if (centers_size() > index) {
    return m_obj->m_centers->at(m_obj->data.centers_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<float> ConstSparseVector::getCenters() const {
  auto begin = m_obj->m_centers->begin();
  std::advance(begin, m_obj->data.centers_begin);
  auto end = m_obj->m_centers->begin();
  std::advance(end, m_obj->data.centers_end);
  return {begin, end};
}


std::vector<float>::const_iterator ConstSparseVector::contents_begin() const {
  auto ret_value = m_obj->m_contents->begin();
  std::advance(ret_value, m_obj->data.contents_begin);
  return ret_value;
}

std::vector<float>::const_iterator ConstSparseVector::contents_end() const {
  auto ret_value = m_obj->m_contents->begin();
  std::advance(ret_value, m_obj->data.contents_end);
  return ret_value;
}

unsigned int ConstSparseVector::contents_size() const {
  return m_obj->data.contents_end - m_obj->data.contents_begin;
}

float ConstSparseVector::getContents(unsigned int index) const {
  if (contents_size() > index) {
    return m_obj->m_contents->at(m_obj->data.contents_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<float> ConstSparseVector::getContents() const {
  auto begin = m_obj->m_contents->begin();
  std::advance(begin, m_obj->data.contents_begin);
  auto end = m_obj->m_contents->begin();
  std::advance(end, m_obj->data.contents_end);
  return {begin, end};
}





bool ConstSparseVector::isAvailable() const {
  if (m_obj) {
    return true;
  }
  return false;
}

const podio::ObjectID ConstSparseVector::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{podio::ObjectID::invalid, podio::ObjectID::invalid};
}

bool ConstSparseVector::operator==(const SparseVector& other) const {
  return m_obj == other.m_obj;
}

} // namespace edm4hep

