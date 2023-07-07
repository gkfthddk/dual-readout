// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "edm4hep/SparseVector.h"
#include "edm4hep/SparseVectorConst.h"
#include "edm4hep/SparseVectorObj.h"
#include "edm4hep/SparseVectorData.h"
#include "edm4hep/SparseVectorCollection.h"


#include <ostream>

namespace edm4hep {


SparseVector::SparseVector() : m_obj(new SparseVectorObj()) {
  m_obj->acquire();
}

SparseVector::SparseVector(float sampling, edm4hep::ObjectID assocObj) : m_obj(new SparseVectorObj()) {
  m_obj->acquire();
  m_obj->data.sampling = sampling;
  m_obj->data.assocObj = assocObj;
}

SparseVector::SparseVector(const SparseVector& other) : m_obj(other.m_obj) {
  m_obj->acquire();
}

SparseVector& SparseVector::operator=(const SparseVector& other) {
  if (m_obj) m_obj->release();
  m_obj = other.m_obj;
  return *this;
}

SparseVector::SparseVector( SparseVectorObj* obj) : m_obj(obj) {
  if (m_obj) m_obj->acquire();
}

SparseVector SparseVector::clone() const {
  return {new SparseVectorObj(*m_obj)};
}

SparseVector::~SparseVector() {
  if (m_obj) m_obj->release();
}
SparseVector::operator ConstSparseVector() const { return ConstSparseVector(m_obj); }

const float& SparseVector::getSampling() const { return m_obj->data.sampling; }
const edm4hep::ObjectID& SparseVector::getAssocObj() const { return m_obj->data.assocObj; }


void SparseVector::setSampling(float value) { m_obj->data.sampling = value; }
void SparseVector::setAssocObj(edm4hep::ObjectID value) { m_obj->data.assocObj = value; }
edm4hep::ObjectID& SparseVector::assocObj() { return m_obj->data.assocObj; }


void SparseVector::addToCenters(float component) {
  m_obj->m_centers->push_back(component);
  m_obj->data.centers_end++;
}

std::vector<float>::const_iterator SparseVector::centers_begin() const {
  auto ret_value = m_obj->m_centers->begin();
  std::advance(ret_value, m_obj->data.centers_begin);
  return ret_value;
}

std::vector<float>::const_iterator SparseVector::centers_end() const {
  auto ret_value = m_obj->m_centers->begin();
  std::advance(ret_value, m_obj->data.centers_end);
  return ret_value;
}

unsigned int SparseVector::centers_size() const {
  return m_obj->data.centers_end - m_obj->data.centers_begin;
}

float SparseVector::getCenters(unsigned int index) const {
  if (centers_size() > index) {
    return m_obj->m_centers->at(m_obj->data.centers_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<float> SparseVector::getCenters() const {
  auto begin = m_obj->m_centers->begin();
  std::advance(begin, m_obj->data.centers_begin);
  auto end = m_obj->m_centers->begin();
  std::advance(end, m_obj->data.centers_end);
  return {begin, end};
}

void SparseVector::addToContents(float component) {
  m_obj->m_contents->push_back(component);
  m_obj->data.contents_end++;
}

std::vector<float>::const_iterator SparseVector::contents_begin() const {
  auto ret_value = m_obj->m_contents->begin();
  std::advance(ret_value, m_obj->data.contents_begin);
  return ret_value;
}

std::vector<float>::const_iterator SparseVector::contents_end() const {
  auto ret_value = m_obj->m_contents->begin();
  std::advance(ret_value, m_obj->data.contents_end);
  return ret_value;
}

unsigned int SparseVector::contents_size() const {
  return m_obj->data.contents_end - m_obj->data.contents_begin;
}

float SparseVector::getContents(unsigned int index) const {
  if (contents_size() > index) {
    return m_obj->m_contents->at(m_obj->data.contents_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<float> SparseVector::getContents() const {
  auto begin = m_obj->m_contents->begin();
  std::advance(begin, m_obj->data.contents_begin);
  auto end = m_obj->m_contents->begin();
  std::advance(end, m_obj->data.contents_end);
  return {begin, end};
}






bool SparseVector::isAvailable() const {
  if (m_obj) {
    return true;
  }
  return false;
}

const podio::ObjectID SparseVector::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{podio::ObjectID::invalid, podio::ObjectID::invalid};
}

bool SparseVector::operator==(const ConstSparseVector& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const ConstSparseVector& value) {
  o << " id: " << value.id() << '\n';
  o << " sampling : " << value.getSampling() << '\n';
  o << " assocObj : " << value.getAssocObj() << '\n';


  o << " centers : ";
  for (unsigned i = 0; i < value.centers_size(); ++i) {
    o << value.getCenters(i) << " ";
  }
  o << '\n';
  o << " contents : ";
  for (unsigned i = 0; i < value.contents_size(); ++i) {
    o << value.getContents(i) << " ";
  }
  o << '\n';

  return o;
}

} // namespace edm4hep

