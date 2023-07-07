// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "edm4hep/CalorimeterHit_TowerCollection.h"


// standard includes
#include <stdexcept>
#include <iomanip>

namespace edm4hep {


CalorimeterHit_TowerCollection::CalorimeterHit_TowerCollection() :
  m_isValid(false), m_isReadFromFile(false), m_collectionID(0), m_entries(),
  m_data(new CalorimeterHit_TowerDataContainer()) {
}

CalorimeterHit_TowerCollection::~CalorimeterHit_TowerCollection() {
  clear();
  if (m_data) delete m_data;
}

const CalorimeterHit_Tower CalorimeterHit_TowerCollection::operator[](unsigned int index) const {
  return CalorimeterHit_Tower(m_entries[index]);
}

const CalorimeterHit_Tower CalorimeterHit_TowerCollection::at(unsigned int index) const {
  return CalorimeterHit_Tower(m_entries.at(index));
}

CalorimeterHit_Tower CalorimeterHit_TowerCollection::operator[](unsigned int index) {
  return CalorimeterHit_Tower(m_entries[index]);
}

CalorimeterHit_Tower CalorimeterHit_TowerCollection::at(unsigned int index) {
  return CalorimeterHit_Tower(m_entries.at(index));
}

size_t CalorimeterHit_TowerCollection::size() const {
  return m_entries.size();
}

CalorimeterHit_Tower CalorimeterHit_TowerCollection::create() {
  auto obj = m_entries.emplace_back(new CalorimeterHit_TowerObj());

  obj->id = {int(m_entries.size() - 1), m_collectionID};
  return CalorimeterHit_Tower(obj);
}

void CalorimeterHit_TowerCollection::clear() {
  m_data->clear();
  for (auto& obj : m_entries) { delete obj; }
  m_entries.clear();
}

void CalorimeterHit_TowerCollection::prepareForWrite() {
  const auto size = m_entries.size();
  m_data->reserve(size);
  for (auto& obj : m_entries) { m_data->push_back(obj->data); }

  // if the collection has been read from a file the rest of the information is
  // already in the correct format and we have to skip it, since the temporary
  // buffers are invalid
  if (m_isReadFromFile) return;
  for (auto& pointer : m_refCollections) { pointer->clear(); }


}

void CalorimeterHit_TowerCollection::prepareAfterRead() {
  int index = 0;
  for (auto& data : *m_data) {
    auto obj = new CalorimeterHit_TowerObj({index, m_collectionID}, data);

    m_entries.emplace_back(obj);
    ++index;
  }

  // at this point we are done with the I/O buffer and can safely clear it to not
  // have a redundant (but now useless) copy of the data
  m_data->clear();
  m_isValid = true;
  m_isReadFromFile = true;
}

bool CalorimeterHit_TowerCollection::setReferences(const podio::ICollectionProvider* collectionProvider) {

  return true; //TODO: check success
}

void CalorimeterHit_TowerCollection::push_back(ConstCalorimeterHit_Tower object) {
  const int size = m_entries.size();
  auto obj = object.m_obj;
  if (obj->id.index == podio::ObjectID::untracked) {
    obj->id = {size, m_collectionID};
    m_entries.push_back(obj);

  } else {
    throw std::invalid_argument("Object already in a collection. Cannot add it to a second collection");
  }
}

void CalorimeterHit_TowerCollection::setBuffer(void* address) {
  if (m_data) delete m_data;
  m_data = static_cast<CalorimeterHit_TowerDataContainer*>(address);
}

const CalorimeterHit_Tower CalorimeterHit_TowerCollectionIterator::operator*() const {
  m_object.m_obj = (*m_collection)[m_index];
  return m_object;
}

const CalorimeterHit_Tower* CalorimeterHit_TowerCollectionIterator::operator->() const {
  m_object.m_obj = (*m_collection)[m_index];
  return &m_object;
}

const CalorimeterHit_TowerCollectionIterator& CalorimeterHit_TowerCollectionIterator::operator++() const {
  ++m_index;
  return *this;
}

std::ostream& operator<<(std::ostream& o, const CalorimeterHit_TowerCollection& v) {
  const auto old_flags = o.flags();
  o << "          id:      energy: energyError:        time:        type:      cellID:          iy:          ix:        numy:        numx:       noPhi:     noTheta:" << '\n';

  for (size_t i = 0; i < v.size(); ++i) {
    o << std::scientific << std::showpos << std::setw(12) << v[i].id() << " "
      << std::setw(12) << v[i].getEnergy() << " "
      << std::setw(12) << v[i].getEnergyError() << " "
      << std::setw(12) << v[i].getTime() << " "
      << std::setw(12) << v[i].getType() << " "
      << std::setw(12) << v[i].getCellID() << " "
      << std::setw(12) << v[i].getIy() << " "
      << std::setw(12) << v[i].getIx() << " "
      << std::setw(12) << v[i].getNumy() << " "
      << std::setw(12) << v[i].getNumx() << " "
      << std::setw(12) << v[i].getNoPhi() << " "
      << std::setw(12) << v[i].getNoTheta() << " "
      << std::endl;




  }

  o.flags(old_flags);
  return o;
}

} // namespace edm4hep

