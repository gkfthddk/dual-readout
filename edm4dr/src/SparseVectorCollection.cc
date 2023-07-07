// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "edm4hep/SparseVectorCollection.h"

#include <numeric>

// standard includes
#include <stdexcept>
#include <iomanip>

namespace edm4hep {


SparseVectorCollection::SparseVectorCollection() :
  m_isValid(false), m_isReadFromFile(false), m_collectionID(0), m_entries(),
  m_data(new SparseVectorDataContainer()) {
  m_vecmem_info.emplace_back("float", &m_vec_centers);
  m_vec_centers = new std::vector<float>();
  m_vecmem_info.emplace_back("float", &m_vec_contents);
  m_vec_contents = new std::vector<float>();
}

SparseVectorCollection::~SparseVectorCollection() {
  clear();
  if (m_data) delete m_data;
  if (m_vec_centers) delete m_vec_centers;
  if (m_vec_contents) delete m_vec_contents;
}

const SparseVector SparseVectorCollection::operator[](unsigned int index) const {
  return SparseVector(m_entries[index]);
}

const SparseVector SparseVectorCollection::at(unsigned int index) const {
  return SparseVector(m_entries.at(index));
}

SparseVector SparseVectorCollection::operator[](unsigned int index) {
  return SparseVector(m_entries[index]);
}

SparseVector SparseVectorCollection::at(unsigned int index) {
  return SparseVector(m_entries.at(index));
}

size_t SparseVectorCollection::size() const {
  return m_entries.size();
}

SparseVector SparseVectorCollection::create() {
  auto obj = m_entries.emplace_back(new SparseVectorObj());
  m_vecs_centers.push_back(obj->m_centers);
  m_vecs_contents.push_back(obj->m_contents);

  obj->id = {int(m_entries.size() - 1), m_collectionID};
  return SparseVector(obj);
}

void SparseVectorCollection::clear() {
  m_data->clear();
  m_vec_centers->clear();
  for (auto& vec : m_vecs_centers) { delete vec; }
  m_vecs_centers.clear();

  m_vec_contents->clear();
  for (auto& vec : m_vecs_contents) { delete vec; }
  m_vecs_contents.clear();

  for (auto& obj : m_entries) { delete obj; }
  m_entries.clear();
}

void SparseVectorCollection::prepareForWrite() {
  const auto size = m_entries.size();
  m_data->reserve(size);
  for (auto& obj : m_entries) { m_data->push_back(obj->data); }

  // if the collection has been read from a file the rest of the information is
  // already in the correct format and we have to skip it, since the temporary
  // buffers are invalid
  if (m_isReadFromFile) return;
  for (auto& pointer : m_refCollections) { pointer->clear(); }

  const int centers_size = std::accumulate(m_entries.begin(), m_entries.end(), 0, [](int sum, const SparseVectorObj* obj) { return sum + obj->m_centers->size(); });
  m_vec_centers->reserve(centers_size);
  int centers_index = 0;
  const int contents_size = std::accumulate(m_entries.begin(), m_entries.end(), 0, [](int sum, const SparseVectorObj* obj) { return sum + obj->m_contents->size(); });
  m_vec_contents->reserve(contents_size);
  int contents_index = 0;
  for (int i = 0, size = m_data->size(); i != size; ++i) {
    (*m_data)[i].centers_begin = centers_index;
    (*m_data)[i].centers_end += centers_index;
    centers_index = (*m_data)[i].centers_end;
    for (const auto& it : (*m_vecs_centers[i])) { m_vec_centers->push_back(it); }

    (*m_data)[i].contents_begin = contents_index;
    (*m_data)[i].contents_end += contents_index;
    contents_index = (*m_data)[i].contents_end;
    for (const auto& it : (*m_vecs_contents[i])) { m_vec_contents->push_back(it); }

  }
}

void SparseVectorCollection::prepareAfterRead() {
  int index = 0;
  for (auto& data : *m_data) {
    auto obj = new SparseVectorObj({index, m_collectionID}, data);

    obj->m_centers = m_vec_centers;
    obj->m_contents = m_vec_contents;
    m_entries.emplace_back(obj);
    ++index;
  }

  // at this point we are done with the I/O buffer and can safely clear it to not
  // have a redundant (but now useless) copy of the data
  m_data->clear();
  m_isValid = true;
  m_isReadFromFile = true;
}

bool SparseVectorCollection::setReferences(const podio::ICollectionProvider* collectionProvider) {

  return true; //TODO: check success
}

void SparseVectorCollection::push_back(ConstSparseVector object) {
  const int size = m_entries.size();
  auto obj = object.m_obj;
  if (obj->id.index == podio::ObjectID::untracked) {
    obj->id = {size, m_collectionID};
    m_entries.push_back(obj);

    m_vecs_centers.push_back(obj->m_centers);
    m_vecs_contents.push_back(obj->m_contents);
  } else {
    throw std::invalid_argument("Object already in a collection. Cannot add it to a second collection");
  }
}

void SparseVectorCollection::setBuffer(void* address) {
  if (m_data) delete m_data;
  m_data = static_cast<SparseVectorDataContainer*>(address);
}

const SparseVector SparseVectorCollectionIterator::operator*() const {
  m_object.m_obj = (*m_collection)[m_index];
  return m_object;
}

const SparseVector* SparseVectorCollectionIterator::operator->() const {
  m_object.m_obj = (*m_collection)[m_index];
  return &m_object;
}

const SparseVectorCollectionIterator& SparseVectorCollectionIterator::operator++() const {
  ++m_index;
  return *this;
}

std::ostream& operator<<(std::ostream& o, const SparseVectorCollection& v) {
  const auto old_flags = o.flags();
  o << "          id:    sampling:assocObj [ index, collectionID]:" << '\n';

  for (size_t i = 0; i < v.size(); ++i) {
    o << std::scientific << std::showpos << std::setw(12) << v[i].id() << " "
      << std::setw(12) << v[i].getSampling() << " "
      << std::setw(12) << v[i].getAssocObj() << " "
      << std::endl;



    o << "      centers : ";
    for (unsigned j = 0, N = v[i].centers_size(); j < N; ++j) {
      o << v[i].getCenters(j) << " ";
    }
    o << std::endl;
    o << "      contents : ";
    for (unsigned j = 0, N = v[i].contents_size(); j < N; ++j) {
      o << v[i].getContents(j) << " ";
    }
    o << std::endl;

  }

  o.flags(old_flags);
  return o;
}

} // namespace edm4hep

