// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_SparseVectorDATA_H
#define EDM4HEP_SparseVectorDATA_H

#include "edm4hep/ObjectID.h"

namespace edm4hep {


/** @class SparseVectorData
 *  User class to store sparse vector
 *  @author: Sanghyun Ko
 */
class SparseVectorData {
public:
  float sampling; ///< size of a bin
  ::edm4hep::ObjectID assocObj; ///< associated object ID

  unsigned int centers_begin;
  unsigned int centers_end;
  unsigned int contents_begin;
  unsigned int contents_end;
};

} // namespace edm4hep


#endif
