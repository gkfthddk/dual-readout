// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "edm4hep/CalorimeterHit_TowerSIOBlock.h"

#include <sio/block.h>
#include <sio/io_device.h>
#include <sio/version.h>


namespace edm4hep {


void CalorimeterHit_TowerSIOBlock::read(sio::read_device& device, sio::version_type vers) {

  auto* dataVec = static_cast<edm4hep::CalorimeterHit_TowerCollection*>(_col)->_getBuffer();

  unsigned size(0);
  device.data( size );
  dataVec->resize(size);

  podio::handlePODDataSIO(device, &(*dataVec)[0], size);

  //---- read ref collections -----
  podio::CollRefCollection* refCols = _col->referenceCollections() ;
  for( auto* refC : *refCols ){
    device.data( size ) ;
    refC->resize(size) ;
    podio::handlePODDataSIO( device ,  &((*refC)[0]), size ) ;
  }

}

void CalorimeterHit_TowerSIOBlock::write(sio::write_device& device) {
  _col->prepareForWrite() ;

  auto * dataVec = static_cast<edm4hep::CalorimeterHit_TowerCollection*>(_col)->_getBuffer() ;

  unsigned size = dataVec->size() ;
  device.data( size ) ;

  podio::handlePODDataSIO( device ,  &(*dataVec)[0], size ) ;

  //---- write ref collections -----
  podio::CollRefCollection* refCols = _col->referenceCollections() ;
  for( auto* refC : *refCols ){
    size = refC->size() ;
    device.data( size ) ;
    podio::handlePODDataSIO( device ,  &((*refC)[0]), size ) ;
  }

}

} // namespace edm4hep

