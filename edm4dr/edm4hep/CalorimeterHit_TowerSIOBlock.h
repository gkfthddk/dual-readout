// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_CalorimeterHit_TowerSIOBlock_H
#define EDM4HEP_CalorimeterHit_TowerSIOBlock_H

#include "edm4hep/CalorimeterHit_TowerCollection.h"

#include "podio/SIOBlock.h"

#include <sio/api.h>
#include <sio/io_device.h>
#include <sio/version.h>

#include <typeindex>
#include <string>

namespace edm4hep {


class CalorimeterHit_TowerSIOBlock: public podio::SIOBlock {
public:
  CalorimeterHit_TowerSIOBlock() :
  SIOBlock("CalorimeterHit_Tower", sio::version::encode_version(0, 1)) {
    podio::SIOBlockFactory::instance().registerBlockForCollection("edm4hep::CalorimeterHit_Tower", this);
  }

  CalorimeterHit_TowerSIOBlock(const std::string& name) :
  // SIOBlock(name + "__CalorimeterHit_Tower", sio::version::encode_version(0, 1)) {}
  SIOBlock(name, sio::version::encode_version(0, 1)) {}

  // Read the particle data from the device
  virtual void read(sio::read_device& device, sio::version_type version) override;

  // Write the particle data to the device
  virtual void write(sio::write_device& device) override;

  virtual void createCollection() override {
    setCollection(new CalorimeterHit_TowerCollection);
  }

  SIOBlock* create(const std::string& name) const override { return new CalorimeterHit_TowerSIOBlock(name); }
};

static CalorimeterHit_TowerSIOBlock _dummyCalorimeterHit_TowerSIOBlock;

} // namespace edm4hep


#endif
