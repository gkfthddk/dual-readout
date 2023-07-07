#ifndef DRcalib2D_h
#define DRcalib2D_h 1

#include "k4FWCore/DataHandle.h"

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHit_TowerCollection.h"

#include "GridDRcalo.h"
#include "k4Interface/IGeoSvc.h"

class IGeoSvc;

class DRcalib2D : public GaudiAlgorithm {
public:
  DRcalib2D(const std::string& name, ISvcLocator* svcLoc);
  virtual ~DRcalib2D() {};

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

  void readCSV(std::string filename);

private:
  edm4hep::Vector3f getPosition(dd4hep::DDSegmentation::CellID& cID);

  ServiceHandle<IGeoSvc> m_geoSvc;
  dd4hep::DDSegmentation::GridDRcalo* pSeg;
  dd4hep::DDSegmentation::DRparamBase* pParamBase;

  DataHandle<edm4hep::RawCalorimeterHitCollection> m_digiHits{"DigiCalorimeterHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::CalorimeterHitCollection> m_caloHits{"DRcalo2dHits", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::CalorimeterHit_TowerCollection> m_caloTowers{"DRcalo2dTowers", Gaudi::DataHandle::Writer, this};

  Gaudi::Property<std::string> m_calibPath{this, "calibPath", "share/calib.csv", "relative path to calibration csv file"};
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "DRcaloSiPMreadout", "readout name of DRcalo"};

  Gaudi::Property<double> m_sampling{this, "sampling", 0.1, "SiPM sampling rate in ns"};

  std::vector<std::pair<float,float>> m_calibs;
};

#endif
