// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EDM4HEP_CalorimeterHit_TowerDATA_H
#define EDM4HEP_CalorimeterHit_TowerDATA_H


namespace edm4hep {


/** @class CalorimeterHit_TowerData
 *  User class to store calorimeterhit
 *  @author: Yunjae Lee
 */
class CalorimeterHit_TowerData {
public:
  float energy; ///< energy of the hit in [GeV].
  float energyError; ///< error of the hit energy in [GeV].
  float time; ///< time of the hit in [ns].
  short type; ///< type of hit. Mapping of integer types to names via collection parameters "CalorimeterHitTypeNames" and "CalorimeterHitTypeValues".
  unsigned long long cellID; ///< detector specific (geometrical) cell id.
  int iy; ///< readout index for theta-wise.
  int ix; ///< readout index for phi-wise.
  int numy; ///< number of readout for phi-wise.
  int numx; ///< number of readout theta-wise.
  int noPhi; ///< tower index for phi-wise.
  int noTheta; ///< tower index for theta-wise.

};

} // namespace edm4hep


#endif
