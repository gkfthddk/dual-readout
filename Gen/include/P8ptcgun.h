#ifndef P8ptcgun_h
#define P8ptcgun_h 1

#include "Pythia8/Pythia.h"
#include "TRandom.h"

using namespace Pythia8;

class P8ptcgun {
public:
  P8ptcgun(int id, double ee1, double ee2, double thetaIn, double phiIn, int seed);
  ~P8ptcgun();

  void fillResonance(Event& event, ParticleData& pdt, Rndm& rndm, bool atRest);
  void fillParton(Event& event, ParticleData& pdt, Rndm& rndm, bool atRest, double scale);

private:
  int fId;
  double fE;
  double fE1;
  double fE2;
  double fThetaIn;
  double fPhiIn;
  int fSeed;

  TRandom* fRandom;
};

#endif
