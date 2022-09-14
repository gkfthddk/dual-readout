#include "P8ptcgun.h"

P8ptcgun::P8ptcgun(int id, double ee1, double ee2, double thetaIn, double phiIn, int seed, int mod)
: fId(id), fE1(ee1), fE2(ee2), fThetaIn(thetaIn), fPhiIn(phiIn), fSeed(seed), fmod(mod){
  fRandom = new TRandom(seed);
}

P8ptcgun::~P8ptcgun() {}

void P8ptcgun::fillResonance(Event& event, ParticleData& pdt, Rndm& /*rndm*/, bool atRest = false) {

  // Reset event record to allow for new event.
  event.reset();
  fE=fRandom->Uniform(fE1,fE2);
  // Select particle mass; where relevant according to Breit-Wigner.
  double mm = pdt.mSel(fId);
  double pp = sqrtpos(fE*fE - mm*mm);

  // Special case when particle is supposed to be at rest.
  if (atRest) {
    fE = mm;
    pp = 0.;
  }

  // Angles as input or uniform in solid angle.
  double cThe, sThe, phi;
  double ThetaIn, PhiIn;
  //1.51529  1.58191 //theta 2 -1
  //-0.0222021  0.0444041//phi 282 2
  if(fmod==0){
    ThetaIn=fThetaIn;
    PhiIn=fPhiIn;
  }
  if(fmod>0){
    ThetaIn=fRandom->Uniform(fThetaIn-0.02221*fmod,fThetaIn+0.02221*(fmod-1));
    PhiIn=fRandom->Uniform(fPhiIn-0.0222021*(fmod-1),fPhiIn+0.0222021*fmod);
  }
  cThe = cos(ThetaIn);
  sThe = sin(ThetaIn);
  phi  = PhiIn;

  // Store the particle in the event record.
  event.append( fId, 1, 0, 0, pp * sThe * cos(phi), pp * sThe * sin(phi), pp * cThe, fE, mm);
}

void P8ptcgun::fillParton(Event& event, ParticleData& pdt, Rndm& /*rndm*/, bool atRest = false, double scale = 20.) {

  // Reset event record to allow for new event.
  event.reset();
  fE=fRandom->Uniform(fE1,fE2);

  // Select particle mass; where relevant according to Breit-Wigner.
  double mm = pdt.mSel(fId);
  double pp = sqrtpos(fE*fE - mm*mm);

  // Special case when particle is supposed to be at rest.
  if (atRest) {
    fE = mm;
    pp = 0.;
  }

  // Angles as input or uniform in solid angle.
  double cThe, sThe, phi;
  cThe = cos(fThetaIn);
  sThe = sin(fThetaIn);
  phi  = fPhiIn;

  int col1, acol1, col2, acol2, aid;
  if (fId==21) { col1 = 101; acol1 = 102; col2 = 102; acol2 = 101; aid = fId; }
  else { col1 = 101; acol1 = 0; col2 = 0; acol2 = 101; aid = -fId; }

  // Store the particle in the event record.
  event.append( fId, 23, col1, acol1, pp * sThe * cos(phi), pp * sThe * sin(phi), pp * cThe, fE, mm, scale);
  event.append( aid, 23, col2, acol2, -pp * sThe * cos(phi), -pp * sThe * sin(phi), -pp * cThe, fE, mm, scale);
}
