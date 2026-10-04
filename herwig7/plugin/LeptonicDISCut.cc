// -*- C++ -*-
//
// LeptonicDISCut: see LeptonicDISCut.h for the physics motivation.
//
#include "LeptonicDISCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/PDT/StandardMatchers.h"

using namespace ThePEG;

IBPtr LeptonicDISCut::clone() const {
  return new_ptr(*this);
}

IBPtr LeptonicDISCut::fullclone() const {
  return new_ptr(*this);
}

bool LeptonicDISCut::check(long idi, long idj) const {
  // idi is the incoming lepton, idj the outgoing one.  Neutral current keeps
  // the same particle; charged current steps one unit in |id| (nu_l <-> l).
  if ( !LeptonMatcher::Check(idi) ) return false;
  if ( chargedCurrent ) {
    if ( abs(idi) % 2 == 0 ) return idj == idi - 1 || idj == idi + 1;
    return idj == idi - 1 || idj == idi + 1;
  }
  return idj == idi;
}

Energy2 LeptonicDISCut::minSij(tcPDPtr, tcPDPtr) const { return ZERO; }
Energy2 LeptonicDISCut::minTij(tcPDPtr, tcPDPtr) const { return ZERO; }
double  LeptonicDISCut::minDeltaR(tcPDPtr, tcPDPtr) const { return 0.0; }
Energy  LeptonicDISCut::minKTClus(tcPDPtr, tcPDPtr) const { return ZERO; }
double  LeptonicDISCut::minDurham(tcPDPtr, tcPDPtr) const { return 0.0; }

bool LeptonicDISCut::
passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
         LorentzMomentum pi, LorentzMomentum pj,
         bool inci, bool incj) const {

  // Only act on the (incoming lepton, outgoing lepton) pair; everything else
  // passes untouched, exactly as SimpleDISCut does.
  LorentzMomentum k, kp;
  if ( inci && !incj ) {
    if ( !check(pitype->id(), pjtype->id()) ) return true;
    k = pi; kp = pj;
  } else if ( incj && !inci ) {
    if ( !check(pjtype->id(), pitype->id()) ) return true;
    k = pj; kp = pi;
  } else {
    return true;
  }

  // Q2 = -(k - k')^2, from the lepton momenta only -- invariant under the
  // dipole mapping, which is why SimpleDISCut's Q2 cut was always safe.
  const LorentzMomentum q = k - kp;
  const Energy2 Q2 = -q.m2();
  if ( Q2 <= theMinQ2 || Q2 >= theMaxQ2 ) return false;

  // y = (P.q)/(P.k).  In the collision c.m. frame the hadron beam is
  // back-to-back with the lepton along z, so P ~ E_B (1,0,0,-sgn) and the
  // E_B cancels in the ratio:  y = (q.e() + sgn*q.z()) / (2 E_A).
  const Energy EA = k.e();
  if ( EA <= ZERO ) return false;
  const double sgn = ( k.z() > ZERO ) ? 1.0 : -1.0;
  const double y = ( q.e() + sgn*q.z() ) / ( 2.0*EA );
  if ( y <= theMiny || y >= theMaxy ) return false;

  // W2 = M_h^2 + y s - Q2, with the hadron mass dropped exactly as it is in y
  // above.  The omission makes this cut tighter than the true W2 by M_h^2, so
  // the card must set MinW2 BELOW the analysis cut (see the header).
  if ( theMinW2 > ZERO && parent ) {
    const Energy2 W2 = y*parent->SMax() - Q2;
    if ( W2 < theMinW2 ) return false;
  }

  return true;
}

void LeptonicDISCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMinQ2, GeV2) << ounit(theMaxQ2, GeV2)
     << theMiny << theMaxy << ounit(theMinW2, GeV2) << chargedCurrent;
}

void LeptonicDISCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMinQ2, GeV2) >> iunit(theMaxQ2, GeV2)
     >> theMiny >> theMaxy >> iunit(theMinW2, GeV2) >> chargedCurrent;
}

ClassDescription<LeptonicDISCut> LeptonicDISCut::initLeptonicDISCut;

void LeptonicDISCut::Init() {

  static ClassDocumentation<LeptonicDISCut> documentation
    ("A DIS Q2/y generation cut built only from the two lepton momenta, so "
     "that it is invariant under the dipole mapping and therefore safe to "
     "apply at NLO -- unlike SimpleDISCut, whose y uses the parton momentum "
     "fraction.");

  static Parameter<LeptonicDISCut,Energy2> interfaceMinQ2
    ("MinQ2", "The minimum \\f$Q^2\\f$.",
     &LeptonicDISCut::theMinQ2, GeV2, 4.0*GeV2, ZERO, ZERO,
     true, false, Interface::lowerlim);

  static Parameter<LeptonicDISCut,Energy2> interfaceMaxQ2
    ("MaxQ2", "The maximum \\f$Q^2\\f$.",
     &LeptonicDISCut::theMaxQ2, GeV2, 1000000.0*GeV2, ZERO, ZERO,
     true, false, Interface::lowerlim);

  static Parameter<LeptonicDISCut,double> interfaceMiny
    ("Miny", "The minimum inelasticity y, built from the lepton momenta.",
     &LeptonicDISCut::theMiny, 0.0, 0.0, 1.0,
     true, false, Interface::limited);

  static Parameter<LeptonicDISCut,double> interfaceMaxy
    ("Maxy", "The maximum inelasticity y, built from the lepton momenta.",
     &LeptonicDISCut::theMaxy, 1.0, 0.0, 1.0,
     true, false, Interface::limited);

  static Parameter<LeptonicDISCut,Energy2> interfaceMinW2
    ("MinW2", "The minimum hadronic invariant mass squared, W2 = y s - Q2 "
     "(the hadron mass is dropped, as in y).  Zero means no W cut.",
     &LeptonicDISCut::theMinW2, GeV2, ZERO, ZERO, ZERO,
     true, false, Interface::lowerlim);

  static Switch<LeptonicDISCut,bool> interfaceCC
    ("ChargedCurrent", "Charged or neutral current.",
     &LeptonicDISCut::chargedCurrent, false, true, false);
  static SwitchOption interfaceCCYes
    (interfaceCC, "Yes", "Charged current.", true);
  static SwitchOption interfaceCCNo
    (interfaceCC, "No", "Neutral current.", false);

}
