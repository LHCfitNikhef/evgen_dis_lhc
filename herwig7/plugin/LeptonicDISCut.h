// -*- C++ -*-
#ifndef FASER_LeptonicDISCut_H
#define FASER_LeptonicDISCut_H
//
// LeptonicDISCut: a DIS Q2/y generation cut built ONLY from the two lepton
// momenta, for the FASER evgen benchmark.
//
// WHY THIS EXISTS.  ThePEG's SimpleDISCut gets Q2 right -- it builds it from
// the incoming and outgoing lepton momenta -- but it builds y as
//
//     y = Q2 / (SMax * x)      with x the PARTON momentum fraction
//
// (SimpleDISCut.cc, "double y = Q2/cuts->SMax()/x;").  At LO x = x_Bjorken and
// that is fine.  At NLO it is not: the Catani-Seymour dipole mapping rescales
// the incoming parton while leaving both lepton momenta untouched, so a y cut
// accepts a real-emission point and rejects its subtraction counterterm (or
// vice versa) and the IR cancellation is destroyed.  Measured symptom in this
// benchmark: with Miny/Maxy = 0.2/0.9 the NLO run returned sigma = -5(2) nb
// with 52% negative weights.
//
// The workaround was to drop the y cut at generation entirely, which cured the
// negative cross-section but opened the y -> 1 corner: 64% of generated events
// then fall outside the fiducial region, and the |w| > 5 rate is 113 per 10k
// events for y > 0.95 against 9.7 in the fiducial region, leaving N_eff = 1%.
//
// THE FIX.  y is an invariant of the LEPTON momenta alone,
//
//     y = (P.q)/(P.k),   q = k - k'
//
// and the dipole mapping does not touch k or k'.  So a y cut built this way is
// the same for a real-emission point and its counterterm, and the cancellation
// survives.  P is the incoming hadron momentum, which this class does not see
// -- but in the collision c.m. frame the beams are back-to-back along z, so
// with the lepton along +z
//
//     P ~ E_B (1, 0, 0, -1)  =>  y = (q.e() + q.z()) / (2 E_A)
//
// and E_A is just the incoming lepton energy, which IS available.  The proton
// mass is neglected here; that is an O(m_p^2/s) = 0.05% effect on a generation
// cut whose only job is to keep the pathological corner out of the sample.
// The ANALYSIS still applies the exact fiducial cut.
//
// Generate slightly WIDER than the fiducial region (e.g. 0.15 < y < 0.93 for a
// 0.2 < y < 0.9 analysis) so that migration across the boundary is simulated in
// both directions and sigma_fid = sigma_gen x (fiducial weight fraction) stays
// unbiased.
//
// MinW2 (added 2026-09-11) IS THE SAME CUT IN THE OTHER VARIABLE, for the
// study of replacing the y window by a W cut.  W2 = M_h^2 + y*s - Q2 and the
// hadron mass is dropped here, as it is in y above, so the cut this class
// applies is TIGHTER than the true W2 by M_h^2 ~ 0.88 GeV2 -- set it BELOW the
// analysis cut (8 GeV2 for a W > 3 GeV analysis) so the generation region
// stays the wider of the two.  s is taken as the cuts' SMax, which for fixed
// beams is the hadronic (k+P)^2 and is exactly what y multiplies.
//
// A W CUT DOES NOT REPLACE Maxy.  y -> 1 is where Herwig's NLO weights blow up
// (measured: |w| > 5 at 113 per 10k events above y = 0.95, against 9.7 in the
// fiducial region, N_eff ~ 1%) and W is LARGEST there, so nothing in a W cut
// keeps that corner out.  Set both when replacing one by the other.
//

#include "ThePEG/Cuts/TwoCutBase.h"

namespace ThePEG {

class LeptonicDISCut: public TwoCutBase {

public:

  LeptonicDISCut()
    : theMinQ2(4.0*GeV2), theMaxQ2(1000000.0*GeV2),
      theMiny(0.0), theMaxy(1.0), theMinW2(ZERO), chargedCurrent(false) {}

public:

  /** Minimum invariant mass squared -- no constraint from this cut. */
  virtual Energy2 minSij(tcPDPtr pi, tcPDPtr pj) const;
  virtual Energy2 minTij(tcPDPtr pi, tcPDPtr po) const;
  virtual double minDeltaR(tcPDPtr pi, tcPDPtr pj) const;
  virtual Energy minKTClus(tcPDPtr pi, tcPDPtr pj) const;
  virtual double minDurham(tcPDPtr pi, tcPDPtr pj) const;

  /** The actual cut. */
  virtual bool passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
                        LorentzMomentum pi, LorentzMomentum pj,
                        bool inci = false, bool incj = false) const;

public:

  void persistentOutput(PersistentOStream & os) const;
  void persistentInput(PersistentIStream & is, int version);
  static void Init();

protected:

  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;

private:

  /** Is this leg the lepton line we cut on? */
  bool check(long idi, long idj) const;

  Energy2 theMinQ2;
  Energy2 theMaxQ2;
  double theMiny;
  double theMaxy;
  Energy2 theMinW2;
  bool chargedCurrent;

  static ClassDescription<LeptonicDISCut> initLeptonicDISCut;
  LeptonicDISCut & operator=(const LeptonicDISCut &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

template <>
struct BaseClassTrait<LeptonicDISCut,1> {
  typedef TwoCutBase NthBase;
};

template <>
struct ClassTraits<LeptonicDISCut>
  : public ClassTraitsBase<LeptonicDISCut> {
  static string className() { return "ThePEG::LeptonicDISCut"; }
  static string library() { return "LeptonicDISCut.so"; }
};

}

#endif
