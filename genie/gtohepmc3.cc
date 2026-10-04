// Convert GENIE GHEP ROOT output to HepMC3 Ascii for the muon-DIS benchmark.
//
// Writes one vertex per event: GHEP slots 0/1 (probe, target) become HepMC3
// status-4 beam records; kIStStableFinalState particles become status 1.
// That is all analyze.py reads. Units GeV/mm, weight left at the default 1
// (gevgen samples are unweighted).
//
// Usage: gtohepmc3 <in.ghep.root> <out.hepmc>

#include <iostream>
#include <memory>

#include <TFile.h>
#include <TTree.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"

// Charm test by PDG quark-content digits, IDENTICAL to analyze.py's
// is_charm_hadron().  This used to be a four-code whitelist
// (411/421/431/4122), which silently DROPPED every other charm hadron GENIE
// decays -- Xi_c and Omega_c decay weakly with ctau of tens of microns, so
// they are decayed records and never reached the HepMC file at all, while the
// analysis was counting charm by quark content and could not see them.  The
// two definitions must be the same one, so this is the same arithmetic:
// a hadron carries charm if any of its three quark digits is 4.
static bool IsCharmHadron(int pdg)
{
  const int a = abs(pdg);
  if (a < 100 || a >= 1000000000) return false;   // quarks/leptons/bosons; nuclei
  return (a / 10) % 10 == 4 || (a / 100) % 10 == 4 || (a / 1000) % 10 == 4;
}

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "usage: gtohepmc3 <in.ghep.root> <out.hepmc>\n";
    return 1;
  }

  TFile fin(argv[1]);
  if (fin.IsZombie()) { std::cerr << "cannot open " << argv[1] << "\n"; return 1; }
  TTree* tree = dynamic_cast<TTree*>(fin.Get("gtree"));
  if (!tree) { std::cerr << "no gtree in " << argv[1] << "\n"; return 1; }

  genie::NtpMCEventRecord* mcrec = nullptr;
  tree->SetBranchAddress("gmcrec", &mcrec);

  HepMC3::WriterAscii writer(argv[2]);
  const long n = tree->GetEntries();
  long nout = 0;

  for (long i = 0; i < n; i++) {
    tree->GetEntry(i);
    genie::EventRecord& ev = *(mcrec->event);

    HepMC3::GenEvent hev(HepMC3::Units::GEV, HepMC3::Units::MM);
    hev.set_event_number(i);
    auto vtx = std::make_shared<HepMC3::GenVertex>();

    // ON A NUCLEAR TARGET THE BEAM RECORD IS THE STRUCK NUCLEON (user,
    // 2026-10-01, the neutrino-generator appendix).  Slot 1 is then the
    // nucleus, and the analysis builds x, W and Q2 from the status-4 records:
    // with a 171 GeV tungsten "beam" x would be ~184 times too small and W
    // meaningless.  The struck nucleon, with its Fermi momentum, is what the
    // NuWro and GiBUU converters write as well.  On a free-nucleon run slot 1
    // IS the nucleon and the output is unchanged.
    const genie::GHepParticle* tgt1 = ev.Particle(1);
    const bool nuclear = tgt1 && tgt1->Pdg() > 1000000000;
    const genie::GHepParticle* hitnuc = nuclear ? ev.HitNucleon() : nullptr;
    if (nuclear && !hitnuc) {
      std::cerr << "event " << i << ": nuclear target without a hit nucleon"
                << " -- refusing\n";
      return 1;
    }

    for (int ip = 0; ip < ev.GetEntries(); ip++) {
      genie::GHepParticle* p = ev.Particle(ip);
      // charm hadrons always decay (ctau << 10 mm) but the analysis counts
      // them at production from decayed (status-2) records
      const bool is_charm_hadron = IsCharmHadron(p->Pdg());
      // decayed charm shows up as kIStDecayedState (GENIE decayer) or as
      // kIStDISPreFragmHadronicState (decayed internally by the Pythia8
      // hadronizers, which park ks<0 records there)
      const bool is_decayed_record =
          p->Status() == genie::kIStDecayedState ||
          p->Status() == genie::kIStDISPreFragmHadronicState;
      int hstat;
      if (nuclear && ip == 1)                              continue;
      if (ip == 0 || ip == 1 || p == hitnuc)               hstat = 4;
      else if (p->Status() == genie::kIStStableFinalState) hstat = 1;
      else if (is_charm_hadron && is_decayed_record)       hstat = 2;
      else continue;

      auto hp = std::make_shared<HepMC3::GenParticle>(
          HepMC3::FourVector(p->Px(), p->Py(), p->Pz(), p->E()),
          p->Pdg(), hstat);
      if (hstat == 4) vtx->add_particle_in(hp);
      else            vtx->add_particle_out(hp);
    }

    // hard-process channel tag from the GENIE interaction summary: struck
    // quark as a Pythia-convention status-21 record, final quark (set by
    // HEDIS channels, e.g. s -> c) as status 23. Momenta are not stored in
    // the summary, so the records are written with p = 0; the analysis only
    // reads their PDG codes.
    const genie::Interaction* summ = ev.Summary();
    if (summ) {
      const int hitq = summ->InitState().Tgt().HitQrkPdg();
      const int finq = summ->ExclTag().FinalQuarkPdg();
      if (hitq) vtx->add_particle_out(std::make_shared<HepMC3::GenParticle>(
          HepMC3::FourVector(0, 0, 0, 0), hitq, 21));
      if (finq) vtx->add_particle_out(std::make_shared<HepMC3::GenParticle>(
          HepMC3::FourVector(0, 0, 0, 0), finq, 23));
    }

    hev.add_vertex(vtx);
    writer.write_event(hev);
    mcrec->Clear();
    nout++;
  }

  writer.close();
  std::cout << "converted " << nout << " / " << n << " events -> " << argv[2] << "\n";
  return 0;
}
