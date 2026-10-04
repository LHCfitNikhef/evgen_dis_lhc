// Pythia8 shower driver for POWHEG-DIS LHE events (muon-DIS benchmark):
// read a .cmnd card + a POWHEG .lhe file, shower with the DIS-tailored
// POWHEG veto hooks, write HepMC3 events plus a JSON sidecar.
// Usage: ./main_powheg card.cmnd events.lhe outprefix [maxEvents]
//
// The writer is closed explicitly and the process exits via _Exit: with the
// conda-built HepMC3/LHAPDF next to Apple clang, static destructors abort at
// teardown (exit 134) even though all output is complete (see main_dis.cc).

#include "Pythia8/Pythia.h"
#include "DISPowhegHooks.h"          // defines Pythia8::PowhegHooks (DIS maps)
#include "Pythia8Plugins/HepMC3.h"
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

using namespace Pythia8;

int main(int argc, char* argv[]) {

  if (argc < 4) {
    cerr << "usage: " << argv[0]
         << " card.cmnd events.lhe outprefix [maxEvents]" << endl;
    return 1;
  }
  std::string card = argv[1], lhef = argv[2], out = argv[3];
  long maxEvents = (argc > 4) ? atol(argv[4]) : 0;   // 0 = all

  Pythia pythia;

  // register the POWHEG:* modes before reading the card (as in main31.cc);
  // POWHEG:dis_map selects the DIS kinematic map in DISPowhegHooks
  pythia.settings.addMode("POWHEG:nFinal",    2, true, false, 1, 10);
  pythia.settings.addMode("POWHEG:veto",      0, true, true,  0, 2);
  pythia.settings.addMode("POWHEG:vetoCount", 10, true, false, 0, 0);
  pythia.settings.addMode("POWHEG:pThard",    0, true, true,  0, 2);
  pythia.settings.addMode("POWHEG:pTemt",     0, true, true,  0, 2);
  pythia.settings.addMode("POWHEG:emitted",   0, true, true,  0, 3);
  pythia.settings.addMode("POWHEG:pTdef",     0, true, true,  0, 2);
  pythia.settings.addMode("POWHEG:MPIveto",   0, true, true,  0, 1);
  pythia.settings.addMode("POWHEG:QEDveto",   0, true, true,  0, 2);
  pythia.settings.addMode("POWHEG:dis_map",   1, true, true,  0, 2);

  pythia.readFile(card);
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + lhef);

  // load the DIS POWHEG veto hooks
  shared_ptr<PowhegHooks> powhegHooks;
  if (pythia.settings.mode("POWHEG:veto") > 0) {
    powhegHooks = make_shared<PowhegHooks>();
    pythia.setUserHooksPtr(powhegHooks);
  }

  if (!pythia.init()) {
    cerr << "Pythia init failed" << endl;
    return 1;
  }

  HepMC3::Pythia8ToHepMC3 converter;
  HepMC3::WriterAscii writer(out + ".hepmc");

  // LHE input: a failed event is consumed from the file and cannot be
  // regenerated; count and continue (benchmark abort budget ~1%)
  long nDone = 0, nAborted = 0;
  double sumW = 0.0;
  bool runInfoNamed = false;             // names captured from the first event
  std::vector<std::string> weightNames;  // -> the JSON sidecar, see below
  while (true) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) break;
      ++nAborted;
      continue;
    }
    HepMC3::GenEvent evt;
    converter.fill_next_event(pythia, &evt);
    // THE LHE EVENT INDEX, so a showered event can be joined back to the
    // <rwgt> block it came from.
    //
    // WHY NOT JUST PROPAGATE PYTHIA'S WEIGHTS, which is done below as well:
    // because Pythia 8.311's LHEF weight handling is not trustworthy for this
    // file.  Measured on powheg-cmass rwgt-nu1TeV/pwgevents-rwgt.lhe, whose
    // header declares 1001-1007 and 2001-2004:
    //   * it exposes only THREE of the four PDF weights -- id 2004 is dropped
    //     silently, and 2004 is ATLASpdf21, the set furthest from the nominal;
    //   * its NAMES are its own canonical (muR, muF) enumeration and do NOT
    //     follow the file.  The values come back in LHE order, so names and
    //     values are MISMATCHED: Pythia calls the second weight
    //     MUR1.0_MUF2.0 while the header says 1002 is renscfact=2d0
    //     facscfact=1d0, i.e. (2, 1).
    // Using those names would assign the wrong scale point to every value and
    // lose a PDF set, and nothing would look wrong.  The index below lets the
    // analysis read the weights with the repository's OWN LHE parser, which
    // is already validated (analysis/powheg_nlo_uncertainty.py) and gets all
    // eleven.
    //
    // THE INDEX COMES FROM PYTHIA, NOT FROM COUNTING next() CALLS.  The
    // obvious counter is wrong: Pythia RETRIES internally, reading further
    // LHE events without next() returning false, so on this sample it
    // recorded 11 skips where the run actually consumed 940 extra events
    // (Pythia's own table: 50000 tried, 49060 accepted).  The join built on
    // it was silently PERMUTED -- the total weight still matched, because it
    // is the same multiset, while every selected subset got the wrong
    // weights.  info.nSelected() is Pythia's own count of LHE events read,
    // so minus one it is the 0-based position of the event just produced.
    evt.add_attribute("lhe_index",
                      std::make_shared<HepMC3::IntAttribute>(
                          (int)(pythia.info.nSelected() - 1)));

    // CARRY THE LHE VARIATION WEIGHTS INTO THE HepMC.
    //
    // POWHEG writes its scale and PDF variations into each <event> as an
    // <rwgt> block, and Pythia's LHEF reader parses them, but the converter
    // keeps only the nominal.  Without them a HADRON-LEVEL selection can
    // never be reweighted: the FASER tiers and the dimuon tag need the
    // shower and the decays, so they cannot be evaluated on the LHE, and the
    // LHE is where the weights were.  Propagating them here is what lets a
    // PDF or scale study be done UNDER a tier.
    //
    // Written only when there is more than one, so a sample generated from an
    // unreweighted LHE keeps exactly the one-weight records it always had and
    // nothing already published changes.  The names go on the run info, which
    // is where HepMC3 expects them and how the analysis maps a weight index
    // back to a PDF set.
    {
      int nw = pythia.info.getWeightsDetailedSize();
      if (nw > 1) {
        std::vector<std::string> names;
        std::vector<double> vals;
        names.reserve(nw); vals.reserve(nw);
        for (int iw = 0; iw < nw; ++iw) {
          names.push_back(pythia.info.weightNameByIndex(iw));
          vals.push_back(pythia.info.weightValueByIndex(iw));
        }
        // THE NAMES GO IN THE JSON SIDECAR, NOT THE HepMC RUN INFO.
        // WriterAscii emits the run-info header before the first event, and
        // the weight structure is only known once an event has been read, so
        // setting names on evt.run_info() here is too late -- it produced a
        // file with eleven unnamed weights, which is worse than none: the
        // analysis would have had to guess the order.  The sidecar is written
        // after the loop and the analysis already reads it.
        if (!runInfoNamed) { weightNames = names; runInfoNamed = true; }
        evt.weights() = vals;
      }
    }
    writer.write_event(evt);
    sumW += pythia.info.weight();
    ++nDone;
    if (maxEvents > 0 && nDone >= maxEvents) break;
  }
  writer.close();
  pythia.stat();

  // sigmaGen is in mb; weights are POWHEG's +-sigma_abs in pb
  {
    std::ofstream js(out + "_xsec.json");
    js << "{ \"sigma_gen_mb\": " << pythia.info.sigmaGen()
       << ", \"sigma_err_mb\": " << pythia.info.sigmaErr()
       << ", \"n_accepted\": " << nDone
       << ", \"n_aborted\": " << nAborted
       << ", \"sum_weights_pb\": " << sumW;
    if (!weightNames.empty()) {
      js << ", \"weight_names\": [";
      for (size_t i = 0; i < weightNames.size(); ++i)
        js << (i ? ", " : "") << "\"" << weightNames[i] << "\"";
      js << "]";
    }
    js << " }" << std::endl;
  }

  std::cout << std::flush;
  std::cerr << std::flush;
  std::_Exit(0);
}
