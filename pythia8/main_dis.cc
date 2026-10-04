// Pythia8 driver for the muon-DIS benchmark: read a .cmnd card, generate,
// write HepMC3 events plus a JSON sidecar with the generated cross-section.
// Usage: ./main_dis card.cmnd outprefix [nEvents] [seed]
//
// The writer is closed explicitly and the process exits via _Exit: with the
// conda-built HepMC3/LHAPDF next to Apple clang, static destructors abort at
// teardown (exit 134) even though all output is complete.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include <cstdlib>
#include <fstream>

using namespace Pythia8;

int main(int argc, char* argv[]) {

  if (argc < 3) {
    cerr << "usage: " << argv[0] << " card.cmnd outprefix [nEvents] [seed]"
         << endl;
    return 1;
  }
  std::string card = argv[1], out = argv[2];

  Pythia pythia;
  pythia.readFile(card);
  if (argc > 3) pythia.readString("Main:numberOfEvents = "
                                  + std::string(argv[3]));
  if (argc > 4) {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::string(argv[4]));
  }
  int nEvent = pythia.mode("Main:numberOfEvents");

  if (!pythia.init()) {
    cerr << "Pythia init failed" << endl;
    return 1;
  }

  HepMC3::Pythia8ToHepMC3 converter;
  HepMC3::WriterAscii writer(out + ".hepmc");

  // tolerate a ~1% failure rate (DIS dipole recoil occasionally fails at
  // parton/hadron level); failed events are regenerated
  int nAbort = std::max(50, nEvent / 100), iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) {
      if (++iAbort > nAbort) { cerr << "too many aborts" << endl; break; }
      --iEvent;
      continue;
    }
    HepMC3::GenEvent evt;
    converter.fill_next_event(pythia, &evt);
    writer.write_event(evt);
  }
  writer.close();
  pythia.stat();

  // sigmaGen is in mb; record it for the analysis normalisation
  {
    std::ofstream js(out + "_xsec.json");
    js << "{ \"sigma_gen_mb\": " << pythia.info.sigmaGen()
       << ", \"sigma_err_mb\": " << pythia.info.sigmaErr()
       << ", \"n_accepted\": " << pythia.info.nAccepted() << " }" << std::endl;
  }

  std::cout << std::flush;
  std::cerr << std::flush;
  std::_Exit(0);
}
