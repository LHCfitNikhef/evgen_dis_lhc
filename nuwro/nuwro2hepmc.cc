// NuWro's ROOT event tree -> HepMC3, so that analysis/analyze_nu.py can read
// it with no change at all.
//
// WHY A CONVERTER AND NOT A SECOND PARSER.  Every number in this benchmark
// passes through `parse_hepmc3`, and CONVENTIONS.md rule 2b plus the parser-drift
// checker exist because that function was written out twice and the copies can
// diverge.  Writing a third reader for NuWro's own format would put a whole
// generator behind code no gate watches.  Converting instead means NuWro's
// events are analysed by exactly the same lines as POWHEG's.
//
// >>> AND THE FILE IS WRITTEN BY HepMC3 ITSELF, NOT BY hand. <<<  The ASCII
// format is ten whitespace-separated fields per particle and it would be easy
// to emit by printf -- and easy to emit subtly wrong, which is the silent
// failure this project keeps meeting.  WriterAscii cannot produce a file its
// own reader would reject.
//
// WHAT IS AND IS NOT CARRIED OVER:
//
//   beams        in[0] is the neutrino and in[1] the struck nucleon, both
//                written with status 4.  On a FREE-nucleon run in[1] is a
//                proton at rest, exactly as every other generator here; on a
//                NUCLEUS run it is a bound nucleon with Fermi momentum, so
//                the DIS invariants the analysis builds are then taken with
//                respect to the struck nucleon -- which is the right
//                definition, and a difference from the free-target samples
//                that has to be remembered when the two are compared.
//   final state  `post`, the particles LEAVING THE NUCLEUS -- i.e. after
//                final-state interactions.  On a free target it equals `out`.
//                This is the whole reason NuWro is in the comparison.
//   weight       NuWro stores the same total cross-section on every event
//                (equally weighted events), in cm^2.  It is written as
//                sigma[pb]/N per event, so summing the weights of the
//                selected events gives sigma_fid directly -- the convention
//                every other sample here uses.
//
//   D MESONS ARE NOT CARRIED OVER, and so no charm observable is available
//   from NuWro yet.  They are counted in this benchmark at PRODUCTION, from
//   status-2 records, and NuWro's decayed particles live in `all` with mother
//   indices rather than in a HepMC graph.  Attaching them to the primary
//   vertex would make the graph say something false; building the real chain
//   is a later job.  Until then the charm figures simply have no NuWro row.
//
// Usage: nuwro2hepmc <in.root> <out.hepmc>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "event1.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"

using namespace HepMC3;

// NuWro works in MeV (src/jednostki.h sets MeV = 1); HepMC3 is asked for GeV.
static const double MEV = 1.0e-3;
// 1 pb = 1e-36 cm^2.
static const double CM2_TO_PB = 1.0e36;

static GenParticlePtr mk(const particle &p, int status) {
    FourVector v(p.x * MEV, p.y * MEV, p.z * MEV, p.t * MEV);
    return std::make_shared<GenParticle>(v, p.pdg, status);
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: nuwro2hepmc <in.root> <out.hepmc>\n");
        return 2;
    }
    TFile *f = TFile::Open(argv[1]);
    if (!f || f->IsZombie()) {
        fprintf(stderr, "cannot open %s\n", argv[1]);
        return 1;
    }
    TTree *t = dynamic_cast<TTree *>(f->Get("treeout"));
    if (!t) {
        fprintf(stderr, "no treeout in %s\n", argv[1]);
        return 1;
    }
    event *e = new event();
    t->SetBranchAddress("e", &e);
    const Long64_t n = t->GetEntries();
    if (n <= 0) {
        fprintf(stderr, "%s holds no events\n", argv[1]);
        return 1;
    }

    // THE PER-EVENT WEIGHT NEEDS THE EVENT COUNT, so the first entry is read
    // before the loop rather than the cross-section being taken from the log.
    t->GetEntry(0);
    const double sigma_pb = e->weight * CM2_TO_PB;
    const double w = sigma_pb / double(n);

    WriterAscii out(argv[2]);
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        // >>> EVERY EVENT MUST CARRY THE SAME CROSS-SECTION. <<<  NuWro's
        // events are equally weighted and `weight` is documented as "set to
        // total cross section on saving to file".  If that ever stops being
        // true the per-event weight computed above is wrong for every event
        // but the first, and nothing downstream could tell -- so it is
        // checked here rather than assumed.
        if (e->weight * CM2_TO_PB != sigma_pb) {
            fprintf(stderr, "event %lld carries sigma %g pb, event 0 carried "
                            "%g -- these events are not equally weighted and "
                            "the conversion would misnormalise them\n",
                    (long long)i, e->weight * CM2_TO_PB, sigma_pb);
            return 1;
        }
        GenEvent evt(Units::GEV, Units::MM);
        evt.set_event_number(int(i));
        evt.add_attribute("alphaQCD", std::make_shared<DoubleAttribute>(0.0));
        evt.weights().push_back(w);
        GenVertexPtr v = std::make_shared<GenVertex>();
        for (size_t k = 0; k < e->in.size(); ++k)
            v->add_particle_in(mk(e->in[k], 4));
        for (size_t k = 0; k < e->post.size(); ++k)
            v->add_particle_out(mk(e->post[k], 1));
        evt.add_vertex(v);
        out.write_event(evt);
    }
    out.close();
    f->Close();

    // >>> AND A SIDECAR THE ANALYSIS CAN READ. <<<  Every other sample here
    // hands its normalisation to analyze_nu.py through a file rather than a
    // log line, and for a good reason: a log is written for a human and its
    // format changes without anyone noticing.
    //
    // THE EVENT COUNT IS IN IT BECAUSE THE CLOSURE GATE CANNOT HELP HERE.
    // That gate compares the delivered sum of weights against the
    // integrator's cross-section -- but this converter DEFINES the weight as
    // sigma/N, so the two agree by construction and the gate is vacuous for
    // NuWro.  What can still go wrong is a truncated file, so the count is
    // recorded and the analysis checks it.
    std::string js(argv[2]);
    const size_t dot = js.rfind(".hepmc");
    js = (dot == std::string::npos ? js : js.substr(0, dot)) + "_xsec.json";
    FILE *jf = fopen(js.c_str(), "w");
    if (!jf) {
        fprintf(stderr, "cannot write %s\n", js.c_str());
        return 1;
    }
    fprintf(jf, "{ \"generator\": \"nuwro\", \"sigma_pb\": %.10g,"
                " \"weight_pb_per_event\": %.10g, \"n_events\": %lld }\n",
            sigma_pb, w, (long long)n);
    fclose(jf);

    printf("%lld events, sigma = %.6g pb, weight %.6g pb/event -> %s\n",
           (long long)n, sigma_pb, w, argv[2]);
    return 0;
}
