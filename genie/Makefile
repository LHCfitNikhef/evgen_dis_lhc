# Build the GHEP -> HepMC3 converter for the muon-DIS benchmark.
# Requires: source genie/genie_setup.sh first (GENIE, ROOT, deps on paths).

include $(dir $(lastword $(MAKEFILE_LIST)))../config.mk
GENIE ?= $(GENIE_DIR)
CONDA  := $(BENCH_DEPS)

CXXFLAGS = $(shell root-config --cflags) -I$(GENIE)/src -I$(CONDA)/include -O2 -Wall
LDFLAGS  = $(shell genie-config --libs) $(shell root-config --libs) \
           -L$(CONDA)/lib -lHepMC3 -Wl,-rpath,$(CONDA)/lib -Wl,-rpath,$(GENIE)/lib

gtohepmc3: gtohepmc3.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f gtohepmc3
