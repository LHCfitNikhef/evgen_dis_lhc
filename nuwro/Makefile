# Build the NuWro -> HepMC3 converter.
#
# It links against NuWro's OWN object files -- the same six that
# src/nuwro2rootracker links -- because the ROOT tree holds a `event` object
# and reading it needs that class and its dictionary.  Nothing is copied out
# of the generator tree (CONVENTIONS.md rule 1); only this source is ours.
include $(dir $(lastword $(MAKEFILE_LIST)))../config.mk
NW    := $(NUWRO_DIR)
CONDA := $(BENCH_DEPS)
NWOBJ := $(NW)/src/event1.o $(NW)/src/pdg.o $(NW)/src/particle.o \
         $(NW)/src/generatormt.o $(NW)/src/dirs.o $(NW)/src/event1Dict.o

nuwro2hepmc: nuwro2hepmc.cc
	$(CXX) $< -o $@ -O2 -std=c++17 -I$(NW)/src -I$(CONDA)/include \
	  $(shell root-config --cflags) $(NWOBJ) \
	  $(shell root-config --libs) -L$(CONDA)/lib -lHepMC3 \
	  -Wl,-rpath,$(CONDA)/lib

clean:
	rm -f nuwro2hepmc
