# Build the POWHEG-LHE shower driver for the muon-DIS benchmark
include $(dir $(lastword $(MAKEFILE_LIST)))../config.mk
PYTHIA := $(PYTHIA8_DIR)
CONDA  := $(BENCH_DEPS)

main_powheg: main_powheg.cc DISPowhegHooks.h
	$(CXX) $< -o $@ -O2 -std=c++11 -pthread -I. \
	  -I$(PYTHIA)/include -I$(CONDA)/include \
	  -L$(PYTHIA)/lib -lpythia8 -L$(CONDA)/lib -lHepMC3 \
	  -Wl,-rpath,$(PYTHIA)/lib -Wl,-rpath,$(CONDA)/lib

clean:
	rm -f main_powheg
