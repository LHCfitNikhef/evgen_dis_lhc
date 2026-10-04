# Build the Pythia8 HepMC3 driver for the muon-DIS benchmark
include $(dir $(lastword $(MAKEFILE_LIST)))../config.mk
PYTHIA := $(PYTHIA8_DIR)
CONDA  := $(BENCH_DEPS)

main_dis: main_dis.cc
	$(CXX) $< -o $@ -O2 -std=c++11 -pthread \
	  -I$(PYTHIA)/include -I$(CONDA)/include \
	  -L$(PYTHIA)/lib -lpythia8 -L$(CONDA)/lib -lHepMC3 \
	  -Wl,-rpath,$(PYTHIA)/lib -Wl,-rpath,$(CONDA)/lib

clean:
	rm -f main_dis
