# Build the LeptonicDISCut ThePEG plugin as a standalone .so.
#
# Deliberately a PLUGIN and not a patch to ThePEG's own SimpleDISCut: the
# working Herwig install is never rebuilt, so a mistake here cannot break the
# generator the rest of the benchmark depends on.  Load it from a card with
#     create ThePEG::LeptonicDISCut /Herwig/Cuts/LeptonicCut LeptonicDISCut.so
include $(dir $(lastword $(MAKEFILE_LIST)))../../config.mk
TPCONFIG ?= $(HERWIG_BIN)/thepeg-config
CXXFLAGS += -O2 -fPIC -std=c++14 $(shell $(TPCONFIG) --cppflags)
TPLIBDIR := $(shell $(TPCONFIG) --libdir)
# A ThePEG plugin is dlopen'd into a host that already provides the ThePEG
# symbols, so link against libThePEG explicitly and add an rpath -- on macOS a
# .so must resolve every symbol at link time, unlike ELF.
LDFLAGS  += -shared $(shell $(TPCONFIG) --ldflags) \
            -L$(TPLIBDIR) -lThePEG -Wl,-rpath,$(TPLIBDIR)

all: LeptonicDISCut.so

LeptonicDISCut.so: LeptonicDISCut.cc LeptonicDISCut.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ LeptonicDISCut.cc

clean:
	rm -f LeptonicDISCut.so *.o
