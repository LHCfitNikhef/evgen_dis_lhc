# Path configuration for the benchmark's Makefiles.
#
# Values come from config.sh -- the single source -- via tools/print_config.sh,
# so a Makefile never re-declares a default that could drift.  Include it as:
#
#     include $(dir $(lastword $(MAKEFILE_LIST)))../config.mk
#
# and override the same way as everywhere else:
#
#     BENCH_DEPS=/cvmfs/... make
BENCH_CONFIG_DIR := $(patsubst %/,%,$(dir $(lastword $(MAKEFILE_LIST))))
BENCH_CFG := $(BENCH_CONFIG_DIR)/tools/print_config.sh

# One variable per call, taking everything after the first "=".
#
# NOT a $(foreach)+$(eval) over all of them: make's $(shell) collapses newlines
# into spaces, so a value containing a space (BENCH_NOSLEEP is "caffeinate -i",
# and a path could have one) silently becomes two make words and the eval dies
# with "missing separator".
bench_cfg = $(shell $(BENCH_CFG) $(1) | cut -d= -f2-)

PYTHIA8_DIR := $(call bench_cfg,PYTHIA8_DIR)
BENCH_DEPS  := $(call bench_cfg,BENCH_DEPS)
GENIE_DIR   := $(call bench_cfg,GENIE_DIR)
HERWIG_BIN  := $(call bench_cfg,HERWIG_BIN)
NUWRO_DIR   := $(call bench_cfg,NUWRO_DIR)

# Compiler: a name, not a path, but just as unportable -- Apple clang on the
# laptop, g++ on the Linux cluster; config.sh picks by uname.
#
# ":=" from the config, NOT "CXX ?=": make's "?=" honours the environment, and
# conda exports CXX=arm64-apple-darwin20.0.0-clang++, whose availability model
# does not survive the macOS 26 SDK (see genie/genie_setup.sh) and which has
# already cost this project a build once by suppressing -std=c++11.
# Override with BENCH_CXX.
CXX := $(call bench_cfg,BENCH_CXX)
