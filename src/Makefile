# See LICENSE.txt for license details.

CXX=g++
CXX_FLAGS += -std=c++11 -Wall -march=native -g -Isrc/lib/ -Lsrc/lib/ -I boost_1_76_0

PAR_FLAG = -fopenmp
OPT_FLAG = -O3
DEBUG_FLAG = -fkeep-inline-functions

ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG = -openmp
endif

ifneq (,$(findstring sunCC,$(CXX)))
	CXX_FLAGS = -std=c++11 -xO3 -m64 -xtarget=native
	PAR_FLAG = -xopenmp
endif
ifneq ($(DEBUG),1)
	CXX_FLAGS += $(OPT_FLAG)
else
	CXX_FLAGS += $(DEBUG_FLAG) 
endif

ifneq ($(SERIAL), 1)
	CXX_FLAGS += $(PAR_FLAG)
endif
KERNELS = tc_doulion tc_colorful tc_bf tc_kh tc_1h tc_base tc_pgp tc_redex tc_auto_approx_v1 tc_auto_approx_v2 jp-cn_bf jp-cn_kh jp-cn_1h jp-cn_base jp-cn_pgp jp-cn_redex jp-jc_bf jp-jc_kh jp-jc_1h jp-jc_base jp-jc_pgp jp-jc_redex jp-ov_bf jp-ov_kh jp-ov_1h jp-ov_base jp-ov_pgp jp-ov_redex 4c_bf 4c_bf 4c_1h 4c_kh 4c_base 4c_pgp 4c_redex

SUITE = $(KERNELS)

.PHONY: all
all: $(SUITE)

tc_set_simd: src/tc_set_simd.cc src/*.h src/*.hpp
	$(CXX) $(CXX_FLAGS) $< src/set_operation.cpp src/util.cpp -o $@

% : src/algorithms/%.cc src/*.h src/*.hpp
	$(CXX) $(CXX_FLAGS) $< src/set_operation.cpp src/sets.hpp src/util.cpp src/MurmurHash3.cpp -o $@

# g
# include test/test.mk

# Benchmark Automation
#include benchmark/bench.mk


.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*
