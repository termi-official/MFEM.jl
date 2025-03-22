JL_SHARE = $(shell julia --project -e 'print(joinpath(Sys.BINDIR, Base.DATAROOTDIR, "julia"))')
CXXFLAGS += $(patsubst -std=gnu%,,$(shell $(JL_SHARE)/julia-config.jl --cflags))
CXXFLAGS += -DJLCXX_USE_TYPE_MAP
CXXFLAGS += -DVERBOSE_IMPORT #To enable the verbose mode of the libray loading
# CXXFLAGS += -Wall -O0 -g     #To compile with debugger infomation
# CXXFLAGS += -O0 -g     #To compile with debugger infomation
CXXWRAP_PREFIX = $(shell julia --project -e "import CxxWrap; print(CxxWrap.prefix_path())")
LDLIBS    	   = $(shell $(JL_SHARE)/julia-config.jl --ldlibs) -L$(CXXWRAP_PREFIX)/lib -lcxxwrap_julia -lcxxwrap_julia_stl
LDLIBS        += $(shell $(JL_SHARE)/julia-config.jl --ldlibs)
LDFLAGS       += $(shell $(JL_SHARE)/julia-config.jl --ldflags)
CXXWRAP_CPPFLAGS=-I $(shell echo $(CXXWRAP_PREFIX)"/include") -I $(shell pwd) -std=c++20

MFEM_LIBS?=$(shell julia --project -e 'using MFEM_seq_jll; print(MFEM_seq_jll.libmfem_path)')
MFEM_INC_DIR?=$(shell julia --project -e 'using MFEM_seq_jll; print(joinpath(MFEM_seq_jll.artifact_dir, "include"))')

CPPFLAGS += -MMD -I $(MFEM_INC_DIR) -I ./libMFEM/src

LINK.o = $(CXX) $(LDFLAGS) $(TARGET_ARCH)

.PHONY: all check_mfem make_lib_from_objs

# .PRECIOUS: $(GENERATED_CXX) MFEM.wit

PRODUCTS=libMFEM/libMFEM.so
GENERATED_CXX:=$(file < libMFEM/src/generated_cxx)
OBJS=$(addprefix libMFEM/build/, $(patsubst %.cxx,%.o, $(GENERATED_CXX)))
# DEPS=$(patsubst %.o,%.d, $(OBJS))

all: $(PRODUCTS)

libMFEM/templates.h:
	-mkdir -p libMFEM/build
	echo "Makefiles are painful."

check_mfem:
ifeq ($(MFEM_LIBS),)
	$(error ERROR: "MFEM_LIBS not set.")
endif
ifeq ($(MFEM_INC_DIR),)
	$(error ERROR: "MFEM_INC_DIR not set.")
endif

libMFEM/build/%.o: libMFEM/src/%.cxx
	$(COMPILE.cc) $(CXXWRAP_CPPFLAGS) -o $@ $<

make_lib_from_objs: $(OBJS)
	$(LINK.o) -o libMFEM/libMFEM.so --shared -fPIC $(OBJS) $(MFEM_LIBS) $(LDLIBS) 

libMFEM/libMFEM.so:
	$(MAKE) check_mfem
	$(MAKE) make_lib_from_objs

# -include $(DEPS)
