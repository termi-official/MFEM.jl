JL_SHARE = $(shell julia -e 'print(joinpath(Sys.BINDIR, Base.DATAROOTDIR, "julia"))')
CXXFLAGS += $(patsubst -std=gnu%,,$(shell $(JL_SHARE)/julia-config.jl --cflags))
# CXXFLAGS += -DVERBOSE_IMPORT #To enable the verbose mode of the libray loading
# CXXFLAGS += -Wall -O0 -g     #To compile with debugger infomation
CXXFLAGS += -O0 -g     #To compile with debugger infomation
LDFLAGS  += $(shell $(JL_SHARE)/julia-config.jl --ldflags)
LDLIBS   += $(shell $(JL_SHARE)/julia-config.jl --ldlibs)

CXXWRAP_CPPFLAGS=-I $(shell julia -e 'using CxxWrap; print(CxxWrap.prefix_path() * "/include")') -I $(shell pwd) -std=c++17


CPPFLAGS += -MMD -I $(MFEM_INC_DIR)

MFEM_LIBS?=/home/dogiermann/Builds/mfem-test/install/lib/libmfem.a
MFEM_INC_DIR?=/home/dogiermann/Builds/mfem-test/install/include/

LINK.o = $(CXX) $(LDFLAGS) $(TARGET_ARCH)

.PHONY: all clean run_demo check_mfem make_lib_from_objs

.PRECIOUS: $(GENERATED_CXX) MFEM-generated.wit

PRODUCTS=libMFEM/libjlMFEM.so
GENERATED_CXX:=$(file < libMFEM/src/generated_cxx)
OBJS=$(addprefix libMFEM/build/, $(patsubst %.cxx,%.o, $(GENERATED_CXX)))
DEPS=$(patsubst %.o,%.d, $(OBJS))

all: $(PRODUCTS)

clean:
	-$(RM) $(PRODUCTS) jlMFEM-report.txt jlMFEM.o jlMFEM.d $(GENERATED_CXX) MFEM-generated.jl MFEM-generated.wit $(OBJS) \
libMFEM/src/generated_cxx

MFEM-generated.wit: MFEM.wit
	$(MAKE) check_mfem
	$(shell echo "#\n# Warning: file generated automatically from $<\n#" > $@)
	$(shell sed "s@%MFEM_INC_DIR%@$(MFEM_INC_DIR)@" $< >> $@ || rm $@)

check_mfem:
ifeq ($(MFEM_LIBS),)
	$(error ERROR: "MFEM_LIBS not set.")
endif
ifeq ($(MFEM_INC_DIR),)
	$(error ERROR: "MFEM_INC_DIR not set.")
endif

run_demo: all
	. ./setup.sh && julia -i demo_MFEM.jl

test: all
	. ./setup.sh && julia demo_MFEM.jl

libMFEM/src/generated_cxx: MFEM-generated.wit jlMFEM-veto.h
	/home/dogiermann/Repos/wrapit/build/wrapit -v 10 --force $<
	$(eval GENERATED_CXX:=$(file < libMFEM/src/generated_cxx))
	$(eval OBJS:=$(addprefix libMFEM/build/, $(patsubst %.cxx,%.o, $(GENERATED_CXX))))

libMFEM/build/%.o: libMFEM/src/%.cxx
	-mkdir -p libMFEM/build
	$(COMPILE.cc) $(CXXWRAP_CPPFLAGS) -o $@ $<

libMFEM/libjlMFEM.so: libMFEM/src/generated_cxx
	$(MAKE) check_mfem
	$(MAKE) make_lib_from_objs

make_lib_from_objs: $(OBJS)
	$(LINK.o) -o libMFEM/libjlMFEM.so --shared -fPIC $(OBJS) $(MFEM_LIBS)

echo_%:
	@echo "$* = $(subst ",\",$($*))"

-include $(DEPS)
