### find machine config file

OPSYS        ?= $(shell uname -s)
MACHINE      ?= $(shell uname -n)
ARCHF        ?= $(shell uname -a)
MACHINEMAKE   = machines/$(MACHINE)
MACHINEEXISTS = $(shell ls machines)
INFO          = machine_info

ifeq ("$(findstring $(MACHINE),$(MACHINEEXISTS))","")
	INFO      += no_machine
else
	include $(MACHINEMAKE)
endif

ifneq (,$(findstring Darwin,$(ARCHF))) #
  MACOSX_DEPLOYMENT_TARGET=10.6
  export MACOSX_DEPLOYMENT_TARGET
endif

### filling missing parts

MPI_CXX  ?= mpic++
MPI_FC   ?= mpif90
HDF5_DIR ?=
INCLUDE  ?=
OPTFLAGS ?= -02
F90FLAGS ?= $(OPTFLAGS)
CXXFLAGS ?= $(OPTFLAGS) -std=c++11 
LDFLAGS  ?= -lhdf5 -lhdf5_hl -lz 
DEBUG    ?=

INFO     += compiler

### define object files list

OFILES = eos_math_special.o eos_eqns.o eos.o cop.o ceos.o model.o io.o clte.o main.o

### compilation

.SUFFIXES: .cc .f .o
.PHONY: clean

.f.o:
	$(MPI_FC) $(F90FLAGS) $(DEBUG) -c $< -o $@
.cc.o:
	$(MPI_CXX) $(CXXFLAGS) $(DEBUG) $(INCLUDE) -c $< -o $@

compile: $(INFO) $(OFILES)
	$(MPI_FC) $(F90FLAGS) $(INCLUDE) -o RTma.x *.o $(LDFLAGS)

clean:
	rm -f *.o *.mod

Default: compile

machine_info:
	@echo
	@echo "### Machine name: $(MACHINE) "
	@echo

no_machine:
	@echo "[!] There is no specific makefile for this machine"
	@echo "[!] Trying with default setup"
	@echo $(MACHINEEXISTS)
	@echo $(MACHINEMAKE)
	@echo $(findstring $(MACHINE),$(MACHINEEXISTS))


compiler:
	@echo " MPI_CXX      = " $(MPI_CXX)
	@echo " MPI_FC       = " $(MPI_FC)
	@echo " HDF5_DIR     = " $(HDF5_DIR)	
	@echo " INCLUDE      = " $(INCLUDE)
	@echo " OPTFLAGS     = " $(OPTFLAGS)
	@echo " F90FLAGS     = " $(F90FLAGS)
	@echo " CXXFLAGS     = " $(CXXFLAGS)
	@echo " LDFLAGS      = " $(LDFLAGS)
	@echo " DEBUG        = " $(DEBUG)
	@echo
