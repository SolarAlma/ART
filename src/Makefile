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

### filling missing parts

MPI_CXX   ?= mpic++
HDF5_DIR  ?=
RTINCLUDE ?=
OPTFLAGS  ?= -O2
CXXFLAGS  ?= $(OPTFLAGS) -std=c++11
LDFLAGS   ?= -lhdf5 -lhdf5_hl -lz
DEBUG     ?=

INFO      += compiler

### define object files list

OFILES = cop.o  eoswrap.o witt.o model.o io.o clte.o main.o


### compilation

.SUFFIXES: .cc .o
.PHONY: clean

%.o: %.cc *.h
	$(MPI_CXX) $(CXXFLAGS) $(DEBUG) $(RTINCLUDE) -c $< -o $@

compile: $(INFO) $(OFILES)
	$(MPI_CXX) $(RTINCLUDE) -o RTma.x $(OFILES) $(LDFLAGS)

test: $(INFO) $(TFILES)
	$(MPI_CXX) $(RTINCLUDE) -o test.x $(TFILES) $(LDFLAGS)

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
	@echo


compiler:
	@echo " MPI_CXX      = " $(MPI_CXX)
	@echo " HDF5_DIR     = " $(HDF5_DIR)
	@echo " RTINCLUDE    = " $(RTINCLUDE)
	@echo " OPTFLAGS     = " $(OPTFLAGS)
	@echo " CXXFLAGS     = " $(CXXFLAGS)
	@echo " LDFLAGS      = " $(LDFLAGS)
	@echo " DEBUG        = " $(DEBUG)
	@echo
