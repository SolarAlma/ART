ARCH     = $(shell uname -m)
ARCHF    = $(shell uname -a)

MDIR = /home/jaime/local/parallel_hdf5/
OMPI_CXX = icpc
OMPI_FC = ifort

CXX = mpicxx
FC = mpif90

INCLUDE = -I$(MDIR)/include
LIBS =  -lstdc++  -L$(MDIR)/lib/ -lhdf5 -lhdf5_hl 

OPT = -Ofast -xCORE-AVX-I # CERTAIN VERSIONS OF IFORT BREAK THE EOS WITH AVX2 optimization!!
CPPFLAGS  = -std=c++11  $(OPT) $(INCLUDE)
FFLAGS = $(OPT) # Do not use ffast-math with the EOS!

OFILES = eos_math_special.o eos_eqns.o eos.o cop.o ceos.o model.o io.o clte.o main.o

ifneq (,$(findstring Darwin,$(ARCHF))) #
  MACOSX_DEPLOYMENT_TARGET=10.6
  export MACOSX_DEPLOYMENT_TARGET

endif 


.SUFFIXES: .o .f90 .cc

%.o: %.cc *.h
	$(CXX) $(CPPFLAGS) -c  $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c  $< -o $@


all: $(OFILES)
	$(FC) $(FFLAGS) -o RTma.x *.o $(INCLUDE) $(LIBS) -nofor_main


clean:
	rm -f *.o *.mod
