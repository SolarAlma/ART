# Setup for MacBook Pro with OsX Sierra + PGI 16.0.1 compilers with build in MPICH 

# HDF5_DIR  = #setup by module 

CXX       = mpic++
FC        = mpifort

MPI_CXX   = mpic++
MPI_FC    = mpifort

RTINCLUDE = 
OPTFLAGS  =  
F90FLAGS  = 
CXXFLAGS  = $(OPTFLAGS) -std=c++11

LDFLAGS   = -lhdf5 -lhdf5_hl -lz -lsz -lstdc++ -lc++
DEBUG	    =

