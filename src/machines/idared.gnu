# Setup for MacBook Pro with OsX Sierra + PGI 16.0.1 compilers with build in MPICH

# HDF5_DIR  = #setup by module

CXX       = mpic++
FC        = mpifort

MPI_CXX   = mpic++
MPI_FC    = mpifort

RTINCLUDE =
OPTFLAGS  =  -I/usr/local/opt/hdf5-parallel/include/
F90FLAGS  =
CXXFLAGS  = $(OPTFLAGS) -std=c++11

LDFLAGS   = -lc++ -L/usr/local/opt/hdf5-parallel/lib -lhdf5 -lhdf5_hl
DEBUG	    =

