# Setup for MacBook Pro with OsX // GNU + OpenMPI installed via home brew.

# HDF5_DIR  = #setup by module

OPTFLAGS  = -O2 -mtune=native

CXX       = mpicxx
MPI_CXX   = mpicxx

RTINCLUDE =
CXXFLAGS  = $(OPTFLAGS) -std=c++11 -I/opt/homebrew/include/

LDFLAGS   = -lstdc++ -lc++ -L/opt/homebrew/lib -lhdf5 -lhdf5_hl
DEBUG	    = -g

