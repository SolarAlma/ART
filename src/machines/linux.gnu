# Setup for vanila Ubuntu with HDF installed via
# sudo apt install libhdf5-mpi-dev
# 
# Note that hd5lib will be placed in: /usr/lib/aarch64-linux-gnu/hdf5/openmpi/ 
#                                 or: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/


# HDF5_DIR  = #setup by module

OPTFLAGS  = -O2 -mtune=native

CXX       = mpicxx
MPI_CXX   = mpicxx

RTINCLUDE =
CXXFLAGS  = $(OPTFLAGS) -std=c++11 -I/usr/include/hdf5/openmpi/

LDFLAGS   = -L/usr/lib/aarch64-linux-gnu/hdf5/openmpi/ -lhdf5 -lhdf5_hl
DEBUG	  = -g

