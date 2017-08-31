/* --------------------------------------------------------------------

   I/O routines that read the input file and HDF5 related routines
   Written by J. de la Cruz Rodriguez (ISP-SU 2016)

   Depends on HDF5 C routines with parallel MPI I/O installed

   -------------------------------------------------------------------- */


#ifndef __IO_H__
#define __IO_H__

#include "mtypes.h"
#include "cmemt2.h"
#include "model.h"

#include <string>
#include <cstdio>

/* ---------------------------------------------------------------------------- */

double air2vac(const double lambda_air);
double vac2air(const double lambda_vac);
void readInput(std::string filename, info &inp, bool verbose=true);
std::vector<std::string> strsplit(std::string &var, std::string token, bool rmspaces = true);
std::string removeSpaces(std::string input);
void initReadIO(info &inp);
void h5err(int myrank, const char *rname);
bool bdir_exists( std::string name);
bool h5varExists(hid_t &fid, std::string vname, hid_t &vid, FILE *log = stderr);
void readAtmosTYX(size_t tt, size_t yy, size_t xx, modl::mdepth &m, info &inp);
void readPix(hid_t &hid, hsize_t msp, hsize_t dsp, double *var, FILE *log);
void closeModel(info &inp);
void closeProfile(info &inp);
void writeProfileTYX(size_t tt, size_t yy, size_t xx, double *sp, info &inp);
void initWriteIO(info &inp, double *lambda);
int readValdLines(std::string filename, info &input);


#endif
