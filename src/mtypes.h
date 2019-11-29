#ifndef __MTYPES_H__
#define __MTYPES_H__

#include <string>
#include <vector>
#include "hdf5.h"
#include "mpi.h"

static const int nlevpf = 6;


/* ---------------------------------------------------------------------------- */

struct species{
  int anum, ion, nLev;
  double n_tot[nlevpf];
  double pf[nlevpf];
  double eion[nlevpf];
};

/* ---------------------------------------------------------------------------- */

struct line{
  char elem[8], label[15];
  double Jup, Jlow, Gup, Glow;
  double w0, nu0, width, eion, nu_min, nu_max;
  double gf, e_low, e_up, amass;
  double g_rad, g_str, g_vdw;
  double b_sig, b_alp, b_vbar, b_gvw;
  int anum, ion,off, idx;
  bool barklem, firsttime;
  
  // Zeeman splitting
  int nZ;
  std::vector<double> strength, splitting;
  std::vector<int> iL;
};
typedef line line_t;


/* ---------------------------------------------------------------------------- */

struct region{
  double w0, dw, cscal;
  int nw, off;
  std::vector<double> wav, nu;
  std::vector<int> idx;
  std::string inst, ifile;
};
typedef region region_t;

/* ---------------------------------------------------------------------------- */

struct h5model{
  std::string filename;
  std::vector<double> x, y, z;
  int nx, ny, ndep, nt;
  hid_t vid[13];
  hid_t fid, did, mid;
};

/* ---------------------------------------------------------------------------- */

struct h5prof{
  std::string filename;
  std::vector<double> x, y, z;
  int nx, ny, ndep, nt;
  int nvid = 2;
  int opt_nvid = 1;
  hid_t vid[5]; // Stokes_I + tau = 1
  hid_t opt_vid[1]; // Optional values, 1 for Contribution 
  // hid_t vid[4]; 4 Stokes parameters + formation heigt
  hid_t fid, did[5], opt_did, mid, opt_mid;
};

/* ---------------------------------------------------------------------------- */

struct info{
  int nx, ny, ndep, nt, nw, nstokes, solver, myrank, nproc, verbose, units, eos_type,
    getContrib;
  double mu, temperature_cut, gravity, dlam;
  size_t ipix;
  std::string hostname;
  std::string lines_file;
  std::vector<bool> vdef;
  std::vector<region> reg;
  std::vector<line> lin;
  h5model m;
  h5prof p;
  MPI_Comm comm;
  MPI_Info info;
  FILE *log;
};

/* ---------------------------------------------------------------------------- */




/* ---------------------------------------------------------------------------- */

#endif
