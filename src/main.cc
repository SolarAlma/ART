/* --------------------------------------------------------------

  Radiative transfer code for ALMA applications
  Written by J. de la Cruz Rodriguez (ISP-SU 2016)


  Uses parallel HDF5 I/O and MPI to compute the continuum
  intensity at a given number of wavelengths


  Changelog:
      2017-08-29, JdlCR: EOS is now a baseclass and there are 
      two options "pisk" and "witt". Added support for VALD 
      input.
  


  TODO: the EOS is slow. It is probably an overkill for
  solar applications because it can handle complex molecules.
  Should change it to Wittmann or similar.

  
  

  -------------------------------------------------------------- */

#include <iostream>
#include <cstdio>
#include <vector>
#include <mpi.h>
//
#include "cmemt2.h"
#include "mtypes.h"
#include "io.h"
#include "model.h"
#include "clte.h"
//
#include "eoswrap.h"
#include "witt.h"
#include "ceos.h"

using namespace std;
using namespace modl;
  
/* ---------------------------------------------------------------------------- */

void processData(info &inp)
{
  
  /* --- open atmos file and init dimensions --- */

  initReadIO(inp);
  bool master = (inp.myrank == 0) ? true : false;


  
  /* --- Init EOS  --- */

  eoswrap *EoS = NULL;
  if(inp.eos_type == 0) EoS = new eos::witt(inp.lin, 0, NULL, inp.gravity);
  else if(inp.eos_type == 1) EoS = new ceos(inp.lin, 0, NULL, inp.gravity);
  else{
    fprintf(inp.log,"main: ERROR, accepted values for eos_type are: [pisk] and [witt]\n");
    exit(0);
  }


  /* --- Init Solver --- */
  
  clte atmos(inp.reg, inp.lin);
  std::vector<double> synthetic;
  synthetic.resize(atmos.nw);

  /* --- Init output file --- */

  initWriteIO(inp, &atmos.lambda[0]);
  
  
  
  /* --- loop pixels and synthesize --- */

  size_t npix = inp.nx * inp.ny * inp.nt, step = inp.nproc, i0 = inp.myrank, per = 0, oper = 0;
  mdepth m(inp.ndep, true);

  
  for(size_t ipix = i0; ipix<npix; ipix += step){

    /* --- Compute pixel offset --- */
    
    size_t tt = ipix / (inp.ny*inp.nx);
    size_t yy = (ipix - tt*inp.nx*inp.ny) / (inp.nx);
    size_t xx = (ipix - tt*inp.nx*inp.ny - yy*inp.nx);
    

    /* --- read Atmos from file --- */
    
    readAtmosTYX(tt, yy, xx, m, inp);
    m.fillDensities(inp, *EoS);

    
    /* --- get spectrum --- */

    atmos.synth_nonpol(m, &synthetic[0], *EoS, inp.mu, inp.solver, (bool)inp.getContrib);
    
    //exit(0);
    
    /* --- Write to disk --- */

    writeProfileTYX(tt, yy, xx, &synthetic[0], inp);
    writeTau1TYX(tt,yy,xx, &atmos.tau_eq_1_z[0], inp);

    if ((bool)inp.getContrib) writeCFuncTYX(tt,yy,xx, &atmos.C[0], inp);

    /* --- printout some info --- */

    per = ipix*100./(std::max(npix-1.0,1.0));
    if(oper != per){
      oper = per;
      if(master) fprintf(stderr,"\rinfo: processing %d%s", int(per), "%");
    }
  }
  if(master) fprintf(stderr,"\rinfo: processing %d%s\n", 100, "%");


  
  /* --- close files --- */

  closeModel(inp);
  closeProfile(inp);

  
}

static void show_usage(std::string name) {
  std::cerr << "\nUsage: " << name << " option(s):\n"
            << "Options:\n"
            << "\t-h,--help\t\tShow this help message\n"
            << "\t-i,--input INPUT\tSpecify the input file\n"
            << std::endl;
}

/* ---------------------------------------------------------------------------- */


int main(int narg, char *argv[]) {

  info inp = {};
  int  hlen = 0;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Init(&narg, &argv);
  inp.comm = MPI_COMM_WORLD;
  inp.info = MPI_INFO_NULL;

  
  /* --- init MPI --- */
  
  MPI_Comm_rank(inp.comm, &inp.myrank);      // Job number
  MPI_Comm_size(inp.comm, &inp.nproc);       // number of processors
  MPI_Get_processor_name(&hostname[0], &hlen); // Hostname
  inp.hostname = string(hostname);
  
  /* --- Read input --- */

  std::string input_file = "input.cfg";

  if (narg < 2) {
    if (inp.myrank == 0) std::cout << "Using default input file: " << input_file << std::endl;
  } else {
    for (size_t i = 0; i < narg; ++i) {
      std::string arg = argv[i];
      if ((arg == "-h") || (arg == "--help")) {
        show_usage(argv[0]);
      } else if ((arg == "-i") || (arg == "--input")) {
        if (i + 1 < narg) { // make sure this is not a last argument
          input_file = argv[++i];
          if (inp.myrank == 0) std::cout << "Using input file: " << input_file << std::endl;
        } else {
          if (inp.myrank == 0) std::cerr << "-i,--input option requires one argument." << std::endl;
          MPI_Abort(MPI_COMM_WORLD, 2);
        }
      }
    }
  }

  readInput(input_file, inp, ((inp.myrank == 0) ? true : false));
  fprintf(stderr,"log target %d\n", inp.log);
  MPI_Barrier(MPI_COMM_WORLD);

  fprintf(inp.log,"test log %d/%d\n", inp.myrank, inp.nproc);

  if(inp.lines_file != "") readValdLines(inp.lines_file, inp);
  
  /* --- Do slave or master based on process number --- */
  processData(inp);
  
  /* --- Finish MPI jobs --- */
  
  MPI_Barrier(MPI_COMM_WORLD); // Wait until all processors reach this point
  MPI_Finalize();

}

/* ---------------------------------------------------------------------------- */
