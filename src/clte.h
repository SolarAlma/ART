/* --------------------------------------------------------------

   LTE class, only continuum opacities for the time being
   Written by J. de la Cruz Rodriguez (ISP-SU 2016)

   NOTES: 
   Unpolarized Bezier-3 solver for the time being.
   
   Requires ceos class to retrieve the partial densities of 
   background absorvers.

   
  -------------------------------------------------------------- */


#ifndef __CLTE_H__
#define __CLTE_H__

#include <vector>
#include "mtypes.h"
#include "model.h"
//#include "cprofiles2.h"
#include "ceos.h"
#include "cmemt2.h"

class clte{
 private:
  size_t nreg, nlin, ndep;
  std::vector<line> lin;
  std::vector<region> reg;
  std::vector<double>  scatt, sf;
  std::vector<size_t> roff;
  std::vector<float> part, frac;
  mat<double> opac;
 public:
  size_t nw;
  std::vector<double> lambda;

  clte(std::vector<region> &reg, std::vector<line> &lin);
  ~clte(){};
  
  double vac2air(double alamb);
  double air2vac(double alamb);
  void synth_cont(modl::mdepth &m, double *syn, ceos &eos, double mu = 1.0, int solver = 0);
  void delobez3_int(int ndep, double *z, double *op, double *sf, double &syn, double mu);
  double plank_nu(const double nu, const double temp);
};

#endif
