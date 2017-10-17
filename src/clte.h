/* --------------------------------------------------------------

   LTE class, only continuum opacities for the time being
   Written by J. de la Cruz Rodriguez (ISP-SU 2016)

   NOTES: 
   Unpolarized Bezier-3 and linear solvers for the time being.


   Requires ceoswrap class to retrieve the partial densities of 
   background absorvers.

   
  -------------------------------------------------------------- */

#ifndef __CLTE_H__
#define __CLTE_H__

#include <vector>
#include "mtypes.h"
#include "model.h"
#include "ceos.h"


class clte{
 private:
  size_t nreg, nlin, ndep;
  line *lin;
  std::vector<region> reg;
  std::vector<size_t> roff;
  bool lines;
 public:
  size_t nw;
  std::vector<double> lambda;
  std::vector<float>  tau_eq_1_z;
  std::vector<double> C;
  
  
  clte(std::vector<region> &reg, std::vector<line> &lin);
 ~clte(){};
  
  double vac2air(double alamb);
  double air2vac(double alamb);
  void synth_nonpol(modl::mdepth &m, double *syn, eoswrap &eos, double mu, int solver, bool getContrib);

  void delobez3_int(int ndep, double *z, double *op, double *sf, double &syn, double mu, float &tau_eq_1);
  void linear_int(int ndep, double *z, double *op, double *sf, double &syn,double mu, float &tau_eq_1);
  double lte_opac(double temp, double n_u, double gf, double elow, double nu0);
};

#endif
