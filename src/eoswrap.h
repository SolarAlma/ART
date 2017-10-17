#ifndef EOS_WRAP_H
#define EOS_WRAP_H

#include <cmath>
#include <cstring>
#include <vector>
#include "mtypes.h"

/* --- 
   Wrapper class to both EOS included in the code
   Author: J. de la Cruz Rodriguez (ISP-SU 2017)
   --- */ 
   


class eoswrap{
 public:
  static const int   MAX_ELEM = 99;

  std::vector<std::vector<species>> spectab;
  size_t nspec, nz;
  float abund[MAX_ELEM];
  double tabund;
  
  eoswrap(){};
  eoswrap(std::vector<line> &lines, int n = 0, float *ab = NULL, double grav = 4.44){};
  ~eoswrap(){};

  void init_Species_table(std::vector<line> &lines);
  virtual void fill_Species_table(int kk, int ndep, double temp, double Pg, double Pe);
  virtual void contOpacity(double temp, double xne, int kk, int n, double *w, double *opac);

  virtual double Pe_from_Pg(double temp, double pg, double *fe_out = NULL) = 0;
  virtual double Pg_from_Pe(double temp, double pe, double * fe_out = NULL) = 0;
  virtual double Pg_from_Rho(double temp, double rho, double &pe) = 0;
};

#endif
