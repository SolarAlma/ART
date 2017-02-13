#ifndef __MDEP_H__
#define __MDEP_H__


#include "cmemt2.h"
#include "ceos.h"

namespace modl{

  enum munit{
    cgs = 0,
    si = 1
  };
  
  
  class mdepth{
  private:
    bool allocated;
  public:
    size_t ndep;
    munit units;	
    mat<double> buf;
    double *temp, *vz, *vx, *vy, *B, *vturb, *inc,
      *azi, *pgas, *rho, *nne, *z, *ltau;
    

    /* ---- Prototypes ---- */

    mdepth();
    mdepth(int ndep_in, bool alloc = true);
    ~mdepth();

    void init(int ndep_in, bool alloc = true, modl::munit unit_in = cgs);
    void to_cgs();
    void to_si();
    void fillDensities(info &inp, ceos &eos);
  };
  
};

#endif
