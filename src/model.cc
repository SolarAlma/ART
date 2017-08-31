#include "model.h"
#include "cmemt2.h"
#include <cmath>
#include "physical_consts.h"


#include "eoswrap.h"
#include "witt.h"
#include "ceos.h"

using namespace std;

/* ---------------------------------------------------------------------------- */

modl::mdepth::mdepth():allocated(false){
  ltau = NULL, z = NULL, temp = NULL, vz = NULL, vx = NULL, vy = NULL,
    vturb = NULL, pgas = NULL, rho = NULL, nne = NULL;
};

/* ---------------------------------------------------------------------------- */

modl::mdepth::mdepth(int ndep, bool alloc):allocated(false)
{
  ltau = NULL, z = NULL, temp = NULL, vz = NULL, vx = NULL, vy = NULL,
    vturb = NULL, pgas = NULL, rho = NULL, nne = NULL;
  init(ndep, alloc);
}

/* ---------------------------------------------------------------------------- */

modl::mdepth::~mdepth()
{   
  allocated = false;
  ltau = NULL, z = NULL, temp = NULL, vz = NULL, vx = NULL, vy = NULL,
    vturb = NULL, pgas = NULL, rho = NULL, nne = NULL;
}

/* ---------------------------------------------------------------------------- */

void modl::mdepth::init(int ndepin, bool alloc, modl::munit unit_in)
{
  ndep = ndepin;
  units = unit_in;
  
  if(alloc){
    buf.rinit(ndep, 13, true);
    allocated = true;

    /* --- Assign pointers --- */

    ltau = &buf(0,0);
    z    = &buf(1,0);
    temp = &buf(2,0);
    vz   = &buf(3,0);
    vx   = &buf(4,0);
    vy   = &buf(5,0);
    vturb= &buf(6,0);
    B    = &buf(7,0);
    inc  = &buf(8,0);
    azi  = &buf(9,0);
    pgas = &buf(10,0);
    rho  = &buf(11,0);
    nne  = &buf(12,0);
  } else allocated = false;
}

/* ---------------------------------------------------------------------------- */

void modl::mdepth::to_si()
{

  /* --- check if current units are cgs and convert variables to si --- */

  if(units == cgs){
    for(size_t zz = 0; zz<ndep; zz++){
      z[zz] *= 0.01; // CM_TO_M
      vz[zz] *= 0.01; // CM_TO_M
      vx[zz] *= 0.01; // CM_TO_M
      vy[zz] *= 0.01; // CM_TO_M
      vturb[zz] *= 0.01; // CM_TO_M
      nne[zz] *= 1.e6; // 1/ CM_TO_M**3
      B[zz] *= 1.e-4; // GAUSS_TO_TESLA
      pgas[zz] *= 0.1; //DYN/CM2_TO_PASCAL
      rho[zz] *= 0.001; // GR/CM3_TO_KG/M3
    }
    units = si;
  } else return; 
}

/* ---------------------------------------------------------------------------- */

void modl::mdepth::to_cgs()
{

  /* --- check if current units are si and convert variables to cgs --- */

  if(units == si){
    for(size_t zz = 0; zz<ndep; zz++){
      z[zz] *= 100.; // 1./CM_TO_M
      vz[zz] *= 100.; // 1./CM_TO_M
      vx[zz] *= 100.; // 1./CM_TO_M
      vy[zz] *= 100.; // 1./CM_TO_M
      vturb[zz] *= 100.; // 1./CM_TO_M
      nne[zz] *= 1.e-6; // CM_TO_M**3
      B[zz] *= 1.e4; // 1./GAUSS_TO_TESLA
      pgas[zz] *= 10.; // 1./ DYN/CM2_TO_PASCAL
      rho[zz] *= 1000.; // 1./(GR/CM3_TO_KG/M3)
    }
    units = cgs;
  } else return; 
}


/* ---------------------------------------------------------------------------- */

void modl::mdepth::fillDensities(info &inp, eoswrap &eos){

  
  /* --- Some checks --- */

  double otau=0, tau=0, kappa, kappa_old, na, ne, Pe;
  vector<float> frac;
  double lambda = 5000.0, scatt = 0;
  float xna;
  
  /* --- Check for T-cut conditions --- */
  
  getTcut(inp);

  
  /* --- Fill densities --- */
  
  for(size_t kk=k0; kk<=k1; kk++){

    /* --- Fill densities using the EOS --- */
    
    double tk = temp[kk] * phyc::BK;
    
    if(inp.vdef[10]){// Pgas 
      if(!inp.vdef[12]) nne[kk] = eos.Pe_from_Pg(temp[kk], pgas[kk], NULL) / tk;
      Pe = nne[kk] * tk; // nne is defined
    }else if(inp.vdef[11]){  // rho
      pgas[kk] = eos.Pg_from_Rho(temp[kk], rho[kk], Pe); // Pe is given as output
      if(!inp.vdef[12]) nne[kk] = Pe / tk;
      nne[kk] = Pe / tk;
    }else if(inp.vdef[12]){ // only nne given ?
      Pe = nne[kk] * tk;
      pgas[kk] = eos.Pg_from_Pe(temp[kk], Pe, NULL);
    }
    
    eos.fill_Species_table((int)kk, (int)ndep, temp[kk], pgas[kk], Pe);

    /* --- Do we need to compute the z-scale from ltau500? --- */
    
    // if(!inp.vdef[1] && inp.vdef[0]){
    //   na = (pgas[kk] - Pe) / tk, ne = nne[kk];
    //   if(kk == 0) z[0] = 0.0;
    //   else z[kk] = z[kk-1] - 2.0 * (tau-otau) / (kappa-kappa_old);
    // }
  }

  
}

/* ---------------------------------------------------------------------------- */

void modl::mdepth::getTcut(info &inp)
{
  k0 = 0, k1 = ndep-1;
  if(inp.temperature_cut < 0.0) return;
  
  for(int ii=0;ii<(ndep-2); ii++){
    if(temp[ii] >= inp.temperature_cut) k0 = ii;
    else break;
  }
}

/* ---------------------------------------------------------------------------- */
