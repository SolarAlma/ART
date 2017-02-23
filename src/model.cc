#include "model.h"
#include "cmemt2.h"
#include "ceos.h"
#include <cmath>
#include "physical_consts.h"

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

void modl::mdepth::fillDensities(info &inp, ceos &eos){

  
  /* --- Some checks --- */

  double otau=0, tau=0, kappa, kappa_old, na, ne;
  vector<float> frac;
  double lambda = 5000.0, scatt = 0;
  float xna;
  
  for(size_t kk=0; kk<ndep; kk++){
    
    /* --- Fill densities using the EOS --- */
    
    if(inp.vdef[10]){
      if(!inp.vdef[12]) nne[kk] = eos.nne_from_T_Pg(temp[kk], pgas[kk], rho[kk]);
      else                        eos.nne_from_T_Pg_nne(temp[kk], pgas[kk], pgas[kk], nne[kk]);
      
    }else if(inp.vdef[11]){ 
      if(!inp.vdef[12])	nne[kk] = eos.nne_from_T_rho(temp[kk], pgas[kk], rho[kk]);
      else                        eos.nne_from_T_rho_nne(temp[kk], pgas[kk], rho[kk], nne[kk]);
    }else               rho[kk] = eos.rho_from_T_nne(temp[kk], pgas[kk], nne[kk]);

    
    /* --- Store partial pressures and partition densities to 
       compute background opacities later --- */
    
    xna =  pgas[kk] / (temp[kk] * phyc::BK);
    eos.store_partial_pressures((int)ndep, (int)kk, xna, (float)nne[kk]);

    /* --- Compute z-scale if ltau500 is the input depth scale --- */
    
    if(inp.vdef[0] && !inp.vdef[1]){
      tau  = pow(10.0, ltau[kk]);
      na = pgas[kk] / (temp[kk] * phyc::BK);
      ne = nne[kk];
      
      eos.contOpacity(temp[kk], 1,  &lambda, &kappa, &scatt, eos.fract, na, ne);
      
      if(kk == 0) {
	z[0] = 0;
      }else{
	z[kk] = z[kk-1] - 2.0 * (tau-otau) / (kappa-kappa_old);
      }
      otau = tau;
      kappa_old = kappa;
    }

    /* --- Compute ltau500 if Z is the input depth scale --- */
    
    if(inp.vdef[1] && !inp.vdef[0]){
      na = pgas[kk] / (temp[kk] * phyc::BK);
      ne = nne[kk];
      eos.contOpacity(temp[kk], 1,  &lambda, &kappa, &scatt, eos.fract, na, ne);
      
      if(kk == 0){
	tau  = 0.0;//;0.5 * kappa * (z[0] - z[1]);
	ltau[kk] = log10(tau);
      }else{
	tau = otau + 0.5 * (kappa_old + kappa) * (z[kk-1] - z[kk]);
        ltau[kk] = log10(tau);
      }
      
      otau = tau;
      kappa_old = kappa;
    }
    
  }// kk

  
  /* --- check ltau at the upper layer --- */
  
  tau = pow(10.0, ltau[1]), otau = pow(10.0, ltau[2]);
  if(inp.vdef[1] && !inp.vdef[0]) ltau[0] = log10(exp(2*log(tau)-log(otau)));
	
  
}


/* ---------------------------------------------------------------------------- */
