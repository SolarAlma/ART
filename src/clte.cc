/* --------------------------------------------------------------

   LTE class, only continuum opacities for the time being
   Written by J. de la Cruz Rodriguez (ISP-SU 2016)

   NOTES: 
   Unpolarized Bezier-3 solver for the time being.
   
   Requires ceos class to retrieve the partial densities of 
   background absorvers.

   
  -------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include "clte.h"
#include "cprofiles2.h"
#include "physical_consts.h"
#include "mtypes.h"
#include "interpol.h"
#include "cmemt2.h"

using namespace std;
using namespace modl;

/* ---------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------- */

clte::clte(std::vector<region> &reg_in, std::vector<line> &lin_in)
{

  lin = lin_in, reg = reg_in;
  nw = 0, ndep = 0, nreg = reg.size(), nlin = lin.size();
  for(auto &it: reg) nw += it.nw;

  lambda.resize(nw);

  
  /* --- loop regions and collect wavs --- */

  size_t idx = 0, kk = 0;
  for(auto &it: reg){
    it.off = idx;
    idx += it.nw;
    it.wav.resize(it.nw);
    it.nu.resize(it.nw);

    for(size_t ii=0; ii<it.nw; ii++){
      it.wav[ii] = air2vac(it.w0 + it.dw*double(ii));
      it.nu[ii] = phyc::CC / (it.wav[ii] * 1.e-8);
      lambda[kk++] = vac2air(it.wav[ii]);
    }
    
  }

  /* --- Add 500 nm to the wav array to compute the reference opacity --- */
  
  lambda.push_back(5000.);
  scatt.resize(nw+1);
}

/* ---------------------------------------------------------------------------- */

double clte::vac2air(double alamb)
{
  if(alamb < 2000.) return alamb;
  else return alamb/(1.0+2.735182e-4+131.4182/alamb/alamb+ 
		     2.76249e8/alamb/alamb/alamb/alamb);
}

/* ---------------------------------------------------------------------------- */

double clte::air2vac(double lambda_air)
{
  // Init values
  double lambda_vacuum;
  if(lambda_air > 2000) lambda_vacuum = lambda_air / 1.00029;
  else lambda_vacuum = lambda_air;

  // Iterate
  double error = 1.0;
  while(error > 1.e-6){
    error = lambda_air - vac2air(lambda_vacuum);
    lambda_vacuum = lambda_vacuum + error / 1.0029;
  }
  return lambda_vacuum;
}

/* ---------------------------------------------------------------------------- */

void clte::delobez3_int(int ndep, double *z, double *op, double *sf, double &syn, double mu)
{
  
  double *dtau = new double [ndep](), itau = 0, *tau = new double [ndep](),
    dzu =  z[1]-z[0], deu =  (op[1] - op[0])/dzu,
    der = 0, oder = deu, dzd = 0, ded = 0;
  int k0 = 1, k1 = ndep-1;
  
  /* --- Integrate opacity to get tau scale --- */
  
  for(size_t k=1; k<(ndep-1); k++){
    size_t kd = k+1, ku = k-1;
    
    
    /* --- Compute side derivatives --- */

    dzd = z[kd] - z[k];
    ded = (op[kd] - op[k]) / dzd;
    
    
    /* --- Centered erivative of the opacity following Fritsch & Butland (1984) --- */
    
    if(deu*ded > 0.0){
      double lambda = (1.0 + dzd / (dzd + dzu)) / 3.0;
      der = (deu / (lambda * ded + (1.0 - lambda) * deu)) * ded;
    } else der = 0.0;
    
    /* --- 
       integrate opacity using cubic bezier splines
       The Bezier3 integral should be:
       dz * (op_0 + op_u + cntrl1 + cntr2) / 4 
       --- */
    
    dtau[k] = fabs(dzu) * ((op[k] - der/3.0 * dzu) +
			   (op[ku] + oder/3.0 * dzu) +
			   op[k] + op[ku]) * 0.25 * mu;
    itau += dtau[k];
    tau[k] = tau[k-1] + dtau[k];


    /* --- Store info --- */
    
    dzu  = dzd;
    deu  = ded;
    oder = der;
    
    
    if(itau <= 1.E-4) k0 = k;
    if(itau <= 100.0) k1 = k;
  }
  
  /* --- did we reach the lower boundary ? Bezier2 wih a single control point --- */
  
  dtau[ndep-1] = fabs(dzu) * ( (op[ndep-2] + oder/3.0 * dzu) + op[ndep-1]
				 + op[ndep-2]) / 3.0 * mu;
  tau[ndep-1] = tau[ndep-1] + dtau[ndep-1];
  if(k1 == (ndep-2)) k1 = ndep-1;
  
  
  /* --- Init the intensity with the value of the source function at the lowest point --- */
  
  syn = sf[k1];
  double *dsf = new double [ndep];
  cent_der<double>(ndep, tau, sf, dsf);
  
  
  
  /* --- Integrate ray --- */
 
  for(size_t k = k1-1; k >= k0; k--){
    
    int ku = k + 1;

    /* --- Integration coeffs. and exponential --- */
    
    double dt = dtau[ku];
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt03 = dt / 3.0;
    double eps, alp, bet, gam, mu;
    //
    if(dt >= 1.e-2){
      //
      eps = exp(-dt);
      //
      alp = (-6.0 + 6.0 * dt - 3.0 * dt2 + dt3 + 6.0 * eps)        / dt3;
      bet = (6.0 + (-6.0 - dt * (6.0 + dt * (3.0 + dt))) * eps)    / dt3;
      gam = 3.0 * (6.0 + (-4.0 + dt)*dt - 2.0 * (3.0 + dt) * eps)  / dt3;
      mu  = 3.0 * ( eps * (6.0 + dt2 + 4.0 * dt) + 2.0 * dt - 6.0) / dt3;
    }else{ // Taylor expansion of the exponential
      //
      double dt4 = dt2 * dt2;
      eps = 1.0 - dt + 0.5 * dt2 - dt3 / 6.0 + dt4 / 24.0;
      //
      alp = 0.25 * dt - 0.05 * dt2 + dt3 / 120.0 - dt4 / 840.0;
      bet = 0.25 * dt - 0.20 * dt2 + dt3 / 12.0  - dt4 / 42.0; 
      gam = 0.25 * dt - 0.10 * dt2 + dt3 * 0.025 - dt4 / 210.0; 
      mu  = 0.25 * dt - 0.15 * dt2 + dt3 * 0.05  - dt4 / 84.0; 
    }        
    
    /* --- integrate source function --- */

    double c_u = sf[ku] - dt03 * dsf[ku];
    double c_0 = sf[k]  + dt03 * dsf[k];
    
    syn = syn * eps + sf[k] * alp + sf[ku] * bet + c_0 * gam + c_u * mu;
  }
  
  delete [] dsf;
  delete [] tau;
  delete [] dtau;
  
 }

/* ---------------------------------------------------------------------------- */    

double clte::plank_nu(const double nu, const double temp){
  
    double c1 = (2.0 * phyc::HH * nu*nu*nu) / (phyc::CC*phyc::CC) ;
    double x = phyc::HH * nu / (phyc::BK * temp);

    if(x < 80.0) return c1 / (exp(x) - 1.0);
    else return         c1  * exp(-x);
  }

/* ---------------------------------------------------------------------------- */    

void clte::synth_cont(mdepth &m, double *syn, ceos &eos, double mu, int solver)
{

  /* --- Init vars --- */
  
  ndep = m.ndep;
  memset(&scatt[0], 0, (nw+1)*sizeof(double));
  double *psi = new double [m.ndep];
  
  if(sf.size() != m.ndep){
    sf.resize(m.ndep);
    opac.rinit(nw+1, m.ndep);
  }
  opac.zero();
  
  
  /* --- Loop depth --- */

  float xna=0, xne=0;
  for(size_t k = 0; k<ndep; k++){
    
    /* --- read partial pressures --- */
    
    eos.read_partial_pressures(k, frac, part, xna, xne);


    
    /* --- Compute cont. opacity for all wavelengths --- */

    eos.contOpacity(m.temp[k], nw+1,  &lambda[0], &opac(k,0), &scatt[0], frac, xna, xne);
	
  }// k


  /* --- Now integrate emerging intensity for each wavelength --- */

  for(size_t ir=0; ir<nreg; ir++){
    for(int w = 0; w< reg[ir].nw; w++){ // Loop lambda

      /* --- Source function at that wavlength --- */

      for(int k = 0; k<m.ndep; k++){
	sf[k] = plank_nu(reg[ir].nu[w], m.temp[k]);
	psi[k] = opac(k,w+reg[ir].off);
      }

    
      /* --- Compute formal solution at this wavelength, select method --- */
      // void delobez3_int(int ndep, double *z, double *op, double *sf, double &syn, double mu);

      delobez3_int(m.ndep, m.z, psi, &sf[0], syn[w + reg[ir].off], mu);
      
    }
  }

  delete [] psi;
  
}
