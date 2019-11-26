/* --------------------------------------------------------------

   LTE class, only continuum opacities for the time being
   Written by J. de la Cruz Rodriguez (ISP-SU 2016)

   NOTES: 
   Unpolarized Bezier-3 and linear solvers for the time being.


   Requires ceoswrap class to retrieve the partial densities of 
   background absorvers.

   
  -------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include "clte.h"
#include "cprofiles.h"
#include "physical_consts.h"
#include "mtypes.h"
#include "interpol.h"
#include "eoswrap.h"

using namespace std;
using namespace modl;
template<typename T> inline T sq(T const &var){return var*var;}



/* ---------------------------------------------------------------------------- */

clte::clte(std::vector<region> &reg_in, std::vector<line> &lin_in)
{

  lin = &lin_in[0], reg = reg_in;
  nw = 0, ndep = 0, nreg = reg.size(), nlin = lin_in.size();
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
      it.wav[ii] = air2vac(vac2air(it.w0) + it.dw*double(ii));
      it.nu[ii] = phyc::CC / (it.wav[ii] * 1.e-8);
      lambda[kk++] = vac2air(it.wav[ii]);
    }
    
  }

  /* --- allocate array to store z(tau=1) --- */

  tau_eq_1_z.resize(nw, 0.0);

  

  lines = ((nlin > 0) ? true : false);
  
  /* --- Init Zeeman components --- */

  if(lines)
    for(size_t ii=0;ii<(size_t)nlin; ii++)
      cprof::init_zeeman_components(lin_in[ii]);

} 

/* ---------------------------------------------------------------------------- */

double clte::vac2air(double alamb)
{
  if(alamb < 2000.) return alamb;
  else return alamb/(1.0+2.735182e-4+131.4182/alamb/alamb+ 
		     2.76249e8/alamb/alamb/alamb/alamb);
}

/* ---------------------------------------------------------------------------- */

double clte::air2vac(double lambda_air) {
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

// -------------------------------------------------------------------------
// LTE opacity, combination of Mihalas (1970), pag. 68 - Eq. 3.4 &
// Rutten (2003) eq. 2.98, pag. 31
// -------------------------------------------------------------------------

double clte::lte_opac(double temp, double n_u, double gf, double elow, double nu0){
  static const double lte_const = phyc::PI*phyc::EE*phyc::EE/(phyc::ME*phyc::CC);
  //
  double tk = phyc::BK * temp;
  return (lte_const * gf * n_u * exp(-elow / tk) * (1.0 - exp( -(phyc::HH * nu0) / tk)));
}

/* ---------------------------------------------------------------------------- */    

void clte::synth_nonpol(mdepth &m, double *syn, eoswrap &eos, double mu, int solver, bool getContrib) {

  /* --- Init vars --- */
  
  ndep = m.ndep;
  double *sf = new double [ndep], *opac = new double [ndep], *pC;
  int k0 = m.k0, k1 = m.k1;

  
  /* --- Check dimensions of array to store contribution functions --- */

  size_t depwavtot = size_t(nw) * size_t(ndep);
  
  if((bool)getContrib){
    if(C.size() != depwavtot) C.resize(depwavtot, 0.0);
    memset(&C[0], 0, depwavtot*sizeof(double));
  }
  
  
  /* --- Now integrate emerging intensity for each wavelength --- */

  for(size_t ir=0; ir<nreg; ir++){
    for(int w = 0; w < reg[ir].nw; w++){ // Loop lambda
      size_t idx = w + reg[ir].off;
      double ilambda = lambda[idx], inu =  reg[ir].nu[w];
      
      for(int kk = k0; kk<=k1; kk++){
        
        /* --- Get background opacity --- */

        eos.contOpacity(m.temp[kk], m.nne[kk], kk, 1, &ilambda, &opac[kk]);

        /* --- Planck function --- */

        sf[kk] = cprof::planck_nu<double>(inu, m.temp[kk]);

        /* --- Compute line opacities if needed --- */

        if (this->lines) {

        /*--- neutral Hydrogen and neutral Helium for Barklem's van der Waals broadening --- */

          double nH1 =  eos.spectab[kk][0].n_tot[0] * eos.spectab[kk][0].pf[0];
          double nHe1 = eos.spectab[kk][1].n_tot[0] * eos.spectab[kk][1].pf[0];

        /* --- Cycle each line and sum the opacity contribution --- */
	  
	        for(size_t ii = 0; ii<nlin; ii++){
            line &li = lin[ii];
	          if(fabs(reg[ir].wav[w] - li.w0) > li.width) continue;
 
	      /* ---  ionization stage population / partition function --- */
	    
	          double n_u = eos.spectab[kk][li.idx].n_tot[li.ion-1];
	    
	      /* --- Doppler width and damping --- */

	          double dlnu = li.nu0 * cprof::get_doppler_factor(m.temp[kk], m.vturb[kk], li.amass);
	          double damp = cprof::damp(li, m.temp[kk], m.nne[kk], nH1, nHe1, dlnu);
	    
	      /* --- LTE opac --- */
	    
	          double lineopac = lte_opac(m.temp[kk], n_u, li.gf, li.e_low, li.nu0);
	    
	      /* ---Normalized Voigt profile * opacity --- */
	    
	          opac[kk] +=  lineopac * cprof::getProfile(inu, li, m.vz[kk], dlnu, damp);
	        } // ii
	      } // if lines
      } // kk

      /* --- Compute formal solution at this wavelength, select method --- */

      if((bool)getContrib){
        pC = &C[idx*ndep+k0];
      } else {
        pC = NULL;
      }
      
      if(solver == 0) cprof::bezier3_int(k1-k0+1, &m.z[k0], &opac[k0], &sf[k0], syn[idx], mu, tau_eq_1_z[idx], pC);
      else             cprof::linear_int(k1-k0+1, &m.z[k0], &opac[k0], &sf[k0], syn[idx], mu, tau_eq_1_z[idx], pC);
    } // w
  } // ir
  
  /* --- Cleanup --- */

  delete [] sf;
  delete [] opac;
  
}



inline void integrate_alma_pol(int const ndep, int const nw, const double *z, const double *op, const double *temp, double *syn,
			double mu, double const iNu, const double *B, const double *inc, const double *Ne)
{
  
  constexpr static const double vE = phyc::EE / (2.0 * phyc::PI * phyc::ME * phyc::CC);
  static const double wP = phyc::EE / sqrt(phyc::PI * phyc::ME);
  static constexpr const double sig[2] = {-1.0, 1.0};
  
  
  
  for(int ii=0; ii<2; ++ii){
    double itau = 0.0;
    double const iSig = sig[ii];
    
    {
      double const dS = fabs(z[0]-z[1]) / mu;
      itau += dS * op[0];
    }

    syn[ii*nw] = 0.0;
    
    for(int kk=1; kk<ndep; ++kk){
      
      double const vB = vE * B[kk];
      double const v = sq(wP * sqrt(Ne[kk]) / iNu); 
      double const u = sq(vB / iNu);
      
      double const cosi = cos(inc[kk]);
      double const sini = sin(inc[kk]);
       
       
       double const D = sq(u) * sq(sq(sini)) + 4*u*sq(1.0-v)*sq(cosi);
       double const sD = sqrt(D);
       
       double const n2 = 1.0 - (2.0*v*(1.0-v)) / (2.0*(1-v) - u*sq(sini) + iSig*sD);
       if(itau > 15. || n2 < 0.0) break;

       
       double const Fs = 2.0 * (iSig*sD*(u*sq(sini)+2.0*sq(1.0-v))-sq(u*sq(sini))) /
	                                (iSig*sD*sq(2.0*(1.0-v)-u*sq(sini)+iSig*sD));
       
       double const opSig = op[kk] * (Fs / sqrt(fabs(n2)));
       
       double const dS = fabs(z[kk-1]-z[kk]) / mu;
       double const eps = exp(-itau);
       
       syn[ii*nw] += (1.0-exp(-opSig*dS))*temp[kk]*eps;
       itau += op[kk] * dS; // eps should not include the local contribution, which we add now for next height
       
       // cerr<<ii<<" "<<kk<<" "<<itau<<" "<<opSig<<" "<<dS*1.e-5<<endl;
     } // kk
   }// ii
   

   syn[2*nw] = (syn[0]-syn[nw]) / (syn[nw] + syn[0]);
 }


void clte::synth_pol_alma_only(modl::mdepth &m, double *syn, eoswrap &eos, double mu, int solver, bool getContrib)
{

  /* --- Init vars --- */
  
  ndep = m.ndep;
  double* __restrict__ opac = new double [ndep], *pC;
  int k0 = m.k0, k1 = m.k1;
  int  nw = 0;

  for(size_t ir=0; ir<nreg; ir++) nw+= int(reg[ir].nw);
  
  /* --- Check dimensions of array to store contribution functions --- */
  
  size_t depwavtot = size_t(nw) * size_t(ndep);
  
  if((bool)getContrib){
    if(C.size() != depwavtot) C.resize(depwavtot, 0.0);
    memset(&C[0], 0, depwavtot*sizeof(double));
  }
  
  
  /* --- Now integrate emerging intensity for each wavelength --- */
  
  for(size_t ir=0; ir<nreg; ir++){
    for(int w = 0; w < reg[ir].nw; w++){ // Loop lambda
      size_t idx = w + reg[ir].off;
      double ilambda = lambda[idx], inu =  reg[ir].nu[w];
      
      for(int kk = k0; kk<=k1; kk++){
	
        /* --- Get background opacity --- */
	
        eos.contOpacity(m.temp[kk], m.nne[kk], kk, 1, &ilambda, &opac[kk]);
	
      } // kk


      integrate_alma_pol(k1-k0+1, nw, &m.z[k0], &opac[k0], &m.temp[k0], &syn[idx], mu,  inu, &m.B[k0], &m.inc[k0], &m.nne[k0]);


    } // w
  } // ir

  delete [] opac;
  
}
