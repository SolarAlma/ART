/* -----------------------------------------------------------
   
   Base class for the equation of state. 
   It allows to use Wittmann-Mihalas or Piskunov's EOS cleanly.
   
   Written by J. de la Cruz Rodriguez (ISP-SU, 2017)
   
   ----------------------------------------------------------*/

#include <cmath>
#include <vector>
#include "mtypes.h"
#include "eoswrap.h"
#include "physical_consts.h"
#include "partition.h"
#include "cop.h"

using namespace std;


inline double isaha(const double t, const double xne, const double u0, const double u1, const double eion)
{
  static const double sfac = 2.0 * pow( (2.0 * phyc::PI * phyc::ME * phyc::BK) / (phyc::HH*phyc::HH), 1.5);
  return sfac * u1 / u0 * sqrt(t) * t * exp(-eion / (t*phyc::BK)) / xne;
}

/* --------------------------------------------------------------------------------------- */


void eoswrap::init_Species_table(std::vector<line> &lines) {

  /* --- Init background absorvers --- */
  
  spectab.resize(1), spectab[0].resize(10);
  
  spectab[0][0].anum = 0,  spectab[0][1].anum = 1, spectab[0][2].anum = 5,
  spectab[0][3].anum = 12, spectab[0][4].anum = 13,
  spectab[0][5].anum = 19, spectab[0][6].anum = 11, spectab[0][7].anum = 25,
  spectab[0][8].anum = 6, spectab[0][9].anum = 7;
  
  /* --- add elements from lines --- */

  size_t nlines = lines.size();
  for(size_t ii = 0; ii<nlines; ii++){

    size_t nn = spectab[0].size();
    size_t idx = 0;
    int anum = lines[ii].anum;
    bool exists = false;
    for(size_t jj=0; jj<nn; jj++){

      /* --- Check if we already have that atom --- */
      
      if(anum == spectab[0][jj].anum){
	exists = true;
	idx = jj;
	lines[ii].idx = (int)jj;
	break;
      }	
    } // jj

    /* --- If we have not allocated that element, just add it to the end of the array --- */
    
    if(!exists){
      species tmp = {};
      tmp.anum = anum;
      spectab[0].push_back(tmp);
      lines[ii].idx = nn; // no need to add +1 because the index number is nn-1
    }
    
  } // ii


  nspec = spectab[0].size(), nz=1;
  
}


/* --------------------------------------------------------------------------------------- */


void eoswrap::fill_Species_table(int kk, int ndep, double temp, double Pg, double Pe)
{

  /* --- need to resize the array ? --- */
  
  if(nz != ndep){
    spectab.resize(ndep);
    for(size_t ii = 1; ii<ndep; ii++) spectab[ii] = spectab[0];
    nz = ndep;
  }


  /* --- Fill species --- */

  
  double tkb = std::max(temp * phyc::BK,1.0e-22), xna = (Pg-Pe)/tkb, xne = Pe/tkb;
  double tmp[6], tmp1[6];
  
  for(size_t ii = 0; ii<nspec; ii++){

    species &ispec = spectab[kk][ii];
    double *xpa = &ispec.n_tot[0];
    
    memset(xpa, 0, 6*sizeof(double));
    memset(tmp,  0, 6*sizeof(double));
    memset(tmp1, 0, 6*sizeof(double));

    /* --- Get partition functions and scaled ionization potential for this atom --- */
    
    int nLev = pfn::partition_f<double>(ispec.anum, temp, xne, xna, tmp, tmp1, true);


    /* --- If hydrogen, append H- at the beginning --- */

    if(ispec.anum == 0){ // add H-
      for(int ii = 2; ii>0; ii--){
	tmp[ii] = tmp[ii-1];
	tmp1[ii] = tmp1[ii-1];
      }

      tmp[0] = 0.7542, tmp1[0] = 1.0;
      nLev = 3;
    }

    ispec.nLev = nLev;

    
    memcpy(ispec.pf,  tmp1, 6*sizeof(double));
    memcpy(ispec.eion, tmp, 6*sizeof(double));


    /* --- Solve partial densities for each ionization stage --- */

    double ntot = xna * abund[ispec.anum] / tabund;
    xpa[0] = 1.0;
    
    for(int ii=1; ii<nLev; ii++) xpa[ii] = isaha(temp, xne, tmp1[ii-1], tmp1[ii], tmp[ii-1]*phyc::EV);
    for(int ii=nLev-1;ii>0;ii--) xpa[0] = 1.0 + xpa[ii]*xpa[0];
    xpa[0] = 1.0/xpa[0];

    for(int ii=1; ii<nLev;ii++) xpa[ii] *= xpa[ii-1]; 
    for(int ii=0; ii<nLev;ii++) xpa[ii] *= ntot / ispec.pf[ii]; // divide by the partition function

    
    /* --- If hydrogen, place H- at the end --- */
    
    if(ispec.anum == 0){
      double dum = tmp[0], dum1 = tmp1[0], dum2 = xpa[0];
      for(int ii=0; ii<2; ii++){
	tmp[ii] = tmp[ii+1];
	tmp1[ii] = tmp1[ii+1];
	xpa[ii] = xpa[ii+1];
      }
      
      tmp[2] = dum, tmp1[2] = dum1, xpa[2] = dum2;
      
      memcpy(ispec.pf,  tmp1, 6*sizeof(double));
      memcpy(ispec.eion, tmp, 6*sizeof(double));

    }  
  }
}

/* --------------------------------------------------------------------------------------- */

void eoswrap::contOpacity(double T, double xne, int kk, int nw, double *w, double *opac)
{

  /* --- Just a wrapper for a call to cop --- */

  double TKEV = 8.6171E-5*T;
  double TK   = 1.38065E-16*T;
  double HTK  = 6.6256E-27/TK;
  double TLOG = log(T);
  int nlines = 0, ntotallist = 0;

  
  cop(T, TKEV, TK, HTK, TLOG, 0.0, xne, w, opac, NULL,
      spectab[kk][0].n_tot[0], // H (neutral)
      spectab[kk][0].n_tot[1], // H+
      spectab[kk][0].n_tot[2], // H-
      spectab[kk][1].n_tot[0], // He
      spectab[kk][1].n_tot[1], // He+
      spectab[kk][1].n_tot[2], // He++
      spectab[kk][2].n_tot[0], // C
      spectab[kk][3].n_tot[0], // Al
      spectab[kk][4].n_tot[0], // Si
      spectab[kk][4].n_tot[1], // Si+
      spectab[kk][5].n_tot[0], // Ca
      spectab[kk][5].n_tot[1], // Ca+
      spectab[kk][6].n_tot[0], // Mg
      spectab[kk][6].n_tot[1], // Mg+
      spectab[kk][7].n_tot[0], // Fe
      spectab[kk][8].n_tot[0], // N
      spectab[kk][9].n_tot[0], // O
      nw, nlines, ntotallist);
}
