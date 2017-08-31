#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include "ceos.h"
#include "mtypes.h"
#include "cop.h"
//
using namespace std;


/* Some definitions */

const float ceos::AMASS[MAX_ELEM]=
  {1.008,  4.003,  6.941,  9.012, 10.811, 12.011, 14.007, 15.999,
   18.998, 20.179, 22.990, 24.305, 26.982, 28.086, 30.974, 32.060,
   35.453, 39.948, 39.102, 40.080, 44.956, 47.900, 50.941, 51.996,
   54.938, 55.847, 58.933, 58.710, 63.546, 65.370, 69.720, 72.590,
   74.922, 78.960, 79.904, 83.800, 85.468, 87.620, 88.906, 91.220,
   92.906, 95.940, 98.906,101.070,102.905,106.400,107.868,112.400,
   114.820,118.690,121.750,127.600,126.905,131.300,132.905,137.340,
   138.906,140.120,140.908,144.240,146.000,150.400,151.960,157.250,
   158.925,162.500,164.930,167.260,168.934,170.040,174.970,178.490,
   180.948,183.850,186.200,190.200,192.200,195.090,196.967,200.590,
   204.370,207.190,208.981,210.000,210.000,222.000,223.000,226.025,
   227.000,232.038,230.040,238.029,237.048,242.000,242.000,245.000,
   248.000,252.000,253.000};



const char ceos::ELEMEN[MAX_ELEM][3]=
  {"H ",  "He",  "Li",  "Be",  "B ",  "C ",  "N ",  "O ",  "F ",  "Ne",
   "Na",  "Mg",  "Al",  "Si",  "P ",  "S ",  "Cl",  "Ar",  "K ",  "Ca",
   "Sc",  "Ti",  "V ",  "Cr",  "Mn",  "Fe",  "Co",  "Ni",  "Cu",  "Zn",
   "Ga",  "Ge",  "As",  "Se",  "Br",  "Kr",  "Rb",  "Sr",  "Y ",  "Zr",
   "Nb",  "Mo",  "Tc",  "Ru",  "Rh",  "Pd",  "Ag",  "Cd",  "In",  "Sn",
   "Sb",  "Te",  "I ",  "Xe",  "Cs",  "Ba",  "La",  "Ce",  "Pr",  "Nd",
   "Pm",  "Sm",  "Eu",  "Gd",  "Tb",  "Dy",  "Ho",  "Er",  "Tm",  "Yb",
   "Lu",  "Hf",  "Ta",  "W ",  "Re",  "Os",  "Ir",  "Pt",  "Au",  "Hg",
   "Tl",  "Pb",  "Bi",  "Po",  "At",  "Rn",  "Fr",  "Ra",  "Ac",  "Th",
   "Pa",  "U ",  "Np",  "Pu",  "Am",  "Cm",  "Bk",  "Cs",  "Es"};



const float ceos::ABUND_default[MAX_ELEM] =
  {-0.04048,-1.07,-10.95,-10.89, -9.44, -3.48, -3.99, 
   -3.11, -7.48, -3.95,  -5.71, -4.46, -5.57, -4.49, -6.59, -4.83, 
   -6.54, -5.48, -6.82,  -5.68, -8.94, -7.05, -8.04, -6.37, -6.65, 
   -4.50, -7.12, -5.79,  -7.83, -7.44, -9.16, -8.63, -9.67, -8.69, 
   -9.41, -8.81, -9.44,  -9.14, -9.80, -9.54,-10.62,-10.12,-20.00, 
   -10.20,-10.92,-10.35, -11.10,-10.18,-10.58,-10.04,-11.04, -9.80, 
   -10.53, -9.81,-10.92,  -9.91,-10.82,-10.49,-11.33,-10.54,-20.00, 
   -11.04,-11.53,-10.92, -11.94,-10.94,-11.78,-11.11,-12.04,-10.96, 
   -11.28,-11.16,-11.91, -10.93,-11.77,-10.59,-10.69,-10.24,-11.03, 
   -10.95,-11.14,-10.19, -11.33,-20.00,-20.00,-20.00,-20.00,-20.00, 
   -20.00,-11.92,-20.00, -12.51,-20.00,-20.00,-20.00,-20.00,-20.00, 
   -20.00,-20.00};

/* End definitions */

string pS(string str, const size_t num, const char paddingChar = ' '){
  if(num > str.size())
    str.resize(num, paddingChar);
  return str;
}

void s2f(string &in, char *out, unsigned n){
  strcpy(out, in.c_str());

  unsigned l = in.length();
  //unsigned npad = n - l;
  for(unsigned ii = l;ii<(n);ii++) out[ii] = ' ';
  // out[n-1] =0; 
  //out[n-1] = '\0';
}


void ceos::initEOS(vector<line_t> &lines){

  const string inam = "ceos::initEOS: ";
  
  //
  // Init list of background absorvers (this could be more elegant...)
  //
  species = {"H",   "H",     "H",     "He",   "He",   "He", 
	     "C",    "Al",   "Si",    "Si",    "Ca",   "Ca", 
	     "Mg",   "Mg",   "N",     "Fe",   "O"};


  ion =
    { 1, 2, 0, 1, 2, 3,
     1, 1, 1, 2, 1, 2,
     1, 2, 1, 1, 1};


  
  
  //
  // Add external line list
  //
  if(lines.size() > 0){
    int kk = (int)(species.size() + lines.size());
    ion.resize(kk);
    
    kk = (int)species.size();
    for(auto &it: lines){
      species.push_back(string(it.elem));
      ion[kk] = it.ion;
      it.off = kk++; 
    }
  }
  
  char cspecies[species.size()][8];
  for(unsigned ii=0;ii<species.size();ii++) {
    s2f(species[ii], cspecies[ii], 8);
  }
  
  
  

  // WARNING, these indexes are passed to fortran routines, where
  // arrays start at index 1. If used inside this wrapper, this
  // should be accounted for!
  //
  IXH1  = 1,  IXH2  = 2, IXHMIN = 3, IXHE1 = 4,  IXHE2 = 5,  IXHE3 = 6, IXC1 = 7;
  IXAL1 = 8,  IXSI1 = 9, IXSI2  = 10, IXCA1 = 11, IXCA2 = 12, IXMG1 = 13;
  IXMG2 = 14, IXN1 = 15, IXFE1 = 16, IXO1 = 17;

  
  //
  // Count elements on the list
  //
  NELEM = MAX_ELEM;
  NLIST = 0; // Load the default list
  NLINES = species.size();
  
  switch(eqcount_(ELEMEN, cspecies, &ion[0], NLINES, NLIST, NELEM))
    {
    case 0: break;
    case 1:
      cerr << inam <<"ERROR, found ilegal species name"<<endl;
      exit(0);
      break;
    case 2:
      cerr << inam << "ERROR, SPLSIZ must be larger than "<<species.size()<<endl;
      exit(0);
      break;
    }



  
  /*
    Merge default list and user provided list
  */
  
  idxspec.resize(NLINES, 0);
  ntotallist = 8*NLIST;
  totallist.resize(ntotallist);
  ntotallist /= 8;
  NLIST = 0;

   fract.resize(ntotallist);
      pf.resize(ntotallist);
  potion.resize(ntotallist);
  xamass.resize(ntotallist);

  switch(eqlist_(&abund[0], ELEMEN, &cspecies[0], &ion[0], &idxspec[0], &totallist[0],
		 NLINES, NLIST, ntotallist, NELEM))
    {
    case 0: break;
    case 1:
      cerr << inam <<"ERROR, found ilegal species name"<<endl;
      exit(0);
    case 2:
      cerr << inam << "ERROR, SPLSIZ must be larger than "<<species.size()<<endl;
      exit(0);
    case 3:
      cerr << inam << "ERROR, Missing ionization stage";
      exit(0);
    case 4:
      cerr << inam << "ERROR, e- is not the las item on the list";
      exit(0);
    case 5:
      cerr << inam << "ERROR, unreasonable abundances";
      exit(0);
    }
  // cout << inam << "EOS initialized with default list + background absorvers "<<endl;



}


ceos::ceos(double grav){

  string inam = "ceos::ceos: ";
  gravity = pow(10.0, grav);


  
  /* --- Init abundances --- */

  vector<iabund> modABUND;
  initAbundances(modABUND);
  
  
  
  /* --- Init arrays for fortran routines --- */
  
  vector<line_t> lines;
  initEOS(lines);

  
  
  /* --- Get unique species in the user line-list --- */
  
  unique();
  
}



ceos::ceos(vector<line_t> &lines, double grav){

  string inam = "ceos::ceos: ";
  gravity = pow(10.0, grav);

  
  /* --- Init abundances --- */

  vector<iabund> modABUND;
  initAbundances(modABUND);

  
  /* --- Init arrays for fortran routines --- */
  
  initEOS(lines);


  /* --- Get unique species in the user line-list --- */
  
  unique();
  
  
}

ceos::ceos(vector<line_t> &lines, vector<iabund> &ab, double grav){

  string inam = "ceos::ceos: ";
  gravity = pow(10.0, grav);

  
  /* --- Init abundances --- */

  initAbundances(ab);

  
  /* --- Init arrays for fortran routines --- */
  
  initEOS(lines);


  /* --- Get unique species in the user line-list --- */
  
  unique();
  
  
}

ceos::ceos(vector<line_t> &lines, int n, float *iab, double grav){

  string inam = "ceos::ceos: ";
  gravity = pow(10.0, grav);

  
  /* --- Init abundances --- */

  
  for(int ii=0; ii<MAX_ELEM; ii++) abund[ii] = pow(10., ABUND_default[ii]);
  
  if(n > 0)
    for(int ii=0; ii<n; ii++) pow(10., iab[ii]);
  
  tabund = 0.0;
  for(int ii=0; ii<MAX_ELEM; ii++) tabund += abund[ii];
  
  wsum = 0;
  asum = 0;
  for(int ii=0; ii<MAX_ELEM; ii++){
    abund[ii] /= tabund;
    
    wsum += AMASS[ii] * abund[ii];
    asum += abund[ii];
  }
  
  tabund = asum;
  avmol = wsum / asum;

  
  
  /* --- Init arrays for fortran routines --- */
  
  initEOS(lines);


  /* --- Get unique species in the user line-list --- */
  
  unique();

  
  /* --- Init species table --- */
  
  init_Species_table(lines);
  
}


// ------------------------------------------------------------------------- 
// Init abundances
// ------------------------------------------------------------------------- 
void ceos::initAbundances(vector<iabund> &ab, bool verbose)
{

  /* --- Init default abundances --- */
  
  for(int ii=0; ii<MAX_ELEM; ii++) abund[ii] = pow(10., ABUND_default[ii]);
  
  
  /* --- replace abunds --- */
  
  for(auto &it: ab){
    
    for(int ii=0; ii<MAX_ELEM; ii++){
      if(!strcmp(ELEMEN[ii], it.elem)){
	abund[ii] = pow(10., it.abund);
	if(verbose) fprintf(stderr, "ceos::initAbundances: Changed [%s] -> %7.3f\n", it.elem, it.abund);
	break;
      }
    } // ii
  } // it


  double sum = 0.0;
  for(int ii = 0; ii<MAX_ELEM; ii++) sum += abund[ii];
  
  for(int ii = 0; ii<MAX_ELEM; ii++) abund[ii] /= sum;
  tabund = 1.0;

   
  //
  // Estimate average molecular weight, excluding electrons?
  //
  wsum = 0;
  asum = 0;
  for(int ii = 0; ii<MAX_ELEM; ii++){
    wsum += AMASS[ii] * abund[ii];
    asum += abund[ii];
  }

  avmol = wsum / asum;

  
}


// -------------------------------------------------------------------------
// Unique species in line list of the EOS
// -------------------------------------------------------------------------
void ceos::unique(void){

  /* --- Temporary list that will contain spec + ion --- */
  vector<string> bla  = species;
  for(int ii=0; ii<(int)species.size();ii++) bla[ii] += std::to_string(ion[ii]);

  /* --- idx & uspec are class variable. idx maps the unique values to 
     elements in the original array uspec[idx[ii]] = spec[ii] --- */
  idx.resize(bla.size());
  uspec.resize(0); // Unique elements will be push_back(ed)
  
  /* --- Search unique elements, quite an inefficient implementation ;-) ---*/
  for(int jj = 0; jj< (int)bla.size(); jj++){
    bool isin = false;

    /* --- if found in the unique list already, 
       just indicate in which element --- */
    for(int ii=0; ii<(int)uspec.size();ii++){
      if(bla[jj] == uspec[ii]) {
	isin = true;
	idx[jj] = ii;
	continue;
      }
    }
    /* --- If not in the unique list, include it and 
       indicate in which element --- */
    if(!isin) {
      uspec.push_back(bla[jj]);
      idx[jj] = uspec.size()-1;

      /* --- idx1 is a class variable and it maps the position of 
	 the unique elements in one of the elements of the original array, 
	 needed to extract ONCE the species from the EOS --- */
      idx1.push_back(jj);
    }
  }
}

/* ---------------------------------------------------------------------------- */

void ceos::store_partial_pressures(int ndep, int k, float na, float ne){

  int nuspec = (int)uspec.size();
  
  /* --- always init the buffer if (k == 0) --- */
  if(k == 0){
    buf.rinit({2, nuspec + 1, ndep});
  }
  
  /* --- Copy partition function and partial pressure 
     (which is divided by the PF) 
     --- */
  for(int ii = 0; ii<nuspec; ii++){
    buf(k,ii,0) = fract[idxspec[idx1[ii]]-1];
    buf(k,ii,1) = pf[idxspec[idx1[ii]]-1];
  }

  /* --- Also store gass particle density and electron density --- */
  buf(k,nuspec,0) = na;
  buf(k,nuspec,1) = ne;
}

/* ---------------------------------------------------------------------------- */

void ceos::read_partial_pressures(int k, std::vector<float> &frac, std::vector<float> &part, float &xa, float &xe){

  int ndep = (int)buf.shape(0);
  
  if(k >= ndep){
    cerr << "ceos::read_partial_pressures: ERROR, ndep["<<ndep<<"] <= k["<<k<<"]"<<endl;
    return;
  }

  int nspec = (int)idx.size();
  int nuspec = (int)buf.shape(1);

  
  frac.resize(nspec);
  part.resize(nspec);


  for(int ii=0;ii<nspec;ii++){
    frac[ii] = buf(k,idx[ii], 0);
    part[ii] = buf(k,idx[ii], 1);
  }

  xa = buf(k,nuspec-1, 0);
  xe = buf(k,nuspec-1, 1);

}

/* ---------------------------------------------------------------------------- */

double ceos::Pe_from_Pg(double T, double Pg, double *dummy){

  //
  // Estimate Pelect from Pgas, assuming ionization fraction
  //
  float iPe = init_pe_from_T_pg((float)T, (float)Pg);
  

  //
  // Call EQSTAT
  //
  int niter;
  int dum = MAX_ELEM;
  int mode = 0;
  float iT = (float)T, iPg = (float)Pg;


  eqstat_(mode, iT, iPg, iPe, &abund[0], ELEMEN, &AMASS[0], dum, &idxspec[0],
	  &totallist[0], &fract[0], &pf[0], &potion[0], &xamass[0],
	  NLINES, NLIST, xne, xna, RHOest, niter);
  
  return (double)xne * bk * T;
}

/* ---------------------------------------------------------------------------- */

double ceos::Pg_from_Rho(double iT, double irho, double &Pe){

  //
  // Estimate Pelect from Pgas, assuming ionization fraction
  //

  float T = (float)iT, rho = (float)irho;
  float Pg = rho * bk * T / (avmol * mp); // Init approximate gas pressure
  float iPe = (double)init_pe_from_T_pg(T, Pg);

  //
  // Call EQSTAT_RHO
  //
  int niter;
  int dum = MAX_ELEM;
  int mode = 0;

  eqstat_rho_(mode, T, Pg, iPe, &abund[0], ELEMEN, &AMASS[0],
	      dum, &idxspec[0], &totallist[0], &fract[0], &pf[0], &potion[0], &xamass[0],
	      NLINES, NLIST, xne, xna, rho, niter);
  
  Pe = (double)xne * iT * bk;
  return (double)Pg;
}

/* ---------------------------------------------------------------------------- */

float ceos::init_pe_from_T_pg(float t, float pg)
{
  /* --- Init Pe assuming everything is Hydrogen --- */
      	
  float nu=0.9091;//       ! assume that only Hydrogen is ionized
  float saha=pow(10.0,-0.4771+2.5*log10(t)-log10(pg)-(13.6*5040./t));
  float aaa=1.0+saha;
  float bbb=-(nu-1.)*saha;
  float ccc=-saha*nu;
  float ybh=(-bbb+sqrt(bbb*bbb-4.*aaa*ccc))/(2.*aaa); //! ionization fraction
  return (float)(pg*ybh/(1.+ybh));
}

/* ---------------------------------------------------------------------------- */

double ceos::Pg_from_Pe(double T, double Pe, double *dumm){
  
  //
  // Estimate Pelect from Pgas, assuming ionization fraction
  //
  float nne = Pe / (bk * T);
  float Pg;
  // float Pg = rho * bk * T / (avmol * mp); // Init approximate gas pressure
  
  if(T > 8000) Pg = 1.0/0.5;
  else if(T > 4000) Pg = 1.0/0.1;
  else if(T > 2000) Pg = 1.0/0.01;
  else Pg = 1.0/0.001;
  
  Pg *= Pe;
  xne = Pe / (bk * T);
  
  //
  // Call EQSTAT and iterate to get consistent rho with Pg
  // (should be changed by EQSTAT_RHO)!!
  //
  int niter;
  int dum = MAX_ELEM;
  int mode = 10;
  int dir = 0;
  float scale = 2.0, dif = 1.e5, rho_est=0.0;
  int myiter = 0;
  
  //
  float iT = T;
  float iPe = (float)Pe;
  
  eqstat_(mode, iT, Pg, iPe, &abund[0], ELEMEN, &AMASS[0], dum, &idxspec[0],
	  &totallist[0], &fract[0], &pf[0], &potion[0], &xamass[0],
	  NLINES, NLIST, xne, xna, rho_est, niter);
  

  float tol =1.e-4;
  
  while(fabs(dif) > tol){
    myiter++;
    
    
    if(dif > tol){
      if(dir != 1) scale = sqrt(scale);
      Pg *= scale;
      dir = 1;
    }else if(-dif > tol){
      if(dir != -1) scale = sqrt(scale);
      Pg /= scale;
      dir = -1;
    }
    
    eqstat_(mode, iT, Pg, iPe, &abund[0], ELEMEN, &AMASS[0], dum, &idxspec[0],
	    &totallist[0], &fract[0], &pf[0], &potion[0], &xamass[0],
	    NLINES, NLIST, xne, xna, rho_est, niter);
    
    dif = ((nne - xne)/nne);
    
  }
  
  
  return (double)Pg;
}
