/*
  CEOS class
  Author: Jaime de la Cruz Rodriguez (ISP-SU 2014)
  Dependencies: Fortran EOS (N. Piskunov, UU)
 */
#ifndef CEOS_H
#define CEOS_H
//
#include <iostream>
#include <vector>
#include <string>
#include "mtypes.h"
#include "cmemt2.h"
#include "eoswrap.h" 


/* Interface for Fortran routines */
extern "C" {
  int eqcount_(const char [][3], const char [][8], int *, int &, int &, int &);

  int eqlist_(const float *, const char [][3], const char [][8], int *, int *, char *, int &,
	      int &, int &,  int&);

  void eqstat_(int &, float &, float &, float &, float *, const char [][3],
	       const float *, int &, int *, char *, float *, float *, float *,
	       float *, int &, int &, float &, float &, float &, int &);
  void eqstat_rho_(int &, float &, float &, float &, float *, const char [][3],
	       const float *, int &, int *, char *, float *, float *, float *,
	       float *, int &, int &, float &, float &, float &, int &);
  //     subroutine eqstat_rho(mode,temp,Pg,Pe,abund,elemen,amass,
  //   &                  ELESIZ,spindx,splist,xfract,poti,atwght,
  //   &                  nlines,nlist,xne,xna,rho,niter)
  //     subroutine     eqstat(mode,temp,Pg,Pe,abund,elemen,amass,
  //   &                  ELESIZ,spindx,splist,xfract,pfunc,poti,atwght,
  //   &                  nlines,nlist,xne,xna,rho,niter)
  
  void check_(const char [][3], char [][8], int *, int &, int &, int &);

  void xsaha_(int &, float &, float &, float &, int &, float *, double *, int &);

  /*
  void contop_(float &T, float &TKEV, float &TK, float &HKT, float &TLOG,
	       float &XNA, float &XNE, double *WLGRID, double *OPACITY,
	       double *SCATTER,  
	       float &H1, float &H2, float &HMIN, float &HE1, float &HE2,
	       float &HE3, float &C1, float &AL1, float &SI1, float &SI2,
	       float &CA1, float &CA2, float &MG1, float &MG2, float &FE1,
	       float &N1, float &O1, int &nWLGRID, int &NLINES, int &NTOTALLIST);
  */
}


struct iabund{
  char elem[3];
  float abund;
};


/* class definition */
class ceos: public eoswrap{
 private:
  
  
 public:
  // Atomic data (see implementation file)
  // static const int   MAX_ELEM = 99;

  static const float AMASS[MAX_ELEM];
  static const char  ELEMEN[MAX_ELEM][3];
  static const float ABUND_default[MAX_ELEM];
  //float ABUND[MAX_ELEM];
  int NELEM, NLIST, NLINES;

  
  // Physical constants
  static constexpr double bk = 1.3806488E-16;
  static constexpr double cc = 2.99792458E10;
  static constexpr double mp = 1.672621777E-24;
  static constexpr double me = 9.10938215E-28;

  // EOS book keeping
  double avmol, wsum, asum, gravity;
  int    IXH1,IXH2,IXHMIN,IXHE1,IXHE2,IXHE3,IXC1,IXAL1,IXSI1, 
         IXSI2,IXCA1,IXCA2,IXMG1,IXMG2,IXFE1,IXN1,IXO1;

  std::vector<float> fract, pf, potion, xamass;
  std::vector<std::string> species;
  std::vector<int>   idxspec, ion;
  std::vector<int> idx, idx1;
  std::vector<std::string> uspec;

  float xne, xna, RHOest;
  std::vector<char> totallist;
  int ntotallist;

  mat<float> buf;
  
  // Functions
  
  ceos(double grav = 4.44);
  ceos(std::vector<line_t> &lines, double grav = 4.44);
  ceos(std::vector<line_t> &lines, std::vector<iabund> &ab, double grav = 4.44);
  ceos(std::vector<line_t> &lines, int n = 0, float *iab = NULL, double grav = 4.44);

  void initAbundances(std::vector<iabund> &ab, bool verbose = false);
  void initEOS(std::vector<line_t> &lines);
  
  
  ~ceos(){
    totallist.clear();
  }
  double Pe_from_Pg    (double T,  double Pg, double *dum = NULL);
  double Pg_from_Pe    (double T,  double Pg, double *dum = NULL);
  double Pg_from_Rho    (double T,  double rho, double &Pe);


  void store_partial_pressures(int ndep, int k, float na, float ne);
  void read_partial_pressures(int k, std::vector<float> &frac, std::vector<float> &part, float &xa, float &xe);
  void unique(void);

  float init_pe_from_T_pg(float t, float pg);
};


#endif
