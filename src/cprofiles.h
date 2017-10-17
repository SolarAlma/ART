/* ---------------------------------------------------

   Auxiliary routines to compute the line profile and opacity
   Formal solvers are also implemented here so we can re-use them
   in the future if needed.

   Damping includes radiative, Stark and vdW (with Barklem)

   J. de la Cruz Rodriguez (ISP-SU 2017)
   
   --------------------------------------------------- */

#ifndef CPROF_H
#define CPROF_H
//
#include <math.h>
#include <stdio.h>
#include <complex>
//
#include "interpol.h"
#include "physical_consts.h"

namespace cprof{  

  /* --- prototypes --- */
  
  template <class T> inline T planck_nu(const double nu, const T temp);
  template <class T> inline T get_doppler_factor(T temp, T vturb, T amass);
  template <class T> inline T damp(line_t &line, T temp, T vturb, T nne, T nh, T nhe, T dlnu);
  template <class T> inline double stark(double g_str, T temp, T nne);
  template <class T> inline double vanderWaals(line_t &line,  T temp, T nh,  T nhe);
  template <class T> inline T getProfile( double nu, const line_t &lin, T vel, T dlnu, T damping);
  template <class T> inline T voigt_complex(const T a, const T v, T &far);
  template <class T> inline T rvoigt(T a, T v);
  
  inline void getTauNu(int n, double *opac, double *z, double *dtau, int &k0, int &k1, float &z_tau, double mu);
  inline void init_zeeman_components(line_t &line);
  inline void linear_int(int ndep, double *z, double *op, double *sf, double &syn,
		  double mu, float &tau_eq_1);
  inline void bezier3_int(int ndep, double *z, double *op, double *sf, double &syn,
		    double mu, float &tau_eq_1);

  
  /* ---------------------------------------------------------------------------- */

  //----------------------------------------------------------------
  // Planck function at a given frequency for a given temperature
  //----------------------------------------------------------------
  template <class T> T planck_nu(const double nu, const T temp){
  
    double c1 = (2.0 * phyc::HH * nu*nu*nu) / (phyc::CC*phyc::CC) ;
    double x = phyc::HH * nu / (phyc::BK * temp);
  
    if(x < 80.0) return c1 / (exp(x) - 1.0);
    else return         c1  * exp(-x);
  }

  //----------------------------------------------------------------
  // Get Doppler factor, without the freq
  //----------------------------------------------------------------
  template <class T> T get_doppler_factor(T temp, T vturb, T amass){
    return (T)sqrt(2.0 * phyc::BK * temp / (amass * phyc::AMU)  + vturb * vturb) / phyc::CC;
  }

  //----------------------------------------------------------------
  // Damping, assuming VALD3 input
  //----------------------------------------------------------------
  template <class T> T damp(line_t &line, T temp, T nne, T nh, T nhe, T dlnu){
  
    /* --- radiative damping --- */
    double adamp = line.g_rad; // in input.cc, if g_rad is 0, then the value is filled with an approx. formula.
    adamp +=  vanderWaals<T>(line, temp, nh, nhe);  
    adamp += stark<T>(line.g_str, temp, nne);
    adamp /= 4.0 * phyc::PI* dlnu; // Norm from Rutten 2003, pag. 59

    
    return (T)adamp;
  }

  //----------------------------------------------------------------
  // Stark effect (collisions with charged particles),
  // Assuming VALD input, we get gamma_4 * nne (at 10000 K) to correct
  // for temperature, use x(temp/10000.)^(1./6.)
  //----------------------------------------------------------------
  template <class T> double stark(double g_str, T temp, T nne){
    return g_str * nne * pow(temp/10000.,0.16666666666666666); // Assuming Vald input!
  }

  //----------------------------------------------------------------
  // van der Waals broadening (collisions with neutral particles, typically H)
  // Includes Barklem formulation
  //----------------------------------------------------------------
  template <class T> double vanderWaals(line_t &line,  T temp, T nh,  T nhe){
    if(line.barklem){
    
      /* --- Barklem (constants initialized in input.cc) --- */
      //     http://www.astro.uu.se/~barklem/howto.html //
    
      if(line.firsttime){
	line.b_sig *= 2.80028e-21;
	double gx = (2.0 - line.b_alp*0.5) - 1.0;
      
	double gammaf =
	  1.0 + (-0.5748646 + (0.9512363 + (-0.6998588 + (0.4245549 - 0.1010678 * gx) * gx) * gx) * gx) * gx;
      
	line.b_gvw = pow((4.0 / phyc::PI), (line.b_alp * 0.5)) * gammaf * line.b_sig * 1.0e4;
	line.b_vbar =  21172.6 * ( 1.0 / 1.008 + 1.0 / line.amass);
	line.firsttime = false;
      }

      
      /* --- Actual calculation of vdW broadening --- */
    
      double vbar = sqrt(temp * line.b_vbar);

      return pow((vbar / 1.E4), (1.0 - line.b_alp)) *  line.b_gvw * (nh + 0.42*nhe) * 2.e6;
    
    }else if(line.g_vdw != 0.0){
    
      /* --- Unsoeld 1955 from Nikolai Piskunov's routines, assuming VALD input --- */
    
      return line.g_vdw * pow(temp / 10000., 0.3) * (nh + 0.42*nhe);
    
    }else return 0.0; 
  }
  
  /* ---------------------------------------------------------------------------- */
  
  template <class T> T getProfile( double nu, const line_t &lin, T vel, T dlnu, T damping)
  {

    double v  = (lin.nu0 - nu)  / dlnu;
    double va = lin.nu0 * vel   / (phyc::CC * dlnu);
    //double kk;

    double prof = rvoigt<double>(damping,  v - va) / (dlnu*phyc::SQPI);
    //double prof = voigt_complex<double>(damping,  v - va, kk) / (dlnu*phyc::SQPI);

    return (T)prof;
  }

  /* ---------------------------------------------------------------------------- */

  //-------------------------------------------------------------------------
  // Compute Zeeman splitting
  //-------------------------------------------------------------------------
  void init_zeeman_components(line_t &line)
  {
    /* --- 
       
       Compute the Zeeman splitting and strength for a line, 
       the Zeeman components are already scaled so the 
       sum{Red}=1, sum{Blue}=1 and sum{Par}=1
       
       --- */
    
    
    line.nZ = 0;
    int nup = int(2 * line.Jup) + 1; 
    int delta_j = (int)(line.Jup - line.Jlow);
    
    
    for(int iup = 1; iup <= nup; iup++){
      float Mup = line.Jup + 1 - iup;
      for(int ilow = 1; ilow <= 3; ilow++){
	float Mlow = Mup - 2 + ilow;
	if(fabs(Mlow) <= line.Jlow){
	  
	  /* --- Compute relative Zeeman strength, 
	     Landi Degl'innocenti & Landolfi (2004), 
	     table 3.1 - pag. 81 --- */
	
	  double strength = 0.0; // If abs(delta_j) > 1 then don't compute Zeeman splitting. 
	  //
	  if(delta_j == 1){
	  
	    if(ilow == 1)      strength = 1.5 * (line.Jlow + Mlow + 1.0) * (line.Jlow + Mlow + 2.0)
				 / ((line.Jlow+1.0)*(2.0*line.Jlow + 1.0) * (2.0 * line.Jlow + 3.0));
	    else if(ilow == 2) strength = 3.0 * (line.Jlow - Mlow + 1.0) * (line.Jlow + Mlow + 1.0)
				 / ((line.Jlow+1.0)*(2.0*line.Jlow + 1.0) * (2.0 * line.Jlow + 3.0));
	    else               strength = 1.5 * (line.Jlow - Mlow + 1.0) * (line.Jlow - Mlow + 2.0)
				 / ((line.Jlow+1.0)*(2.0*line.Jlow + 1.0) * (2.0 * line.Jlow + 3.0));
	  
	  } else if(delta_j == 0){
	  
	    if(ilow == 1)      strength = 1.5 * (line.Jlow - Mlow) * (line.Jlow + Mlow + 1.0)
				 / (line.Jlow * (line.Jlow + 1.0) * (2.0 * line.Jlow + 1.0));
	    else if(ilow == 2) strength = 3.0 * Mlow * Mlow
				 / (line.Jlow * (line.Jlow + 1.0) * (2.0 * line.Jlow + 1.0));
	    else               strength = 1.5 * (line.Jlow + Mlow) * (line.Jlow - Mlow + 1.0)
				 / (line.Jlow * (line.Jlow + 1.0) * (2.0 * line.Jlow + 1.0));
	  
	  } else if(delta_j == -1){
	  
	    if(ilow == 1)      strength = 1.5 * (line.Jlow - Mlow) * (line.Jlow - Mlow - 1.0)
				 / (line.Jlow * (2.0 * line.Jlow - 1.0) * (2.0 * line.Jlow + 1.0));
	    else if(ilow == 2) strength = 3.0 * (line.Jlow - Mlow) * (line.Jlow + Mlow)
				 / (line.Jlow * (2.0 * line.Jlow - 1.0) * (2.0 * line.Jlow + 1.0));
	    else               strength = 1.5 * (line.Jlow + Mlow) * (line.Jlow + Mlow - 1.0)
				 / (line.Jlow * (2.0 * line.Jlow - 1.0) * (2.0 * line.Jlow + 1.0));
	  
	  }
	
	
	  /* --- Zeeman splitting and strength ---*/
	
	  double splitting = line.Gup*Mup - line.Glow*Mlow;
	  line.splitting.push_back(splitting);
	  line.strength.push_back(strength);
	  line.iL.push_back(ilow-1);
	  line.nZ++;
	}
      }
    }
  }
  
  /* ---------------------------------------------------------------------------- */

  //----------------------------------------------------------------
  // Computation of the Voigt profile
  //----------------------------------------------------------------

  template <class T> T voigt_complex(const T a, const T v, T &far){
    //T vgt;
    T sav = fabs(v) + a;
    std::complex<T> tav(a, -v);
    std::complex<T> uav = tav*tav;
    std::complex<T> w4;
  
    /* --- HUMLICEK'S APPROXIMATION --- */
    if(sav >= 15.0){
      w4 = tav * 0.5641896 / (0.5 + uav);
    } else if(sav >= 5.5){
      w4 = tav * (1.410474 + uav * 0.5641896) / (0.75+uav * (3.0 + uav));
    } else if(a >= (0.195 * fabs(v) - 0.176)){
      w4 = (16.4955 + tav * (20.20933 + tav * (11.96482 + tav * (3.778987 +  tav * 0.5642236)))) / (16.4955 + tav * (38.82363 + tav * (39.27121 + tav * (21.69274 + tav * (6.699398 + tav)))));
    } else{
      w4 = tav * (36183.31 - uav * (3321.9905 - uav * (1540.787 - uav * (219.0313 - uav * (35.76683 - uav * (1.320522 -  uav * 0.56419))))));
      
      std::complex<T> v4 = (32066.6 - uav * (24322.84 - uav * (9022.228 -   uav * (2186.181 - uav * (364.2191 - uav * (61.57037 - uav * (1.841439 - uav)))))));
      w4 = exp(uav) - w4 / v4;
    }
    
    /* ---  Note that FVGT below is in fact 2 * (Faradey-Voigt function) ---*/
    //vgt = w4.real();
    far = 0.5 * w4.imag();

    return (T)w4.real();
  }



  /* ---------------------------------------------------------------------------- */
 
  template <class T> T rvoigt(T a, T v)
  {
    
    /* --- Humlicek's approximation with real numbers, only Voigt --- */
    
    T x1, y1, x2, y2, ti, ui, tr, ur, uu, vv, xx, yy, sav;
    T a2, v2;

    tr = a;
    ti = -v;
    a2 = a*a, v2 = v*v;

    ur = a2 - v2;
    ui = -2.0 * a  * v;
    sav = fabs(v) + a;
    if (sav >= 15.) {
      ur += .5;
      xx = ((a2 > v2)? a2 : v2);      
      tr /= xx;
      ti /= xx;
      ur /= xx;
      ui /= xx;
      return (tr * ur + ti * ui) * .5641896 / (ur * ur + ui * ui);
    } else if (sav >= 5.5) {
      x1 = ur * .5641896 + 1.410474;
      y1 = ui * .5641896;
      xx = x1 * tr - y1 * ti;
      yy = x1 * ti + y1 * tr;
      x1 = ur + 3.;
      y1 = ui;
      uu = x1 * ur - y1 * ui + .75;
      vv = x1 * ui + y1 * ur;
      return (xx * uu + yy * vv) / (uu * uu + vv * vv);
    } else if (a >= fabs(v) * .195 - .176) {
      x1 = tr * .5642236 + 3.778987;
      y1 = ti * .5642236;
      x2 = x1 * tr - y1 * ti + 11.96482;
      y2 = x1 * ti + y1 * tr;
      x1 = x2 * tr - y2 * ti + 20.20933;
      y1 = x2 * ti + y2 * tr;
      xx = x1 * tr - y1 * ti + 16.4955;
      yy = x1 * ti + y1 * tr;
      x1 = tr + 6.699398;
      y1 = ti;
      x2 = x1 * tr - y1 * ti + 21.69274;
      y2 = x1 * ti + y1 * tr;
      x1 = x2 * tr - y2 * ti + 39.27121;
      y1 = x2 * ti + y2 * tr;
      x2 = x1 * tr - y1 * ti + 38.82363;
      y2 = x1 * ti + y1 * tr;
      uu = x2 * tr - y2 * ti + 16.4955;
      vv = x2 * ti + y2 * tr;
      return (xx * uu + yy * vv) / (uu * uu + vv * vv);
    } else {
      x1 = 1.320522 - ur * .56419;
      y1 = -ui * .56419;
      x2 = 35.76683 - (x1 * ur - y1 * ui);
      y2 = -(x1 * ui + y1 * ur);
      x1 = 219.0313 - (x2 * ur - y2 * ui);
      y1 = -(x2 * ui + y2 * ur);
      x2 = 1540.787 - (x1 * ur - y1 * ui);
      y2 = -(x1 * ui + y1 * ur);
      x1 = 3321.9905 - (x2 * ur - y2 * ui);
      y1 = -(x2 * ui + y2 * ur);
      x2 = 36183.31 - (x1 * ur - y1 * ui);
      y2 = -(x1 * ui + y1 * ur);
      xx = x2 * tr - y2 * ti;
      yy = x2 * ti + y2 * tr;
      x1 = 1.841439 - ur;
      y1 = -ui;
      x2 = 61.57037 - (x1 * ur - y1 * ui);
      y2 = -(x1 * ui + y1 * ur);
      x1 = 364.2191 - (x2 * ur - y2 * ui);
      y1 = -(x2 * ui + y2 * ur);
      x2 = 2186.181 - (x1 * ur - y1 * ui);
      y2 = -(x1 * ui + y1 * ur);
      x1 = 9022.228 - (x2 * ur - y2 * ui);
      y1 = -(x2 * ui + y2 * ur);
      x2 = 24322.84 - (x1 * ur - y1 * ui);
      y2 = -(x1 * ui + y1 * ur);
      uu = 32066.6 - (x2 * ur - y2 * ui);
      vv = -(x2 * ui + y2 * ur);
      return exp(ur) * cos(ui) - (xx * uu + yy * vv) / (uu * uu + vv * vv);
    }
  }

  /* ---------------------------------------------------------------------------- */

  inline void getContributionFunction(int n, double *sf, double *tau, double *op, double *C)
  {
    
    /* --- Now compute the contribution function C = S * exp(-tau) * alpha --- */
    
    for(int kk = 0; kk < n; kk++) C[kk] = sf[kk] * exp(-tau[kk]) * op[kk];
    
    
  }
  
  
  /* ---------------------------------------------------------------------------- */

  void getTauNu(int n, double *opac, double *z, double *dtau, int &k0, int &k1, float &z_tau, double mu)
  {

    /* --- Cubic Bezier interpolation of the opacity --- */
    
    k0 = 0, k1 = n-1;

    dtau[0] = 0.;
    double odz = (z[0] - z[1]), oder = (opac[0] - opac[1]) / odz;
    double dz = 0.0, der = 0.0, dop=0.0, odop = oder, tau = 0.0, otau = 0.0;

    bool reached = false;
    
    for(int kk=1; kk<(n-1);kk++){
      int ku = kk-1, kd=kk+1;
      dz = (z[kk]-z[kd]), der = (opac[kk]-opac[kd])/dz;

      
      if(oder*der > 0.0){
	double lambda = (1.0 + dz / (dz + odz)) / 3.0;
	dop = (oder / (lambda * der + (1.0 - lambda) * oder)) * der;
      } else der = 0.0;
      
	/* --- 
	   integrate opacity using cubic bezier splines
	   The Bezier3 integral should be:
	   dz * (op_0 + op_u + cntrl1 + cntr2) / 4 
	   --- */
      
      dtau[kk] = odz * ((opac[kk] - dop/3.0 * odz) +
		       (opac[ku] + odop/3.0 * odz) +
		       opac[kk] + opac[ku]) * 0.25 / mu;
      
      
      tau = otau + dtau[kk];

      /* --- Get z at tau = 1 --- */
      
      if((otau <= 1.0) && (tau >=1.0)){
	double xu = log10(otau), x0 = log10(tau), u = (0.0-xu)/(x0-xu);
	z_tau = float((1.0-u) * z[ku] + u * z[kk]);  
      }

      oder = der, odz = dz, odop = dop, otau = tau;

      /* --- Stop integrating if tau is too large --- */
      
      if(tau <= 1.e-4) k0 = kk;
      if(tau <= 100.0) k1 = kk;
      else{
	reached = true;
	break;
      }
    }

    /* --- integrate last point? --- */

    if(!reached){
      int ku = n-2, kk = n-1;
      dtau[kk] = odz * ((opac[kk] - oder/3.0 * odz) +
			(opac[ku] + odop/3.0 * odz) +
			opac[kk] + opac[ku]) * 0.25 / mu;
      
      tau = otau+dtau[kk];
      
      k1 = kk;
      if((otau <= 1.0) && (tau >=1.0)){
	double xu = log10(otau), x0 = log10(tau), u = (0.0-xu)/(x0-xu);
	z_tau = float((1.0-u) * z[ku] + u * z[kk]);  
      }
    }
  }
  
/* ---------------------------------------------------------------------------- */

  void linear_int(int ndep, double *z, double *op, double *sf, double &syn,
		  double mu, float &tau_eq_1, double *C)
{

  double *dtau = new double [ndep](), dz = 0.0, eps, c0, cu, u0, u1;
  int kdep = ndep-1, kup = 0;

  
  /* --- Integrate tau_nu scale --- */
  
  getTauNu(ndep, op, z, dtau, kup, kdep, tau_eq_1, mu);

  
  /* --- Init emerging intensity --- */
  
  syn = sf[kdep] - (sf[kdep-1] - sf[kdep]) / dtau[kdep];


  /* --- Integrate linearly --- */
  
  for(int k=kdep-1; k>=kup; k--){
    int ku = k+1;
    double dt = dtau[ku];
    
    // Get coeffs

    if(dt >= 1.e-2){
      eps = exp(-dt);
      u0 = 1.0 - eps;
      u1 = dt - u0;
    } else{
      double dt2 = dt*dt, dt3 = dt2*dt;
      eps = 1.0 - dt + dt2 * 0.5 - dt3 / 6.0;
      u1 = (dt2 * 0.5) - dt3/6.0;
      u0 = dt-u1;
    }
    c0 = u1 / dt, cu = u0 - c0;

    syn = syn*eps + sf[k]*c0 + sf[ku]*cu;
  }


  /* --- Compute contribution function --- */
  
  if(C){
    
    double *tau = new double [ndep];
    tau[0] = 1.e-10;
    for(int kk=1; kk<ndep; kk++) tau[kk] = tau[kk-1] + dtau[kk];
    
    getContributionFunction(kdep-kup+1, &sf[kup], &dtau[kup], &op[kup], &C[kup]);
    delete [] tau;
    
  }
  
  /* --- Cleanup --- */
  
  delete [] dtau;
  
}

/* ---------------------------------------------------------------------------- */

void bezier3_int(int ndep, double *z, double *op, double *sf, double &syn,
			double mu, float &tau_eq_1, double *C)
{
  
  double *dtau = new double [ndep]();
  int kup = 0, kdep = ndep-1;
  
  
  /* --- Integrate tau_nu scale --- */
  
  getTauNu(ndep, op, z, dtau, kup, kdep, tau_eq_1, mu);

  
  /* --- Init the intensity with the value of the source function at the lowest point --- */
  
  syn = sf[kdep] - (sf[kdep-1] - sf[kdep]) / dtau[kdep];


  /* --- Get the derivatives of the source function with heigt --- */
  
  double *dsf = new double [ndep], *tau = new double [ndep];
  tau[0] = 1.e-10;
  for(int kk=1; kk<ndep; kk++) tau[kk] = tau[kk-1] + dtau[kk];
  
  cent_der<double>(ndep, tau, sf, dsf);

  
  /* --- Integrate ray --- */
 
  for(int k=kdep-1; k>=kup; k--){

    int ku = k + 1;

    /* --- Integration coeffs. and exponential --- */
    
    double dt = dtau[ku];
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt03 = dt / 3.0;
    double eps, alp, bet, gam, mu;
    //
    if(dt >= 1.e-3){
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

  /* --- Compute contribution function --- */
  fprintf(stderr,"%d %d %d %e\n",kdep-kup+1, kup, kdep, z[kdep]*1.e-5);
  if(C) getContributionFunction(kdep-kup+1, &sf[kup], &tau[kup], &op[kup], &C[kup]);
  

  /* --- Cleanup --- */

  delete [] dsf;
  delete [] dtau;
  delete [] tau;

 }
  
};

#endif
