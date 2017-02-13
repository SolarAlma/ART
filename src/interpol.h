/*
  Various interpolation routines

  Coded by J. de la Cruz Rodriguez (ISP-SU 2015)
 */
#ifndef INTERPOL_H
#define INTERPOL_H
//
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <vector>

/* --------------------------------------------------------------------- */

/*
  Piece-wise linear interpolation, with optinal 
  linear extrapolation outside the input range.
 */ 
template <class T, class U> void linpol(size_t ni, T *x, T *y, size_t nni, U *xx, U *yy, bool extrapolate = false){

  unsigned n = (unsigned)ni, nn = (unsigned)nni;
  
  //
  // increasing or decreasing x?
  //
  bool dir = true;
  if( (x[1]-x[0]) < 0 ) dir = false;

  
  if(dir){

    unsigned off = 0;
    for(unsigned k = 1; k<n; k++){

      // Coeffs.
      double a = (y[k] - y[k-1]) / (x[k] - x[k-1]);
      double b = y[k-1] - a * x[k-1]; 

      // Check if there are points to compute
      for(unsigned j = off; j<nn; j++){

	if((xx[j] >= x[k-1]) && (xx[j] < x[k])){
	  yy[j] = a * xx[j] + b;
	  off++;
	}
	
      } // j
    } // k
    
  }else{
    
    unsigned off = 0;
    for(unsigned k = 1; k<n; k++){
      
      // Coeffs.
      double a = (y[k] - y[k-1]) / (x[k] - x[k-1]);
      double b = y[k-1] - a * x[k-1]; 
      
      // Check if there are points to compute
      for(unsigned j = off; j<nn; j++){
	
	if((xx[j] <= x[k-1]) && (xx[j] > x[k])){
	  yy[j] = a * xx[j] + b;
	  off++;
	}
	
      } // j 
    } // k
  } // Else

  
  double a0, a1, b0, b1;
  if(extrapolate){
    a0 = (y[1] - y[0]) / (x[1] - x[0]);
    b0 = y[0] - a0 * x[0];
    a1 = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
    b1 = y[n-1] - a1 * x[n-1];
  }else{
    a0 = 0;
    a1 = 0;
    b0 = y[0];
    b1 = y[n-1];
  }
  
  if(dir) for(unsigned k = 0; k<nn; k++){
      if(xx[k] <= x[0])   yy[k] = a0 * xx[k] + b0;
      if(xx[k] >= x[n-1]) yy[k] = a1 * xx[k] + b1;
    }
  else for(unsigned k = 0; k<nn; k++){
      if(xx[k] >= x[0])   yy[k] = a0 * xx[k] + b0;
      if(xx[k] <= x[n-1]) yy[k] = a1 * xx[k] + b1;
    }
} // linpol


/* --------------------------------------------------------------------- */


/*
  Piece-wise Hermitian interpolation, with optional 
  linear extrapolation outside the input range.
  
  Hermitian interpolant from Auer (2003), Formal Solution: Explicit answers.

 */ 
template <class T1, class T2> void hermpol(size_t ni, T1 *x, T1 *y, size_t nni, T2 *xx, T2 *yy, bool extrapolate = false){

  unsigned n = (unsigned)ni, nn = (unsigned)nni;
  
  // Increasing x?
  bool sign = true;
  if( (x[1] - x[0]) < 0 ) sign = false;

  
  // Init derivatives
  double odx = x[1] - x[0];
  double oder = (y[1] - y[0]) / odx;
  double ody = oder;
  double dy = 0;
  double der = 0;
  double dx = 0;
  unsigned off = 0;
  
  if(sign){

    for(unsigned k = 1; k<n; k++){
      
      // Derivatives
      if(k<(n-1)){
	dx =  x[k+1] - x[k];
	dy = (y[k+1] - y[k]) / dx;

	//if(dy*ody > 0) der = (dx * ody + odx * dy) / (dx + odx);
	double lambda = (1.0 + dx / (dx + odx)) / 3.0;
	if(dy*ody > 0) der = (dy*ody) / ((1.0-lambda)*ody + lambda * dy);
	else der = 0;
      } else der = ody;
      
      
      // Check if there are points to compute
      for(unsigned j = off; j<nn; j++){      
	
	if( (xx[j] >= x[k-1]) && (xx[j] < x[k]) ){
	  // Normalize interval units
	  double u  = (xx[j] - x[k-1]) / odx;
	  double uu = u*u;

	  // Hermitian interpolant
	  yy[j] = y[k-1] * (1.0 - 3*uu + 2*uu*u) + (3*uu - 2*uu*u) * y[k] +
	    (uu*u - 2*uu + u) * odx * oder + (uu*u - uu) * odx * der;
	  
	  off++;
	}
	
      } // j

      // Store values
      odx = dx;
      ody = dy;
      oder = der;
      
    } // k
  }else{
    for(unsigned k = 1; k<n; k++){
      
      // Derivatives
      if(k<(n-1)){
	dx =  x[k+1] - x[k];
	dy = (y[k+1] - y[k]) / dx;

	//if(dy*ody > 0) der = (dx * ody + odx * dy) / (dx + odx);
	double lambda = (1.0 + dx / (dx + odx)) / 3.0;
	if(dy*ody > 0) der = (dy*ody) / ((1.0-lambda)*ody + lambda * dy);
	else der = 0;
      } else der = ody;
      
      
      // Check if there are points to compute
      for(unsigned j = off; j<nn; j++){      
	
	if( (xx[j] <= x[k-1]) && (xx[j] > x[k]) ){

	  double u  = (xx[j] - x[k-1]) / odx;
	  double uu = u*u;
	  double uuu = uu*u;
	  
	  yy[j] = y[k-1] * (1.0 - 3*uu + 2*uuu) + (3*uu - 2*uuu) * y[k] +
	    (uuu - 2*uu + u) * odx * oder + (uuu - uu) * odx * der;
	  
	  off++;
	}
	
      } // j
      
      odx = dx;
      ody = dy;
      oder = der;
      
  } // k
  }


  //
  // Points outside the x[0], x[n-1]?
  //
  double a0, a1, b0, b1;
  if(extrapolate){
    a0 = (y[1] - y[0]) / (x[1] - x[0]);
    b0 = y[0] - a0 * x[0];
    a1 = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
    b1 = y[n-1] - a1 * x[n-1];
  }else{
    a0 = 0;
    a1 = 0;
    b0 = y[0];
    b1 = y[n-1];
  }
  
  if(sign){
    for(unsigned k = 0; k<nn; k++){
      if(xx[k] <= (x[0]))   yy[k] = a0 * xx[k] + b0;
      if(xx[k] >= (x[n-1])) yy[k] = a1 * xx[k] + b1;
    }
  } else {
    for(unsigned k = 0; k<nn; k++){
      if(xx[k] >= (x[0]))   yy[k] = a0 * xx[k] + b0;
      if(xx[k] <= (x[n-1])) yy[k] = a1 * xx[k] + b1;
    }
  }
  
} // hermpol

/* --------------------------------------------------------------------- */

template <class T1, class T2> void cpol( T1 y, int nn, T2 *yy){

  T2 dum = (T2)y;
  for(int kk = 0; kk<nn; kk++) yy[kk] = dum;
  
  
}

/* --------------------------------------------------------------------- */

template <class T> T sqr(T var){
  return var*var;
}

/* --------------------------------------------------------------------- */

template <class T> std::vector<double> parab_fit(T d, T e, T f, T yd, T ye, T yf){


  std::vector<double> cf;
  cf.resize(3,0.0);

  cf[1] = ((yf - yd) - (f*f - d*d) * ((ye - yd) / (e*e - d*d)))/ 
    ((f - d) - (f*f - d*d) * ((e - d) / (e*e - d*d)));
  cf[2] = ((ye - yd) - cf[1] * (e - d)) / (e*e - d*d);
  cf[0] = yd - cf[1] * d - cf[2] * d*d;

  return cf;
}

/* --------------------------------------------------------------------- */

template <class T1, class T2> void bezpol2(size_t ni, T1 *x, T1 *y, size_t nni, T2 *xx, T2 *yy, bool extrapolate = false){

  unsigned n = (unsigned)ni, nn = (unsigned)nni;

  // Increasing x?
  bool sign = true;
  if( (x[1] - x[0]) < 0 ) sign = false;

  
  // Init derivatives
  double odx = x[1] - x[0];
  double oder = (y[1] - y[0]) / odx;
  double ody = oder;
  double dy = 0;
  double der = 0;
  double dx = 0;
  unsigned off = 0;
  
  if(sign){

    for(unsigned k = 1; k<n; k++){
      
      /* --- Compute centered derivatives --- */
      
      if(k<(n-1)){
	dx =  x[k+1] - x[k];
	dy = (y[k+1] - y[k]) / dx;

	double lambda = (1.0 + dx / (dx + odx)) / 3.0;
	if(dy*ody > 0) der = (dy*ody) / ((1.0-lambda)*ody + lambda * dy);
	else der = 0;
	//der = (dy*ody) / ((1.0-lambda)*ody + lambda * dy);
      } else der = ody;
      


      
      /* --- Make a combined control point from the upwind and central points, 
	 makes a very tight spline compared to only using one of the control points
	 --- */
      double cntrl = 0.5 * (  (y[k-1] + 0.5*odx*oder)  +  (y[k] - 0.5*odx*der)  );

      
      /* ---  Check if there are points to compute --- */
      for(unsigned j = off; j<nn; j++){      
	
	if( (xx[j] >= x[k-1]) && (xx[j] < x[k]) ){
	  // Normalize interval units
	  double u  = (xx[j] - x[k-1]) / odx;
	  double u1 = 1.0 - u;
	  
	  // Quadratic Bezier interpolant
	  yy[j] = y[k-1]*u1*u1 + y[k]*u*u + 2.0*cntrl*u*u1; 
	  
	  off++;
	}
	
      } // j

      /* --- Store values for next interval --- */
      
      odx = dx;
      ody = dy;
      oder = der;
      
    } // k
  }else{
    for(unsigned k = 1; k<n; k++){
      
      // Derivatives
      if(k<(n-1)){
	dx =  x[k+1] - x[k];
	dy = (y[k+1] - y[k]) / dx;

	double lambda = (1.0 + dx / (dx + odx)) / 3.0;
	if(dy*ody > 0) der = (dy*ody) / ((1.0-lambda)*ody + lambda * dy);
	else der = 0;
	//der = (dy*ody) / ((1.0-lambda)*ody + lambda * dy);
      } else der = ody;



      
      /* --- Make a combined control point from the upwind and central points, 
	 makes a very tight spline compared to only using one of the control points
	 --- */
      double cntrl = 0.5 * (  (y[k-1] + 0.5*odx*oder)  +  (y[k] - 0.5*odx*der)  );

      
      
      // Check if there are points to compute
      for(unsigned j = off; j<nn; j++){      
	
	if( (xx[j] <= x[k-1]) && (xx[j] > x[k]) ){

	  double u  = (xx[j] - x[k-1]) / odx;
	  double u1 = 1.0 - u;
	  
	  // Quadratic Bezier interpolant
	  yy[j] = y[k-1]*u1*u1 + y[k]*u*u + 2.0*cntrl*u*u1; 
	  
	  off++;
	}
	
      } // j
      
      odx = dx;
      ody = dy;
      oder = der;
      
  } // k
  }


  //
  // Points outside the x[0], x[n-1]?
  //
  double a0, a1, b0, b1;
  if(extrapolate){
    a0 = (y[1] - y[0]) / (x[1] - x[0]);
    b0 = y[0] - a0 * x[0];
    a1 = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
    b1 = y[n-1] - a1 * x[n-1];
  }else{
    a0 = 0;
    a1 = 0;
    b0 = y[0];
    b1 = y[n-1];
  }
  
  if(sign){
    for(unsigned k = 0; k<nn; k++){
      if(xx[k] <= (x[0]))   yy[k] = a0 * xx[k] + b0;
      if(xx[k] >= (x[n-1])) yy[k] = a1 * xx[k] + b1;
    }
  } else {
    for(unsigned k = 0; k<nn; k++){
      if(xx[k] >= (x[0]))   yy[k] = a0 * xx[k] + b0;
      if(xx[k] <= (x[n-1])) yy[k] = a1 * xx[k] + b1;
    }
  }
  
} // bezpol2


/* --------------------------------------------------------------------- */

template <class T1, class T2> void vlint(size_t n, T1 *x, T1 *y, size_t n1, T2 *xx, T2 *yy)
{

  double x0=0.0, y0=0.0, x1=0.0, y1=0.0, dxu=0.0, dxd=0.0,
    r=0.0, phi=0.0, der=0.0, u = 0.0;
  size_t off = 0, off0 = 0, off1 = 0;
  
  for(size_t ii=1;ii<(n-1); ii++){

    /* --- Compute the derivative --- */
    
    dxu = (x[ii]-x[ii-1])*0.5;
    dxd = (x[ii+1]-x[ii])*0.5;
    //
    r   = (y[ii]-y[ii-1])/(y[ii+1]-y[ii]);
    if(!std::isfinite(r)) r = 0.0;
    //
    phi = (r + fabs(r)) / (1.0 + fabs(r));
    der = phi * (y[ii+1]-y[ii]) / (dxd+dxu);

    /* --- Now check if there are points in this interval --- */

    bool first = true;
    for(size_t jj=off; jj<n1; jj++){
      if( (xx[jj] <= (x[ii]+dxd)) && (xx[jj] > (x[ii]-dxu)) ){
	yy[jj] = der * (xx[jj]-x[ii]) + y[ii];

	if((ii==1) && first) x0 = xx[jj], y0 = yy[jj];
	x1 = xx[jj], y1 = yy[jj];
	
	first = false;
	off++;
      } // if
    }// jj
    
  }// ii


    /* --- Second interval ---*/
  
  dxu = (x[n-1] - x[n-2]) * 0.5;
  for(size_t jj=off;jj<n1; jj++){
    if(xx[jj] >= (x[n-1]-dxu)) yy[jj] = (y[n-1]-y1)/(x[n-1]-x1) * (xx[jj]-x[n-1]) + y[n-1];
  }
  
  
  /* --- first interval --- */

  dxd = (x[1] - x[0]) * 0.5;
  der = (y[0]-y0)/(x[0]-x0);
  for(size_t jj=0;jj<n1; jj++){
    if(xx[jj] <= (x[0]+dxd)) yy[jj] =  der * (xx[jj]-x[0]) + y[0];
    else break;
  }
  
}

/* --------------------------------------------------------------------- */

template <class T> void cent_der(size_t n, T *x, T *y, T *yp)
{

  double odx = x[1]-x[0], dx = 0;
  double oder = (y[1]-y[0])/odx, der = 0;
  yp[0] = oder;

  for(size_t ii=1;ii<(n-1); ii++){
    dx = x[ii+1] - x[ii];
    der = (y[ii+1] - y[ii]) / dx;

    if(der*oder > 0.0){
      double lambda = (1.0 + dx / (dx + odx)) / 3.0;
      yp[ii] = (oder / (lambda * der + (1.0 - lambda) * oder)) * der;
    } else yp[ii] = 0.0; // Set der to zero at extrema;

    /* --- Save values --- */
    
    oder = der, odx = dx;
  }

  yp[n-1] = oder;
}



#endif
