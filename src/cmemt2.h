#ifndef __CMEMT2_H__
#define __CMEMT2_H__

/* ------------------------------------------------------------------------------------
   
  CMEMT2 template class
  
  PURPOSE:      Multi-dimensional dynamic arrays with C++ operators (up to 6 dimensions).
  AUTHOR:       Jaime de la Cruz Rodriguez (ISP-SU 2016).
  
  MODIFICATIONS:
           2016.12.23, JdlCR: Fixed indexing bug and made indexing methods const. 
	                      Does it help when dealing with multiple threads (?).
  
 ------------------------------------------------------------------------------------ */

#include <iostream>
#include <vector>
#include <cstring>
#include <algorithm>

template <class T> class mat {
 protected:
  int ndim;
  long unsigned nelem;
  unsigned n[6];
  size_t nc[6];
  bool alloc;
 public:
  T *d;  
  
  
  /* ---- Constructors & Destructors --- */
  
 mat():d(NULL),ndim(0),nelem(0),alloc(false){};
 mat(int n0, bool allocate = true):d(NULL),ndim(0),nelem(0),alloc(false)
    {init(1,0,0,0,0,0,n0, allocate);};
 mat(int n0, int n1, bool allocate = true):d(NULL),ndim(0),nelem(0),alloc(false)
    {init(2,0,0,0,0,n0,n1, allocate);};
 mat(int n0, int n1, int n2, bool allocate = true):d(NULL),ndim(0),nelem(0),alloc(false)
    {init(3,0,0,0,n0,n1,n2, allocate);};
 mat(int n0, int n1, int n2, int n3, bool allocate = true):d(NULL),ndim(0),nelem(0),alloc(false)
    {init(4,0,0,n0,n1,n2,n3, allocate);};
 mat(int n0, int n1, int n2, int n3, int n4, bool allocate = true):d(NULL),ndim(0),nelem(0),alloc(false)
    {init(5,0,n0,n1,n2,n3,n4, allocate);};
 mat(int n0, int n1, int n2, int n3, int n4, int n5, bool allocate = true):d(NULL),ndim(0),nelem(0),alloc(false)
    {init(6,n0,n1,n2,n3,n4,n5, allocate);};
 mat( const mat<T> &other):d(NULL),ndim(0),nelem(0),alloc(false) // = copy operator at initialization.
    {
      std::vector<unsigned> di = other.dims();
      size_t ndims = (size_t)di.size();
      const bool allocate = true;
      
      if(ndims == 1) init(1,0,0,0,0,0,di[0], allocate);
      else if(ndims == 2) init(2,0,0,0,0,di[0],di[1], allocate);
      else if(ndims == 3) init(3,0,0,0,di[0],di[1],di[2], allocate);
      else if(ndims == 4) init(4,0,0,di[0],di[1],di[2],di[3], allocate);
      else if(ndims == 5) init(5,0,di[0],di[1],di[2],di[3], di[4], allocate);
      else if(ndims == 6) init(6,di[0],di[1],di[2],di[3], di[4], di[5], allocate);
      else return;
      
      memcpy(this->d, other.d, this->nelem*sizeof(T));
      
    }
  
  //mat(mat<T> &m):d(NULL),ndim(0),nelem(0),alloc(false){
  //  unsigned dum[6];
  //}
 mat(std::vector<size_t> &di, bool allocate = true):d(NULL),ndim(0),nelem(0),alloc(false){
    size_t ndims = (size_t)di.size();
    if(ndims == 1) init(1,0,0,0,0,0,di[0], allocate);
    else if(ndims == 2) init(2,0,0,0,0,di[0],di[1], allocate);
    else if(ndims == 3) init(3,0,0,0,di[0],di[1],di[2], allocate);
    else if(ndims == 4) init(4,0,0,di[0],di[1],di[2],di[3], allocate);
    else if(ndims == 5) init(5,0,di[0],di[1],di[2],di[3], di[4], allocate);
    else if(ndims == 6) init(6,di[0],di[1],di[2],di[3], di[4], di[5], allocate);
  };

  /* ------------------------------------------------------------------------------------ */

  void cleanup()
  {
    if(alloc && (d != NULL)) delete [] d;
    d = NULL, ndim = 0, nelem = 0;
    memset(&n[0], 0, 6*sizeof(unsigned));
  }

  /* ------------------------------------------------------------------------------------ */

  ~mat(){
    cleanup();
  }
  
  /* ------------------------------------------------------------------------------------ */

  void rinit(int n0, int n1 = 0, int n2 = 0, int n3 = 0, int n4 = 0,
	     int n5 = 0, bool allocate = true)
  {

    int myndim = 0;
    if(n0 > 0) myndim++;
    if(n1 > 0) myndim++;
    if(n2 > 0) myndim++;
    if(n3 > 0) myndim++;
    if(n4 > 0) myndim++;
    if(n5 > 0) myndim++;

    init(myndim, n5, n4, n3, n2, n1, n0, allocate);
  }
  
  /* ------------------------------------------------------------------------------------ */

  void rinit(std::vector<int> dim, bool allocate = true)
  {

    int myndim = (int)dim.size();
    dim.resize(6, 0);
    init(myndim, dim[5], dim[4], dim[3], dim[2], dim[1], dim[0], allocate);
    
  }

  
  /* ------------------------------------------------------------------------------------ */

  void init(size_t ndimin, size_t n0, size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, bool allocate){
    
    memset(&n[0], 0, 6*sizeof(unsigned));
    memset(&nc[0], 0, 6*sizeof(size_t));

    if(ndimin == 0) return;
    
    n[0] = (unsigned)n0, n[1] = (unsigned)n1, n[2] = (unsigned)n2,
      n[3] = (unsigned)n3, n[4] = (unsigned)n4, n[5] = (unsigned)n5;

    nelem = 1;
    for(int ii=0; ii<6; ii++) if(n[ii] > 0) nelem *= n[ii];

    nc[5] = (size_t)n[5];
    for(int ii=4; ii>=0; ii--)nc[ii] = (size_t)n[ii]*nc[ii+1];
    
    
    ndim = ndimin;

    if(allocate){
      if(d != NULL) delete [] d;
      alloc = true;
      d = new T [nelem];
    } else{
      alloc = false;
      d = NULL;
    }
    
    
  }

  /* ------------------------------------------------------------------------------------ */

  inline void zero(){ memset(d, 0, nelem*sizeof(T));};
  
  /* ------------------------------------------------------------------------------------ */
  
  inline T &operator[](const size_t x) const { return d[x]; };
  
  /* ------------------------------------------------------------------------------------ */
  
  inline T &operator()(const size_t x) const { return d[x]; };
  
  /* ------------------------------------------------------------------------------------ */

  inline T &operator()(const size_t y, const size_t x) const {
    return d[ y * nc[5] + x];
  }

  /* ------------------------------------------------------------------------------------ */
  
  inline T &operator()(const size_t z, const size_t y, const size_t x) const{
    return d[z * nc[4] + y * nc[5] + x];
  }

  /* ------------------------------------------------------------------------------------ */

  inline T &operator()(const size_t t, const size_t z, const size_t y, const size_t x) const{
    return d[t * nc[3]  + z * nc[4] + y * nc[5] + x];
  }

  /* ------------------------------------------------------------------------------------ */

  inline T &operator()(const size_t r, const size_t t, const size_t z, const size_t y, const size_t x) const{
    return d[r * nc[2] + t * nc[3]  + z * nc[4] + y * nc[5] + x];
  }

  /* ------------------------------------------------------------------------------------ */

  inline T &operator()(const size_t w, const size_t r, const size_t t, const size_t z, const size_t y, const size_t x) const{
    return d[w * nc[1] + r * nc[2] + t * nc[3]  + z * nc[4] + y * nc[5] + x];
  }

  /* ------------------------------------------------------------------------------------ */

  inline mat<T> &operator= (const mat<T> &other) { // assignment operator
    
      std::vector<unsigned> di = other.dims();
      size_t ndims = (size_t)di.size();
      const bool allocate = true;
      
      if(ndims == 1) init(1,0,0,0,0,0,di[0], allocate);
      else if(ndims == 2) init(2,0,0,0,0,di[0],di[1], allocate);
      else if(ndims == 3) init(3,0,0,0,di[0],di[1],di[2], allocate);
      else if(ndims == 4) init(4,0,0,di[0],di[1],di[2],di[3], allocate);
      else if(ndims == 5) init(5,0,di[0],di[1],di[2],di[3], di[4], allocate);
      else if(ndims == 6) init(6,di[0],di[1],di[2],di[3], di[4], di[5], allocate);
      else return *this;
      
      memcpy(this->d, other.d, this->nelem*sizeof(T));
    
    return *this;
  }

  /* ------------------------------------------------------------------------------------ */

  inline void operator+= (const mat<T> &m){
    for(size_t ii = 0; ii < this->nelem; ii++) this->d[ii] += m.d[ii];
  }

  /* ------------------------------------------------------------------------------------ */

  inline mat<T> operator+ (const mat<T> &m){
    for(size_t ii=0;ii<this->nelem;ii++) this->d[ii] += m.d[ii];
    return *this;
  }

 /* ------------------------------------------------------------------------------------ */

  inline void operator-= (const mat<T> &m){
    for(size_t ii = 0; ii < this->nelem; ii++) this->d[ii] -= m.d[ii];
  }

  /* ------------------------------------------------------------------------------------ */

  inline mat<T> operator- (const mat<T> &m){
    for(size_t ii=0;ii<this->nelem;ii++) this->d[ii] -= m.d[ii];
    return *this;
  }

 /* ------------------------------------------------------------------------------------ */  

  inline void operator*= (const mat<T> &m){
    for(size_t ii = 0; ii < this->nelem; ii++) this->d[ii] *= m.d[ii];
  }

  /* ------------------------------------------------------------------------------------ */

  inline mat<T> operator* (const mat<T> &m){
    for(size_t ii=0;ii<this->nelem;ii++) this->d[ii] *= m.d[ii];
    return *this;
  }

 /* ------------------------------------------------------------------------------------ */
  
  inline void operator/= (const mat<T> &m){
    for(size_t ii = 0; ii < this->nelem; ii++) this->d[ii] /= m.d[ii];
  }

  /* ------------------------------------------------------------------------------------ */

  inline mat<T> operator/ (const mat<T> &m){
    for(size_t ii=0;ii<this->nelem;ii++) this->d[ii] /= m.d[ii];
    return *this;
  }

 /* ------------------------------------------------------------------------------------ */
  
  unsigned long size() const{
    return nelem;
  }

  /* ------------------------------------------------------------------------------------ */

  T sum() const{
    T suma = 0;
    for(size_t ii = 0; ii<nelem; ii++) suma += d[ii];
    return suma;
  }

  /* ------------------------------------------------------------------------------------ */

  double dsum() const{
    long double suma = 0.0, c = 0.0;
    for(size_t kk = 0; kk<nelem; kk++){
      long double y = d[kk] - c;
      long double t = suma + y;
      c = (t - suma) - y;
      suma = t;
    }
    return (double)suma;
  }
  
  /* ------------------------------------------------------------------------------------ */

  double mean() const{
    return dsum() / (double)nelem;
  }

  /* ------------------------------------------------------------------------------------ */

  std::vector<unsigned> dims() const{

    std::vector<unsigned> idims;
    idims.resize(ndim);
    for(int ii=0; ii<ndim; ii++)
      idims[ii] = (unsigned)n[6-ndim+ii];
    
    return idims;
  }

  /* ------------------------------------------------------------------------------------ */

  inline unsigned shape(size_t idim) const{
    int idim1 = (int)idim;
    idim1 += 6 - ndim;
    if(idim1 > 5) idim1 = 5;
    if(idim1 < 0) idim1 = 0;
    return n[idim1];
  }
  
  /* ------------------------------------------------------------------------------------ */

  inline int ndims() const{
    return ndim;
  }
  
  /* ------------------------------------------------------------------------------------ */

  inline void operator+= (T c){
    for(size_t ii = 0; ii < this->nelem; ii++) this->d[ii] += c;
  }

  /* ------------------------------------------------------------------------------------ */

  inline void operator-= (T c){
    for(size_t ii = 0; ii < this->nelem; ii++) this->d[ii] -= c;
  }

  /* ------------------------------------------------------------------------------------ */

  inline void operator*= (T c){
    for(size_t ii = 0; ii < this->nelem; ii++) this->d[ii] *= c;
  }

  /* ------------------------------------------------------------------------------------ */
  
  inline void operator/= (T c){
    for(size_t ii = 0; ii < this->nelem; ii++) this->d[ii] /= c;
  }

  /* ------------------------------------------------------------------------------------ */
  
  inline void operator= (T c){
    for(size_t ii = 0; ii < this->nelem; ii++) this->d[ii] = c;
  }
  
};

#endif
