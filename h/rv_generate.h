/*===================================================================
FUNCTION: rv_generate.h

DESCRIPTION: random number generator

NOTES: refered from 
       Numerical Recipes 3rd Edition: The Art of Scientific Computing (3/E) / William H. Press
       / Ch.7 Random numbers / p340~410
   
===================================================================*/
#include <complex>
#include <cmath>

#ifndef _Rand_
#define _Rand_

struct Rand {
  unsigned long long int p,q,r;
  Rand(unsigned long long int _seed) : q(4101842887655102017LL), r(1) {
    p = _seed ^ q; int64();
    q = p; int64();
    r = q; int64();
  }
  inline unsigned long long int int64() {
    p = p * 2862933555777941757LL + 7046029254386353087LL;
    q ^= q >> 17; q ^= q << 31; q ^= q >> 8;
    r = 4294957665U*(r & 0xffffffff) + (r >> 32);
    unsigned long long int x = p ^ (p << 21); x ^= x >> 35; x ^= x << 4;
    return (x + q) ^ r;
  }
  inline double u() { return 5.42101086242752217E-20 * int64(); }
  inline unsigned int i() { return (unsigned int) int64(); }
  inline double n() 
  {
    double nu,nv,nx,ny,nz;
    do
    {
      nu = u();
      nv = 1.7156*(u()-0.5);
      nx = nu - 0.449871;
      ny = fabs(nv) + 0.386595;
      nz = nx*nx + ny*(0.19600*ny-0.25472*nx);
    }
    while (nz > 0.27597 && (nz > 0.27846 || nv*nv > -4.*log(nu)*nu*nu));
    return nv/nu;
  }
} ;
 

#endif


