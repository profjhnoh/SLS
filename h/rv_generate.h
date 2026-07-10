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

// ====================================================================
// Deterministic per-work-item RNG seeding (reproducibility fix)
// ====================================================================
// A parallel work item (e.g. one UE's HARQ decode in a given drop/slot) must
// draw the SAME random stream no matter which OpenMP thread runs it or in what
// order. Seeding by thread-id (get_thread_local_rng below) does NOT achieve this:
// OpenMP work distribution varies per run, so the same UE gets different coins
// and the PF feedback loop diverges (run-to-run SE noise). Seeding instead from a
// hash of the work item's identity coordinates makes the stream a pure function
// of (base_seed, drop, slot, ue, ...) → bit-reproducible across runs AND across
// thread counts (single-thread == multi-thread). Rand construction is cheap, so a
// fresh per-work-item generator per call is fine; draw sequentially within the call.
inline unsigned long long splitmix64_mix(unsigned long long x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
inline Rand make_workitem_rng(unsigned long long a, unsigned long long b,
                              unsigned long long c, unsigned long long d) {
  unsigned long long s = splitmix64_mix(a);
  s = splitmix64_mix(s ^ (b + 0x9E3779B97F4A7C15ULL));
  s = splitmix64_mix(s ^ (c + 0x9E3779B97F4A7C15ULL));
  s = splitmix64_mix(s ^ (d + 0x9E3779B97F4A7C15ULL));
  return Rand(s);
}

// ====================================================================
// Thread-local RNG support for OpenMP parallel regions
// ====================================================================
// Returns a thread-local RNG instance to avoid race conditions.
// Each thread gets its own RNG seeded with (base_seed + thread_id).
// Usage in OpenMP parallel regions:
//   Rand& local_rng = get_thread_local_rng();
//   double value = local_rng.u();
// ====================================================================
#ifdef ENABLE_MULTITHREADING
#include <omp.h>

inline Rand& get_thread_local_rng() {
    static thread_local Rand* local_rng = nullptr;
    static thread_local bool initialized = false;

    if (!initialized) {
        // Seed with base seed + thread ID for reproducibility
        extern unsigned long _seed;
        int thread_id = omp_get_thread_num();
        local_rng = new Rand((unsigned long long int)(_seed + thread_id * 1000000));
        initialized = true;
    }
    return *local_rng;
}
#endif

#endif


