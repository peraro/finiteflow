/*
 * We use 63-bit primes, satisfying 2p < 2^64 < 3p.  This assumption
 * is made everywhere.
 */

#ifndef FF_MOD_H
#define FF_MOD_H

#include <stdint.h>
#include <fflow/config.hh>

#ifdef __cplusplus
extern "C" {
#endif

  typedef uint64_t FFU64;
  typedef int64_t FFI64;

  typedef unsigned __int128 FFU128;
  typedef union {
    struct {
#if FFLOW_BIG_ENDIAN
      FFU64 high, low;
#else
      FFU64 low, high;
#endif
    };
    FFU128 z;
  } FFU128HL;


  typedef struct {
    // value
    FFU64 n;
    // precomputed reciprocal of n<<1
    FFU64 v;
  } FFMod;

  // n must satisfy 2n < 2^64 < 3n.
  inline FFMod ffPrecomputedReciprocalMod(FFU64 n);

  // z mod p
  inline FFU64 ffDiv2Mod1(FFU128 z, FFMod p);

  // z = quo * p + rem
  inline void ffDiv2By1(FFU128 z, FFMod p, FFU64 * quo, FFU64 * rem);

  inline FFU64 ffAddMod(FFU64 a, FFU64 b, FFMod p);
  inline FFU64 ffSubMod(FFU64 a, FFU64 b, FFMod p);
  inline FFU64 ffMulMod(FFU64 a, FFU64 b, FFMod p);

  // a + b*c mod p
  inline FFU64 ffAPBCMod(FFU64 a, FFU64 b, FFU64 c, FFMod p);

  inline FFU64 ffPrecomputedMulShoup(FFU64 w, FFMod p);
  inline FFU64 ffMulShoup(FFU64 t, FFU64 w, FFU64 wp, FFMod p);

  FFU64 ffMulInverse(FFU64 z, FFMod p);


  // Implementation of inline functions


  /*
   * The next three functions are based on
   *
   *   Improved Division by Invariant Integers,
   *    by Niels Moller and Torbjorn Granlund
   *
   * The main difference is that, in that paper, they divide by a
   * number `n` satisfying n >= 2^63.  Hence, in our division
   * algorithm, we first multiply `n` and the `dividend` by 2 (with a
   * shift) to satisfy this requirement, and then we shift back the
   * reminder.
   */

  inline FFMod ffPrecomputedReciprocalMod(FFU64 n)
  {
    FFU64 p = n << 1;
    FFU128HL hl;
    hl.z = ~(FFU128)(0);
    hl.high -= p;
    hl.z /= p;
    FFMod res = {n, hl.low};
    return res;
  }

  inline FFU64 ffDiv2Mod1(FFU128 z, FFMod p)
  {
    FFU128HL u, q;
    FFU64 n = p.n << 1;
    u.z = z << 1;
    q.z = (FFU128)(p.v) * u.high + u.z;
    ++q.high;
    FFU64 r = u.low - q.high * n;
    if (r > q.low) {
      r += n;
    }
    if (r >= n) {
      r -= n;
    }
    return r >> 1;
  }

  inline void ffDiv2By1(FFU128 z, FFMod p, FFU64 * quo, FFU64 * rem)
  {
    FFU128HL u, q;
    FFU64 n = p.n << 1;
    u.z = z << 1;
    q.z = (FFU128)(p.v) * u.high + u.z;
    ++q.high;
    FFU64 r = u.low - q.high * n;
    if (r > q.low) {
      --q.high;
      r += n;
    }
    if (r >= n) {
      ++q.high;
      r -= n;
    }
    *quo = q.high;
    if (rem)
      *rem = r >> 1;
  }

  inline FFU64 ffMulMod(FFU64 a, FFU64 b, FFMod p)
  {
    FFU128 z = (FFU128)(a) * b;
    return ffDiv2Mod1(z, p);
  }

  // TODO: test add/sub with full flint

  inline FFU64 ffAddMod(FFU64 a, FFU64 b, FFMod p)
  {
#if !defined(__arm64__) && !defined(__arm64)
    FFU64 ret = a + b - p.n;
    return (FFI64)(ret) < 0 ? ret + p.n : ret;
#else
    FFU64 ret = a + b;
    return ret < p.n ? ret : ret - p.n;
#endif
  }

  inline FFU64 ffSubMod(FFU64 a, FFU64 b, FFMod p)
  {
#if !defined(__arm64__) && !defined(__arm64)
    FFI64 ret = a - b;
    return (FFI64)(ret) < 0 ? ret + p.n : ret;
#else
    return a >= b ? a - b : a + p.n - b;
#endif
  }

  // a + b*c mod p
  inline FFU64 ffAPBCMod(FFU64 a, FFU64 b, FFU64 c, FFMod p)
  {
    FFU128 z = (FFU128)(b) * c + a;
    return ffDiv2Mod1(z, p);
  }


  // Shoup’s modular multiplication algorithm.
  //
  // The following routines are taken directly from
  // http://web.maths.unsw.edu.au/~davidharvey/talks/fastntt-2-talk.pdf
  //
  // Given w in Z_p, with p<2^63, if we need to compute
  //
  //   t * w mod p
  //
  // for several values of t in Z_p, we can precompute
  //
  //   wp = quotient(2^64 * w, p)
  //
  // with ffPrecomputedMulShoup, and then use the fast multiplication
  // algoritm ffMulShoup.

  inline FFU64 ffPrecomputedMulShoup(FFU64 w, FFMod p)
  {
    FFU128HL hl;
    hl.high = w;
    hl.low = 0;
    FFU64 wp;
    ffDiv2By1(hl.z, p, &wp, 0);
    //wp = hl.z/p.n;
    return wp;
  }

  inline FFU64 ffMulShoup(FFU64 t, FFU64 w, FFU64 wp, FFMod p)
  {
    FFU128HL q;
    q.z = (FFU128)(wp) * t;

    FFU64 r = w*t - q.high*p.n;

    if (r >= p.n)
      r -= p.n;

    return r;
  }


#ifdef __cplusplus
}
#endif

#endif // FF_MOD_H
