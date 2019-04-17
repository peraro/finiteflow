#ifndef FFLOW_MP_GCD_HH
#define FFLOW_MP_GCD_HH

#include <fflow/mp_common.hh>

namespace fflow {

  inline bool mul_inv(const MPInt & z, const MPInt & mod,
                      MPInt & res)
  {
    return mpz_invert(res.get(), z.get(), mod.get());
  }

  // get an integer such that res = q mod p
  bool rat_mod(const MPRational & q, const MPInt & p,
               MPInt & res);

  // rational reconstruction in multiple precision
  bool rat_rec(const MPInt & z, const MPInt & mod,
               MPRational & res);

  // return c1, c2, mod12 such that the chinese remainder of (a%mod1)
  // and (b%mod2) is given by
  //
  //   (a*c1 + b*c2) % mod12
  //
  // for generic a,b
  void chinese_remainder_coeffs(const MPInt & mod1, UInt mod2,
                                MPInt & c1, MPInt & c2, MPInt & mod12);

  void chinese_remainder_from_coeffs(const MPInt & a, UInt b,
                                     const MPInt & c1, const MPInt & c2,
                                     const MPInt & mod12,
                                     MPInt & res);

} // namespace fflow


#endif // FFLOW_MP_GCD_HH
