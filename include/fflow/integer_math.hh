#ifndef FFLOW_INTEGER_MATH_HH
#define FFLOW_INTEGER_MATH_HH

#include <fflow/common.hh>

namespace fflow {

  const UInt FACTORIAL[21] = {1ull,
                              1ull,
                              2ull,
                              6ull,
                              24ull,
                              120ull,
                              720ull,
                              5040ull,
                              40320ull,
                              362880ull,
                              3628800ull,
                              39916800ull,
                              479001600ull,
                              6227020800ull,
                              87178291200ull,
                              1307674368000ull,
                              20922789888000ull,
                              355687428096000ull,
                              6402373705728000ull,
                              121645100408832000ull,
                              2432902008176640000ull};

  UInt binomial(UInt a, UInt b);
  UInt binomial(UInt a, UInt b, Mod mod);
  UInt power(UInt a, unsigned b, Mod mod);


  namespace detail {
    UInt binomial_impl(UInt a, UInt b);
    UInt binomial_impl(UInt a, UInt b, Mod mod);
  }

  inline UInt binomial(UInt a, UInt b)
  {
    return a <= 20 ? FACTORIAL[a]/(FACTORIAL[a-b]*FACTORIAL[b])
      : detail::binomial_impl(a,b);
  }

  inline UInt binomial(UInt a, UInt b, Mod mod)
  {
    // note: if a<=20 FACTORIAL[a] is smaller than our smallest prime,
    // hence mod p is not needed
    return a <= 20 ? (FACTORIAL[a]/(FACTORIAL[a-b]*FACTORIAL[b]))
      : detail::binomial_impl(a,b,mod);
  }

} // namespace fflow

#endif // FFLOW_INTEGER_MATH_HH
