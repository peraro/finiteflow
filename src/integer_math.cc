#include <fflow/integer_math.hh>
#include <fflow/gcd.hh>
#include <algorithm>

namespace fflow {

  namespace detail {

    UInt binomial_impl(UInt a, UInt b)
    {
      UInt res = 1;
      for (UInt i=1; i<=std::min(b,a-b); ++i) {
        res *= (a+1-i);
        res /= i;
      }
      return res;
    }

    UInt binomial_impl(UInt a, UInt b, Mod mod)
    {
      UInt num = 1, den = 1;
      for (UInt i=1; i<=std::min(b,a-b); ++i) {
        num = mul_mod(num, (a+1-i), mod);
        den = mul_mod(den, i, mod);
      }
      return div_mod(num, den, mod);
    }

  } // namespace detail


  UInt power(UInt x, unsigned y)
  {
    UInt res = 1;
    while (y) {
      if (y & 1)
        res *= x;
      x *= x;
      y >>= 1;
    }
    return res;
  }

  UInt power(UInt x, unsigned y, Mod mod)
  {
    UInt res = 1;
    while (y) {
      if (y & 1)
        res = mul_mod(res, x, mod);
      x = mul_mod(x, x, mod);
      y >>= 1;
    }
    return res;
  }


} // namespace fflow
