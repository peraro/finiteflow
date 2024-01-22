#ifndef FFLOW_GCD_HH
#define FFLOW_GCD_HH

#include <stdexcept>
#include <fflow/common.hh>

namespace fflow {

#ifdef FFLOW_USE_EXCEPTIONS
  struct DivisionByZero : public std::invalid_argument {
    DivisionByZero() : std::invalid_argument("Division by zero") {}
  };
#endif

  // gcd
  inline UInt gcd(UInt a, UInt b)
  {
    return n_gcd(a,b);
  }

  // multiplicative inverse a^(-1) mod n
  inline UInt mul_inv(UInt a, Mod n)
  {
    return n_invmod(a, n.n());
  }

  // a/b mod n
  inline UInt div_mod(UInt a, UInt b, Mod n)
  {
    b = mul_inv(b, n);
    return mul_mod(a, b, n);
  }

  // rational reconstruction: num,den such that num/den = a mod b
  void rat_rec(UInt a, UInt b, Int & num, Int & den);


  // Rational

  inline Rational rat_normal(Rational z)
  {
    UInt n = iabs(z.num), d = iabs(z.den);
    if (!n)
      return Rational{0,1};
    UInt g = gcd(n,d);
    return Rational{sign(z.num)*sign(z.den)*Int(n/g), Int(d/g)};
  }

  inline Rational rat_rec(UInt a, UInt b)
  {
    Rational res;
    rat_rec(a,b,res.num,res.den);
    return res;
  }

} // namespace fflow

#endif // FFLOW_GCD_HH
