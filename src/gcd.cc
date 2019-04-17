#include <fflow/gcd.hh>
#include <fflow/mp_gcd.hh>

namespace fflow {

  void rat_rec(UInt a, UInt b, Int & num, Int & den)
  {
    MPInt mpa(a), mpb(b);
    MPRational res;
    if (rat_rec(mpa, mpb, res))
      res.to_int(num, den);
    else
      num = den = 0;
  }

} // namespace fflow
