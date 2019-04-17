#include <fflow/spinor.hh>

namespace fflow {

  class Basis6D {
  public:
    Basis6D ()
      : e1(), e2(), e3(), e4(), e5(), e6(),
        s12(0), s34(0), s56(0) {}

    Basis6D (const Basis6D & oth)
      : e1(oth.e1), e2(oth.e2), e3(oth.e3), e4(oth.e4),
        e5(oth.e5), e6(oth.e6),
        s12(oth.s12), s34(oth.s34), s56(oth.s56) {}

    Basis6D (const Spinor4D & sp1, const Spinor4D & sp2, Mod mod)
      : e1(Momentum4D(sp1,mod)), e2(Momentum4D(sp2,mod)),
        e3(Momentum4D(sp1,sp2,mod)), e4(Momentum4D(sp2,sp1,mod)),
        e5(0,0,0,0,1,0), e6(0,0,0,0,0,1),
        s12(sij(sp1,sp2,mod)), s34(), s56(mod-1)
    {
      s34 = neg_mod(s12,mod);
    }

  public:
    Momentum6D e1, e2, e3, e4, e5, e6;
    UInt s12, s34, s56;
  };

} // namespace fflow
