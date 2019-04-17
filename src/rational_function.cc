#include <fflow/rational_function.hh>
#include <algorithm>

namespace fflow {

  void URationalFunction::normalize(Mod mod)
  {
    std::size_t i = 0, inum=0;
    UInt norm;

    // look for non-zero entry in denominator
    while (i<=getDenDegree() && den_[i] == 0)
      ++i;

    // if den_=0, return 0/0
    if (i == getDenDegree()+1) {
      clear();
      return;
    }

    // look for nonzero entry in numerator
    while (inum<=getNumDegree() && num_[inum] == 0)
      ++inum;

    // if num_=0, return 0/1
    if (inum == getNumDegree()+1) {
      clear();
      den_[0] = 1;
      return;
    }

    // if both num and den have a factor x^k, simplify it out
    if (inum && i) {
      unsigned r = std::min(inum,i);
      num_.shift_left(r);
      den_.shift_left(r);
      i -= r;
      //inum -= r;
    }

    norm = mul_inv(den_[i], mod);
    den_[i] = 1;

    for (unsigned j=i+1; j<=getDenDegree(); ++j)
      den_[j] = mul_mod(den_[j], norm, mod);
    for (unsigned j=0; j<=getNumDegree(); ++j)
      num_[j] = mul_mod(num_[j], norm, mod);

    num_.remove_high_zeroes();
    den_.remove_high_zeroes();
  }


  std::string URationalFunction::to_str(const std::string & var) const
  {
    if (den_.size() == 1 && den_[0] == 1)
      return num_.to_str(var);
    return format("({})/({})",num_.to_str(var),den_.to_str(var));
  }


  UInt ThieleURationalFunction::eval(UInt x, Mod mod) const
  {
    int terms = size()-1;
    UInt res = ai_[terms];
    for (int i=terms-1; i>=0; --i) {
      res = div_mod(sub_mod(x, getXi(i), mod), res, mod);
      res = add_mod(res, ai_[i], mod);
    }
    return res;
  }


  UInt ThieleURationalFunction::safe_eval(UInt x, Mod mod) const
  {
    int terms = size()-1;
    UInt num = ai_[terms];
    UInt den = 1;
    for (int i=terms-1; i>=0; --i) {
      std::swap(num,den);
      num = mul_mod(num, sub_mod(x, getXi(i), mod), mod);
      addmul_mod(num, ai_[i], den, mod);
    }
    if (!den)
      return FAILED;
    return div_mod(num, den, mod);
  }


  UInt ThieleURationalFunction::evalThieleSystemCoeff(UInt x, UInt y,
                                                     Mod mod,
                                                     unsigned i) const
  {
    if (i==0)
      return y;

    UInt num = sub_mod(x, getX0(), mod);
    UInt den = sub_mod(y, ai_[0], mod);

    for (unsigned j=1; j<=i-1; ++j) {
      addmul_mod(num, ai_[j], neg_mod(den, mod), mod);
      if (!num)
        return FAILED;
      std::swap(num,den);
      num = mul_mod(num, sub_mod(x, getXi(j), mod), mod);
    }

    if (!den)
      return FAILED;

    return div_mod(num, den, mod);
  }


  void ThieleURationalFunction::toStdRep(Mod mod, URationalFunction & f) const
  {
    unsigned numdeg = getNumDegree();
    unsigned dendeg = getDenDegree();

    // this accounts for the number of times the algorithm swaps num
    // and den
    if (!(size() % 2))
      std::swap(numdeg,dendeg);

    f.clear(); // needed?
    f.setNumDegree(numdeg);
    f.setDenDegree(dendeg);

    int maxterm = size()-1;

    f.num(numdeg) = ai_[maxterm];
    f.den(dendeg) = 1;

    // these count the current degree of num,den
    unsigned dn = 0, dd = 0;

    for (int i=maxterm-1; i>=0; --i) {
      f.invert();
      std::swap(dn,dd);
      std::swap(numdeg,dendeg);
      for (unsigned j = numdeg-dn-1; j<numdeg; ++j)
        addmul_mod(f.num(j), neg_mod(getXi(i), mod), f.num(j+1), mod);
      ++dn;
      for (unsigned j=0; j<=dd; ++j)
        addmul_mod(f.num(numdeg-dn+j), ai_[i], f.den(dendeg-dd+j), mod);
    }

    // not needed but a good idea
    f.normalize(mod);
  }


  void SparseRationalFunction::normalize(Mod mod)
  {
    // if den_ == 0, return 0/0
    if (den_.is_zero()) {
      num_.clear();
      return;
    }

    UInt norm = mul_inv(den_.lowest().coeff(), mod);

    for (auto & m : num_)
      m.coeff() = mul_mod(m.coeff(), norm, mod);

    for (auto & m : den_)
      m.coeff() = mul_mod(m.coeff(), norm, mod);

  }

  std::string SparseRationalFunction::to_str(const std::string vars[]) const
  {
    if (den_.size() == 1 && den_.degree() == 0 && den_.constTerm() == 1)
      return num_.to_str(vars);
    return format("({})/({})",num_.to_str(vars),den_.to_str(vars));
  }


  UInt horner_ratfun_eval(const UInt * __restrict numd,
                          const UInt * __restrict dend,
                          unsigned nvars,
                          UInt * __restrict workspace,
                          const UInt * __restrict x,
                          const UInt * __restrict xp,
                          Mod mod)
  {
    UInt num, den;

    if (!numd)
      return 0;

    if (!dend || !(den = horner_eval(dend, nvars, workspace, x, xp, mod)))
      return FAILED;

    num = horner_eval(numd, nvars, workspace, x, xp, mod);
    if (den == 1)
      return num;
    return div_mod(num, den, mod);
  }


  // Laurent expansion

  int laurent_expansion_learn(const UInt * num, unsigned num_deg,
                              const UInt * den, unsigned den_deg)
  {
    // factor out x^pref_exp
    int pref_exp_num=0, pref_exp_den=0;
    for (unsigned j=0; j<num_deg && (num[j]==0); ++j)
      ++pref_exp_num;
    for (unsigned j=0; j<den_deg && (den[j]==0); ++j)
      ++pref_exp_den;
    int pref_exp = pref_exp_num - pref_exp_den;
    return pref_exp;
  }

  void laurent_expansion(const UInt * num, int num_deg,
                         const UInt * den, int den_deg,
                         int order, Mod mod, int pref_exp,
                         UInt coeff[])
  {
    if (pref_exp > order)
      return;

    if (pref_exp > 0) {
      num += pref_exp;
      num_deg -= pref_exp;
    } else {
      den += -pref_exp;
      den_deg -= -pref_exp;
    }

    // note: we assume conventional normalization den[0] = 1

    coeff[0] = num[0];

    for (int i=1; i<=order-pref_exp; ++i) {

      UInt sub = 0;
      int jmax = std::min<int>(i, den_deg);

      for (int j=1; j<=jmax; ++j)
        sub = apbc_mod(sub, coeff[i-j], den[j], mod);

      coeff[i] = i<=num_deg ? sub_mod(num[i], sub, mod) : neg_mod(sub, mod);

    }
  }


} // namespace fflow
