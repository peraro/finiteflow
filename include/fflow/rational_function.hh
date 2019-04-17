#ifndef FFLOW_RATIONAL_FUNCTION_HH
#define FFLOW_RATIONAL_FUNCTION_HH

#include <fflow/common.hh>
#include <fflow/gcd.hh>
#include <fflow/polynomial.hh>

namespace fflow {

  ///////////////////////////////////
  // Univariate rational functions //
  ///////////////////////////////////

  // Univariate rational function
  class URationalFunction {
  public:
    explicit URationalFunction(std::size_t degree=0)
      : num_(degree), den_(degree) {}

    URationalFunction(std::size_t num_degree, std::size_t den_degree)
      : num_(num_degree), den_(den_degree) {}

    void setNumDegree(std::size_t degree)
    {
      num_.setDegree(degree);
    }

    void setDenDegree(std::size_t degree)
    {
      den_.setDegree(degree);
    }

    void setDegree(std::size_t degree)
    {
      setNumDegree(degree);
      setDenDegree(degree);
    }

    void setDegree(std::size_t num_degree, std::size_t den_degree)
    {
      setNumDegree(num_degree);
      setDenDegree(den_degree);
    }

    std::size_t getNumDegree() const
    {
      return num_.getDegree();
    }

    std::size_t getDenDegree() const
    {
      return den_.getDegree();
    }

    UPoly & num()
    {
      return num_;
    }

    UPoly & den()
    {
      return den_;
    }

    const UPoly & num() const
    {
      return num_;
    }

    const UPoly & den() const
    {
      return den_;
    }

    UInt & num(std::size_t i)
    {
      return num_[i];
    }

    UInt & den(std::size_t i)
    {
      return den_[i];
    }

    UInt num(std::size_t i) const
    {
      return num_[i];
    }

    UInt den(std::size_t i) const
    {
      return den_[i];
    }

    UInt operator() (UInt x, Mod mod) const
    {
      return eval(x,mod);
    }

    UInt eval(UInt x, Mod mod) const
    {
      return div_mod(num_.eval(x,mod), den_.eval(x,mod), mod);
    }

    // remove higher degree terms which are zeroes
    void remove_high_zeroes()
    {
      num_.remove_high_zeroes();
      den_.remove_high_zeroes();
    }

    // same as remove_high_zeroes, but also requesting a reduction of
    // the capacity of the internal vectors
    void shrink_to_fit()
    {
      num_.shrink_to_fit();
      den_.shrink_to_fit();
    }

    // swaps num and den
    void invert()
    {
      num_.swap(den_);
    }

    // reset to 0/0
    void clear()
    {
      setDegree(0);
      num_[0] = den_[0] = 0;
    }

    // divide all coefficients by the lowest degree non-zero entry of
    // the denominator
    void normalize(Mod mod);

    // only makes sense after a normalize()
    bool is_undefined() const
    {
      return (getNumDegree() == 0) && (getDenDegree() == 0) && (den_[0] == 0);
    }

    // only makes sense after a normalize()
    bool is_equal(const URationalFunction & other) const
    {
      return num_.is_equal(other.num_) && den_.is_equal(other.den_);
    }

    std::string to_str(const std::string & var) const;

  private:
    UPoly num_, den_;
  };


  // Univariate rational function as continued fraction (used for
  // Thiele's interpolation).  Numerator and denominator have same
  // degree.
  class ThieleURationalFunction {
  public:

    explicit ThieleURationalFunction(std::size_t terms=1)
      : ai_(terms ? terms : 1,0), xi_(terms ? terms : 1,1) {}

    void resize(std::size_t terms)
    {
      xi_.resize(terms);
      ai_.resize(terms);
    }

    std::size_t size() const
    {
      return ai_.size();
    }

    // get degree of num,den (possibly exceeding the latter by one)
    std::size_t getDegree() const
    {
      return ai_.size()/2;
    }

    // the following are more accurate than the previous one
    std::size_t getNumDegree() const
    {
      return ai_.size()/2;
    }
    std::size_t getDenDegree() const
    {
      return (ai_.size()-1)/2;
    }

    void setX0(UInt x0)
    {
      xi_[0] = x0;
    }

    UInt getX0() const
    {
      return xi_[0];
    }

    void setXi(std::size_t i, std::size_t xi)
    {
      xi_[i] = xi;
    }

    UInt getXi(std::size_t i) const
    {
      return xi_[i];
    }

    UInt operator() (UInt x, Mod mod) const
    {
      return eval(x,mod);
    }

    const UInt * coeffs() const
    {
      return ai_.data();
    }

    UInt * coeffs()
    {
      return ai_.data();
    }

    // reset to 0
    void clear()
    {
      resize(1);
      ai_[0] = 0;
    }

    // resize to n and reset entries to 0
    void clear(std::size_t n)
    {
      resize(n);
      for (auto & a : ai_)
        a = 0;
    }

    void push_back(UInt xi, UInt a)
    {
      xi_.push_back(xi);
      ai_.push_back(a);
    }

    // remove higher degree terms which are zeroes
    void remove_high_zeroes()
    {
      while (ai_.size()>1 && ai_.back() == 0)
        ai_.pop_back();
      xi_.resize(ai_.size());
    }

    // same as remove_high_zeroes, but also requesting a reduction of
    // the capacity of the internal vector
    void shrink_to_fit()
    {
      remove_high_zeroes();
      ai_.shrink_to_fit();
      xi_.shrink_to_fit();
    }

    UInt eval(UInt x, Mod mod) const;
    UInt safe_eval(UInt x, Mod mod) const;

    // evaluate the i-th coefficient a[i] assuming to know all the
    // previous coefficients and knowing the additional point f(x)=y.
    UInt evalThieleSystemCoeff(UInt x, UInt y, Mod mod, unsigned i) const;

    // writes a std rep. of the same function.  mod must be consistent
    // with the current entries of the instance.  It also sets the
    // right degrees for f
    void toStdRep(Mod mod, URationalFunction & f) const;

  private:
    std::vector<UInt> ai_;
    std::vector<UInt> xi_;
  };



  /////////////////////////////////////
  // Multivariate rational functions //
  /////////////////////////////////////

  class SparseRationalFunction {
  public:

    SparseRationalFunction() : num_(), den_() {}

    explicit SparseRationalFunction(std::size_t nvars, bool init_to_zero=false)
      : num_(nvars), den_(nvars)
    {
      if (init_to_zero)
        den_.fromConst(1);
    }

    SparseRationalFunction(const SparsePoly & num, const SparsePoly & den)
      : num_(num), den_(den) {}

    SparseRationalFunction(SparsePoly && num, SparsePoly && den)
      : num_(std::move(num)), den_(std::move(den)) {}

    SparseRationalFunction(const SparseRationalFunction & oth)
      : num_(), den_()
    {
      copy(oth);
    }

    SparseRationalFunction(SparseRationalFunction && oth)
      : num_(), den_()
    {
      copy(std::move(oth));
    }

    void swap(SparseRationalFunction & oth)
    {
      num_.swap(oth.num_);
      den_.swap(oth.den_);
    }

    SparseRationalFunction & operator=(const SparseRationalFunction & oth)
    {
      SparseRationalFunction newfun(oth);
      newfun.swap(*this);
      return *this;
    }

    SparseRationalFunction & operator=(SparseRationalFunction && oth)
    {
      oth.swap(*this);
      return *this;
    }

    const SparsePoly & numerator() const
    {
      return num_;
    }

    SparsePoly & numerator()
    {
      return num_;
    }

    const SparsePoly & denominator() const
    {
      return den_;
    }

    SparsePoly & denominator()
    {
      return den_;
    }

    void neg(Mod mod)
    {
      num_.neg(mod);
    }

    void copy(const SparseRationalFunction & oth)
    {
      num_.copy(oth.num_);
      den_.copy(oth.den_);
    }

    void copy(SparseRationalFunction && oth)
    {
      num_.copy(std::move(oth.num_));
      den_.copy(std::move(oth.den_));
    }

    UInt eval(const UInt x[], Mod mod) const
    {
      UInt num = num_.eval(x,mod);
      UInt den = den_.eval(x,mod);
      return div_mod(num, den, mod);
    }

    // build from a univariate function in the variable var
    void fromUnivariate(const URationalFunction & uf, std::size_t var = 0)
    {
      num_.fromUnivariate(uf.num(), var);
      den_.fromUnivariate(uf.den(), var);
    }

    std::string to_str(const std::string vars[]) const;

    // normalize such that lower degree term in denominator has
    // coefficient=1
    void normalize(Mod mod);

  private:
    SparsePoly num_, den_;
  };



  // Horner representation
  UInt horner_ratfun_eval(const UInt * __restrict num,
                          const UInt * __restrict den,
                          unsigned nvars,
                          UInt * __restrict workspace,
                          const UInt * __restrict x,
                          const UInt * __restrict xp_shoup,
                          Mod mod);

  // wrapper owning the pointers to num and den data
  class HornerRatFunPtr {
  public:
    HornerRatFunPtr() : num_(nullptr), den_(nullptr) {}

    HornerRatFunPtr(const HornerRatFunPtr & oth) = delete;

    HornerRatFunPtr(HornerRatFunPtr && oth)
      : num_(std::move(oth.num_)), den_(std::move(oth.den_)) {}

    HornerRatFunPtr & operator= (const HornerRatFunPtr & oth) = delete;

    HornerRatFunPtr & operator= (HornerRatFunPtr && oth)
    {
      num_ = std::move(oth.num_);
      den_ = std::move(oth.den_);
      return *this;
    }

    HornerPtr & num_ptr()
    {
      return num_;
    }

    HornerPtr & den_ptr()
    {
      return den_;
    }

    const UInt * num() const
    {
      return num_.get();
    }

    const UInt * den() const
    {
      return den_.get();
    }

    UInt eval(unsigned nvars, UInt * __restrict workspace,
              const UInt * __restrict x,
              const UInt * __restrict xp_shoup,
              Mod mod) const
    {
      return horner_ratfun_eval(num_.get(), den_.get(), nvars,
                                workspace, x, xp_shoup, mod);
    }

    void from_sparse_poly(const Monomial * mnum, std::size_t numsize,
                          const Monomial * mden, std::size_t densize,
                          unsigned nvars, unsigned first_var=0,
                          unsigned * num_pos = 0, unsigned * den_pos = 0)
    {
      num_ = hornerptr_from_sparse_poly(mnum, numsize, nvars,
                                        first_var, num_pos);
      den_ = hornerptr_from_sparse_poly(mden, densize, nvars,
                                        first_var, den_pos);
    }

    std::size_t required_workspace()
    {
      return std::max(horner_required_workspace(num()),
                      horner_required_workspace(den()));
    }

  private:
    HornerRatFunPtr(HornerRatFunPtr & oth) = delete;

  private:
    HornerPtr num_, den_;
  };


  int laurent_expansion_learn(const UInt * num, unsigned num_deg,
                              const UInt * den, unsigned den_deg);

  void laurent_expansion(const UInt * num, int num_deg,
                         const UInt * den, int den_deg,
                         int order, Mod mod, int pref_exp,
                         UInt coeff[]);

  inline int laurent_expansion_learn(const URationalFunction & ufun)
  {
    return laurent_expansion_learn(ufun.num().coeffs(), ufun.getNumDegree(),
                                   ufun.den().coeffs(), ufun.getDenDegree());
  }

  inline void laurent_expansion(const URationalFunction & ufun,
                                int order, Mod mod, int pref_exp,
                                UInt coeff[])
  {
    return laurent_expansion(ufun.num().coeffs(), ufun.getNumDegree(),
                             ufun.den().coeffs(), ufun.getDenDegree(),
                             order, mod, pref_exp, coeff);
  }

} // namespace fflow


#endif // FFLOW_RATIONAL_FUNCTION_HH
