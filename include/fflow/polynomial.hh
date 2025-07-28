#ifndef FFLOW_POLYNOMIAL_HH
#define FFLOW_POLYNOMIAL_HH

#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>
#include <fflow/common.hh>
#include <fflow/gcd.hh>
#include <fflow/integer_math.hh>

namespace fflow {


  ////////////////////////////
  // Univariate polynomials //
  ////////////////////////////

  // Univariate polynomial in a standard representation
  class UPoly {
  public:

    explicit UPoly(std::size_t degree=0) : coeff_(degree+1) {}

    void reserveDegree(std::size_t degree)
    {
      coeff_.reserve(degree+1);
    }

    void setDegree(std::size_t degree)
    {
      coeff_.resize(degree+1);
    }

    std::size_t getDegree() const
    {
      return coeff_.size()-1;
    }

    std::size_t getMinDegree() const
    {
      for (unsigned i=0; i<coeff_.size(); ++i)
        if (coeff_[i])
          return i;
      return 0;
    }

    std::size_t size() const
    {
      return coeff_.size();
    }

    // reset to 0
    void clear()
    {
      setDegree(0);
      coeff_[0] = 0;
    }

    // reset to zero but with degree=deg
    void clear(std::size_t deg)
    {
      setDegree(deg);
      for (unsigned i=0; i<=deg; ++i)
        coeff_[i] = 0;
    }

    UInt operator() (UInt x, Mod mod) const
    {
      return eval(x,mod);
    }

    const UInt * coeffs() const
    {
      return coeff_.data();
    }

    UInt * coeffs()
    {
      return coeff_.data();
    }

    UInt & operator[] (std::size_t i)
    {
      return coeff_[i];
    }

    UInt operator[] (std::size_t i) const
    {
      return coeff_[i];
    }

    // leading term
    UInt lt() const
    {
      return coeff_.back();
    }

    // remove higher degree terms which are zeroes
    void remove_high_zeroes()
    {
      while (coeff_.size()>1 && coeff_.back() == 0)
        coeff_.pop_back();
    }

    // note: must call remove_high_zeroes before calling this
    bool is_equal(const UPoly & other) const;

    // remove entries of degree<n and divide by x^n.  Entries of
    // degree=n+1 and above are filled by zeroes (call
    // remove_high_zeroes to get rid of them)
    void shift_left(std::size_t n);

    // same as remove_high_zeroes, but also requesting a reduction of
    // the capacity of the internal vector
    void shrink_to_fit()
    {
      remove_high_zeroes();
      coeff_.shrink_to_fit();
    }

    UInt eval(UInt x, Mod mod) const;

    void terms(UInt x, Mod mod, UInt res[]) const;

    void swap(UPoly & other)
    {
      coeff_.swap(other.coeff_);
    }

    std::string to_str(const std::string & var) const;

  private:
    std::vector<UInt> coeff_;
  };


  //////////////////////////////
  // Multivariate polynomials //
  //////////////////////////////


  // class for a monomial, which stores a chunk of contigous data on
  // the heap
  class Monomial {
  public:

    typedef std::int16_t VarExponent;

    static_assert(sizeof(UInt) % sizeof(VarExponent) == 0,
                  "Error: integer types have incompatible size.");

    Monomial() : data_() {}

    explicit Monomial(std::size_t nvars)
      : data_(new UInt[Monomial::uint_size_(nvars)]()) {}

    template <typename T>
    Monomial(UInt c, const std::initializer_list<T> & exponents)
      : data_()
    {
      allocate(c,exponents);
    }

    Monomial(const Monomial & oth, std::size_t nvars)
      : data_()
    {
      const std::size_t size = Monomial::uint_size_(nvars);
      data_.reset(new UInt[size]);
      from_same_vars_size_(oth,size);
    }

    // a monomial doesn't store the number of variables, so there
    // would be no way to copy the data in a copy constructor (use the
    // previous one instead)
    Monomial(const Monomial & oth) = delete;

    // move constructors are instead okay (nvars is there for
    // compatibility with trhe previous one, but it's ignored)
    Monomial(Monomial && oth, std::size_t nvars=0) noexcept
      : data_(std::move(oth.data_))
    {
      (void)(nvars);
    }

    // copy another monomial assuming they are both already allocated
    // with the same number nvars of variables
    void from_same_vars(const Monomial & oth, std::size_t nvars)
    {
      const std::size_t size = Monomial::uint_size_(nvars);
      from_same_vars_size_(oth,size);
    }

    void swap(Monomial & oth)
    {
      data_.swap(oth.data_);
    }

    // return a copy
    Monomial copy(std::size_t nvars) const
    {
      Monomial m(*this,nvars);
      return m;
    }

    // use copy instead of this (same rationale as for the deleted
    // copy constructor)
    Monomial & operator=(const Monomial & oth) = delete;

    // cannot be copied with operator= (see above) but can be moved
    Monomial & operator=(Monomial && oth) noexcept
    {
      data_ = std::move(oth.data_);
      return *this;
    }

    UInt coeff() const
    {
      return data_[0];
    }

    UInt & coeff()
    {
      return data_[0];
    }

    VarExponent degree() const
    {
      return *reinterpret_cast<VarExponent *>(data_.get() + 1);
    }

    VarExponent & degree()
    {
      return *reinterpret_cast<VarExponent *>(data_.get() + 1);
    }

    VarExponent exponent(std::size_t var) const
    {
      return reinterpret_cast<VarExponent *>(data_.get() + 1)[1+var];
    }

    VarExponent & exponent(std::size_t var)
    {
      return reinterpret_cast<VarExponent *>(data_.get() + 1)[1+var];
    }

    const VarExponent * exponents() const
    {
      return reinterpret_cast<const VarExponent *>(data_.get() + 1) + 1;
    }

    VarExponent * exponents()
    {
      return reinterpret_cast<VarExponent *>(data_.get() + 1) + 1;
    }

    template <typename T>
    void set(UInt c, const std::initializer_list<T> & exponents)
    {
      coeff() = c;
      std::size_t i = 0;
      unsigned deg = 0;
      for (const auto j : exponents) {
        exponent(i++) = j;
        deg += j;
      }
      degree() = deg;
    }

    template <typename T>
    void allocate(UInt c, const std::initializer_list<T> & exponents)
    {
      Monomial m(exponents.size());
      m.swap(*this);
      set(c,exponents);
    }

    operator bool() const
    {
      return data_.get();
    }

  private:

    // copy another monmial assuming they are both already allocated
    // data_ with the same size given as argument
    void from_same_vars_size_(const Monomial & oth, std::size_t size)
    {
      std::copy(oth.data_.get(),oth.data_.get()+size,data_.get());
    }

    static std::size_t uint_size_(std::size_t nvars)
    {
      const std::size_t max_padding = sizeof(UInt)/sizeof(VarExponent)-1;
      return (sizeof(UInt) + sizeof(VarExponent)*(max_padding+1+nvars)
              )/sizeof(UInt);
    }

    friend struct MonomialOp;

  private:
    // we store all the data (coeff,degree,exponent) in a cotiguous
    // array
    typedef std::unique_ptr<UInt[]> Ptr_;
    Ptr_ data_;
  };



  struct MonomialCompareLex {
    int operator() (const Monomial & m1, const Monomial & m2) const
    {
      for (unsigned i=0; i<nvars; ++i) {
        int cmp = compare(m1.exponent(i),m2.exponent(i));
        if (cmp)
          return cmp;
      }
      return 0;
    }
    std::size_t nvars;
  };

  struct MonomialCompare {
    int operator() (const Monomial & m1, const Monomial & m2) const
    {
      int cmp = compare(m1.degree(),m2.degree());
      return cmp ? cmp : MonomialCompareLex{nvars}(m1,m2);
    }
    std::size_t nvars;
  };

  struct MonomialIsLess {

    bool operator() (const Monomial & m1, const Monomial & m2) const
    {
      return MonomialCompare{nvars}(m1,m2) < 0;
    }

    std::size_t nvars;
  };


  // operations on monomials
  struct MonomialOp {

    Monomial copy(const Monomial & oth)
    {
      return oth.copy(n_vars);
    }

    int compare(const Monomial & m1, const Monomial & m2) const
    {
      MonomialCompare cmp{n_vars};
      return cmp(m1,m2);
    }

    // assumes all monomials are already allocated with same number of
    // vars

    void add(const Monomial & m1, const Monomial & m2, Mod mod,
             Monomial & res) const
    {
      res.from_same_vars(m1,n_vars);
      res.coeff() = add_mod(res.coeff(), m2.coeff(), mod);
    }

    void sub(const Monomial & m1, const Monomial & m2, Mod mod,
             Monomial & res) const
    {
      res.from_same_vars(m1,n_vars);
      res.coeff() = sub_mod(res.coeff(), m2.coeff(), mod);
    }

    void mul(const Monomial & m1, const Monomial & m2, Mod mod,
             Monomial & res) const
    {
      std::size_t size = Monomial::uint_size_(n_vars);
      for (unsigned i=1; i<size; ++i)
        res.data_[i] = m1.data_[i] + m2.data_[i];
      res.coeff() = mul_mod(m1.coeff(), m2.coeff(), mod);
    }

    void div(const Monomial & m1, const Monomial & m2, Mod mod,
             Monomial & res) const
    {
      std::size_t size = Monomial::uint_size_(n_vars);
      for (unsigned i=1; i<size; ++i)
        res.data_[i] = m1.data_[i] - m2.data_[i];
      res.coeff() = div_mod(m1.coeff(), m2.coeff(), mod);
    }

    Monomial add(const Monomial & m1, const Monomial & m2, Mod mod) const
    {
      Monomial res(n_vars);
      add(m1,m2,mod,res);
      return res;
    }

    Monomial sub(const Monomial & m1, const Monomial & m2, Mod mod) const
    {
      Monomial res(n_vars);
      sub(m1,m2,mod,res);
      return res;
    }

    Monomial mul(const Monomial & m1, const Monomial & m2, Mod mod) const
    {
      Monomial res(n_vars);
      mul(m1,m2,mod,res);
      return res;
    }

    Monomial div(const Monomial & m1, const Monomial & m2, Mod mod) const
    {
      Monomial res(n_vars);
      div(m1,m2,mod,res);
      return res;
    }

    Monomial copy(const Monomial & m) const
    {
      Monomial ret(m.copy(n_vars));
      return ret;
    }

    UInt eval(const Monomial & m, const UInt x[], Mod mod) const
    {
      const Monomial::VarExponent * exp = m.exponents();
      UInt res = m.coeff();
      for (unsigned i=0; i<n_vars; ++i)
        res = mul_mod(res, power(x[i],exp[i],mod), mod);
      return res;
    }

    void print(MemoryWriter & o, const Monomial & m,
               const std::string vars[],
               bool print_coeff = true) const;

    // data: only n_vars
    std::size_t n_vars;
  };


  // sparse multivariate polynomial.  ordered from higher to lower in
  // graded lexicographic
  class SparsePoly {
  public:

    typedef std::vector<Monomial>::iterator iterator;
    typedef std::vector<Monomial>::const_iterator const_iterator;

    SparsePoly() : monomials_(), op_{0} {}

    explicit SparsePoly(std::size_t nvars) : monomials_(), op_{nvars} {}

    SparsePoly(const SparsePoly & oth) : monomials_(), op_{0}
    {
      copy(oth);
    }

    SparsePoly(SparsePoly && oth) noexcept : monomials_(), op_{0}
    {
      copy(std::move(oth));
    }

    SparsePoly & operator=(const SparsePoly & oth)
    {
      SparsePoly newpoly(oth);
      newpoly.swap(*this);
      return *this;
    }

    SparsePoly & operator=(SparsePoly && oth) noexcept
    {
      oth.swap(*this);
      return *this;
    }

    void swap(SparsePoly & oth)
    {
      std::swap(op_,oth.op_);
      monomials_.swap(oth.monomials_);
    }

    void clear()
    {
      monomials_.clear();
    }

    std::size_t nvars() const
    {
      return op_.n_vars;
    }

    std::size_t size() const
    {
      return monomials_.size();
    }

    bool is_zero() const
    {
      return monomials_.empty();
    }

    void copy(const SparsePoly & oth)
    {
      op_ = oth.op_;
      monomials_.resize(oth.monomials_.size());

      std::size_t i = 0;
      for (const auto & m : oth.monomials_)
        monomials_[i++] = op_.copy(m);
    }

    void copy(SparsePoly && oth)
    {
      op_ = oth.op_;
      monomials_ = std::move(oth.monomials_);
    }

    // build from a univariate polynomial in the variable var
    void fromUnivariate(const UPoly & up, std::size_t var = 0);

    // constant polynomial = c
    void fromConst(UInt c);

    // from poly = x_var + c
    void fromUnivariateLinear(std::size_t var, UInt c);

    // make polynomial of total degree "deg" and all coefficients
    // equal to 1 (useful for mapping a position to a list of
    // exponents, which is the inverse of PolyTerms::pos)
    void fromDegree(std::size_t deg);

    void fromMonomials(const Monomial * mons, unsigned nmons,
                       unsigned first_var=0);

    // multiply by constant
    void mul(UInt c, Mod mod);

    // multiply by monomial
    void mul(const Monomial & m, Mod mod);

    // multiply by (x_var + c)
    void mulUnivariateLinear(std::size_t var, UInt c, Mod mod);

    // make the polynomial homogeneous of degree degree by multiplying
    // by suitable powers of the variables var, assuming the input is
    // independent of var
    void homogenize(std::size_t var, std::size_t degree);

    // add/sub other polynomials
    void add(const SparsePoly & oth, Mod mod);
    void sub(const SparsePoly & oth, Mod mod);
    void add(SparsePoly && oth, Mod mod);
    void sub(SparsePoly && oth, Mod mod);

    // add a monomial
    template <typename T>
    void add(UInt c, const std::initializer_list<T> & exponents, Mod mod)
    {
      SparsePoly p(nvars());
      p.monomials_.push_back(Monomial(c,exponents));
      add(std::move(p),mod);
    }

    // multiply by -1
    void neg(Mod mod)
    {
      for (auto & m : monomials_)
        m.coeff() = neg_mod(m.coeff(), mod);
    }

    // p(x) = p(x+c)
    void shift(const UInt c[], Mod mod);

    void shift(const std::initializer_list<UInt> & c, Mod mod)
    {
      std::vector<UInt> v(c);
      shift(v.data(),mod);
    }

    // return p(x+c) - p(x)
    SparsePoly getShiftedTerms(const UInt c[], Mod mod) const;
    SparsePoly getShiftedTerms(const std::initializer_list<UInt> & c,
                               Mod mod) const
    {
      std::vector<UInt> v(c);
      return getShiftedTerms(v.data(),mod);
    }

    UInt eval(const UInt x[], Mod mod) const;

    // evaluate only terms of total degree deg
    UInt evalDegree(unsigned deg, const UInt x[], Mod mod) const;

    UInt constTerm() const
    {
      if (monomials_.empty())
        return 0;

      if (monomials_.back().degree() == 0)
        return monomials_.back().coeff();
      else
        return 0;
    }

    unsigned degree() const
    {
      if (monomials_.empty())
        return 0;
      return monomials_[0].degree();
    }

    MonomialOp op() const
    {
      return op_;
    }


    iterator begin()
    {
      return monomials_.begin();
    }
    iterator end()
    {
      return monomials_.end();
    }
    const_iterator begin() const
    {
      return monomials_.begin();
    }
    const_iterator end() const
    {
      return monomials_.end();
    }
    const_iterator cbegin() const
    {
      return monomials_.cbegin();
    }
    const_iterator cend() const
    {
      return monomials_.cend();
    }

    // access array of monomials
    const Monomial * data() const
    {
      return monomials_.data();
    }
    Monomial * data()
    {
      return monomials_.data();
    }

    const Monomial & monomial(std::size_t i) const
    {
      return monomials_[i];
    }

    // undefined if *this.is_zero()
    const Monomial & lowest()
    {
      return monomials_.back();
    }

    // if you set by hand some monomial coefficients to zero, then
    // call this
    void remove_zeroes()
    {
      remove_zeroes_();
    }

    std::string to_str(const std::string vars[]) const;


  private:

    // static utility methods distinguishing btw. lvalues and rvalues

    // This copies the i-th monomial if p is an lvalue, or return an
    // already existing pointer to that if p is an rvalue.  Notice
    // that for rvalues we move the original data, which sets the
    // monomial to nullptr, therefore all methods calling this should
    // properly clear p before returning, otherwise p is left in an
    // invalid state.  For this purpose one can use clear_if_rv_.
    static Monomial monomial_copy_(const SparsePoly & p, std::size_t i)
    {
      return p.op_.copy(p.monomials_[i]);
    }
    static Monomial monomial_copy_(SparsePoly && p, std::size_t i)
    {
      return std::move(p.monomials_[i]);
    }

    static void clear_if_rv_(SparsePoly && p)
    {
      p.clear();
    }
    static void clear_if_rv_(const SparsePoly &)
    {
    }

  private:

    // utility for removing monomials with zero coefficients, which
    // indeed should never be there (call this before returning from a
    // method which might introduce zeroes)
    void remove_zeroes_();

    // utility for sorting the monomials, assuming they are all
    // distinct
    void sort_()
    {
      std::sort(monomials_.begin(), monomials_.end(),
                [this](const Monomial & m1, const Monomial & m2) -> bool
                {
                  return op_.compare(m1,m2)>0;
                });
    }

    template <typename SparsePolyT, typename AddOP>
    void add_(SparsePolyT && oth, Mod mod);

  private:
    friend class MPReconstructedPoly;

  private:
    std::vector<Monomial> monomials_;
    MonomialOp op_;
  };



  inline void SparsePoly::fromUnivariateLinear(std::size_t var, UInt c)
  {
    clear();
    monomials_.reserve(2);
    Monomial m0(nvars());
    Monomial m1(nvars());
    m0.coeff() = c;
    m0.degree() = 0;
    m1.coeff() = 1;
    m1.degree() = 1;
    m1.exponent(var) = 1;
    monomials_.push_back(std::move(m0));
    monomials_.push_back(std::move(m1));
  }

  inline void SparsePoly::fromConst(UInt c)
  {
    clear();
    if (c) {
      Monomial m0(nvars());
      m0.coeff() = c;
      m0.degree() = 0;
      monomials_.push_back(std::move(m0));
    }
  }

  inline void SparsePoly::mul(UInt c, Mod mod)
  {
    if (!c) {
      clear();
      return;
    }
    for (auto & m : monomials_)
      m.coeff() = mul_mod(m.coeff(), c, mod);
    // note: if mod i prime, no zeroes are introduced
  }

  inline void SparsePoly::mul(const Monomial & moth, Mod mod)
  {
    if (!moth.coeff()) {
      clear();
      return;
    }
    for (auto & m : monomials_)
      op_.mul(m, moth, mod, m);
  }

  inline void SparsePoly::fromMonomials(const Monomial * mons,
                                        unsigned nmons,
                                        unsigned first_var)
  {
    monomials_.resize(nmons);

    if (first_var == 0) {
      for (unsigned i=0; i<nmons; ++i)
        monomials_[i] = op_.copy(mons[i]);
    } else {
      for (unsigned i=0; i<nmons; ++i) {
        Monomial & m = monomials_[i] = Monomial(nvars());
        const Monomial & n = mons[i];
        m.coeff() = n.coeff();
        m.degree() = n.degree();
        std::copy(n.exponents(), n.exponents()+nvars()-first_var,
                  m.exponents()+first_var);
      }
    }

    remove_zeroes_();
    sort_();
  }


  inline
  SparsePoly SparsePoly::getShiftedTerms(const UInt c[], Mod mod) const
  {
    SparsePoly oth(*this);
    oth.shift(c,mod);
    oth.sub(*this,mod);
    return oth;
  }



  ////////////////////////
  // Newton polynomials //
  ////////////////////////


  // Generic abstract interface to Newton polynomials
  class NewtonPoly {
  public:

    typedef std::unique_ptr<NewtonPoly> Ptr;

    NewtonPoly() : xi_() {}

    virtual ~NewtonPoly() {}

    virtual void setDegree(std::size_t degree) = 0;
    virtual std::size_t getDegree() const = 0;
    virtual std::size_t size() const = 0;

    virtual unsigned nvars() const = 0;

    // returns a newly-allocated and unique copy, wrapped by a
    // uniqe_ptr
    virtual NewtonPoly::Ptr copy() const = 0;

    virtual bool is_univariate() const
    {
      return false;
    }

    virtual bool is_zero() const = 0;

    void setX0(UInt x0)
    {
      xi_[0] = x0;
    }

    std::size_t getX0() const
    {
      return xi_[0];
    }

    UInt getXi(std::size_t i) const
    {
      return xi_[i];
    }

    void setXi(std::size_t i, UInt x)
    {
      xi_[i] = x;
    }

    // set x0 for all variables
    virtual void setAllX0(UInt x0) = 0;

    // polynomial coefficients
    virtual const NewtonPoly::Ptr * pcoeffs() const = 0;
    virtual NewtonPoly::Ptr * pcoeffs() = 0;

    // reset to 0
    virtual void clear() = 0;

    // remove higher degree terms which are zeroes
    virtual void remove_high_zeroes() = 0;

    // same as is_zero(), but only guaranteed to work after
    // remove_high_zeroes()
    virtual bool is_zero_quick() const = 0;

    // same as remove_high_zeroes, but also requesting a reduction of
    // the capacity of the internal vector
    virtual void shrink_to_fit() = 0;

    // evaluate the polynomial
    virtual UInt eval(UInt x, UInt xi[], Mod mod) const = 0;

    void toSparsePoly(Mod mod, SparsePoly & res) const
    {
      toSparsePolyVar(0,mod,res);
    }

    // Convert to sparse polynomial in nvars variables assuming
    // first_var is the first.  This is called recursively, ending
    // with the univariate case.
    virtual void toSparsePolyVar(unsigned first_var, Mod mod,
                                 SparsePoly & res) const = 0;

    static void writeSparsePolyVar(const NewtonPoly & poly,
                                   unsigned first_var, Mod mod,
                                   SparsePoly & res)
    {
      poly.toSparsePolyVar(first_var,mod,res);
    }

    static NewtonPoly::Ptr create(std::size_t nvars, UInt x0 = 1);

  protected:
    explicit NewtonPoly(std::size_t degree) : xi_(degree+1) {}
    explicit NewtonPoly(std::size_t degree, UInt val) : xi_(degree+1, val) {}
    NewtonPoly(const NewtonPoly & oth) : xi_(oth.xi_) {}
    NewtonPoly(NewtonPoly && oth) : xi_(std::move(oth.xi_)) {}

  protected:
    std::vector<UInt> xi_;
  };


  // this is thrown if a multivariate method is called for the
  // univariate implementation case
#ifdef FFLOW_USE_EXCEPTIONS
  struct UnivariateException: public std::runtime_error {
    UnivariateException()
      : std::runtime_error("Calling multivariate method "
                           "on univariate polynomial") {}
  };
#endif
  inline void univariate_error()
  {
    FFLOW_THROW(UnivariateException());
  }


  // Univariate Newton polynomial with spacing = 1
  class NewtonUPoly final : public NewtonPoly {
  public:
    NewtonUPoly(std::size_t degree=0, UInt x0=1) : NewtonPoly(degree, x0), ai_(degree+1) {}
    NewtonUPoly(const NewtonUPoly & oth) : NewtonPoly(oth), ai_(oth.ai_) {}
    NewtonUPoly(NewtonUPoly && oth) : NewtonPoly(std::move(oth)), ai_(std::move(oth.ai_)) {}

    void swap(NewtonUPoly & oth)
    {
      xi_.swap(oth.xi_);
      ai_.swap(oth.ai_);
    }

    NewtonUPoly & operator=(const NewtonUPoly & oth)
    {
      NewtonUPoly p(oth);
      p.swap(*this);
      return *this;
    }

    NewtonUPoly & operator=(NewtonUPoly && oth)
    {
      oth.swap(*this);
      return *this;
    }

    virtual void setDegree(std::size_t degree) override final
    {
      xi_.resize(degree+1);
      ai_.resize(degree+1);
    }

    virtual std::size_t getDegree() const override final
    {
      return ai_.size()-1;
    }

    virtual std::size_t size() const override final
    {
      return ai_.size();
    }

    virtual unsigned nvars() const override final
    {
      return 1;
    }

    void push_back(UInt xi, UInt a)
    {
      xi_.push_back(xi);
      ai_.push_back(a);
    }

    virtual NewtonPoly::Ptr copy() const override final
    {
      return NewtonPoly::Ptr(new NewtonUPoly(*this));
    }

    virtual bool is_univariate() const override final
    {
      return true;
    }

    virtual bool is_zero() const override final
    {
      return std::all_of(ai_.begin(),ai_.end(),[](UInt i) { return i==0; });
    }

    UInt operator() (UInt x, Mod mod) const
    {
      return eval(x,mod);
    }

    virtual const NewtonPoly::Ptr * pcoeffs() const override final
    {
      univariate_error();
      return nullptr;
    }
    virtual NewtonPoly::Ptr * pcoeffs() override final
    {
      univariate_error();
      return nullptr;
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
    virtual void clear() override final
    {
      setDegree(0);
      ai_[0] = 0;
    }

    // reset to zero but with degree=deg
    virtual void clear(std::size_t deg) //override final
    {
      setDegree(deg);
      for (unsigned i=0; i<=deg; ++i)
        ai_[i] = 0;
    }

    // remove higher degree terms which are zeroes
    virtual void remove_high_zeroes() override final
    {
      while (ai_.size()>1 && ai_.back() == 0)
        ai_.pop_back();
    }

    virtual bool is_zero_quick() const override final
    {
      return (ai_.size() == 1) && (ai_[0] == 0);
    }

    // same as remove_high_zeroes, but also requesting a reduction of
    // the capacity of the internal vector
    virtual void shrink_to_fit() override final
    {
      remove_high_zeroes();
      ai_.shrink_to_fit();
    }

    virtual UInt eval(UInt x, UInt *, Mod mod) const override final
    {
      return eval(x,mod);
    }

    UInt eval(UInt x, Mod mod) const;

    virtual void setAllX0(UInt x0) override final
    {
      setX0(x0);
    }

    // utility for newton interpolation system
    //
    // returns = \sum_i^{row-1} a_i eta_i and writes eta_row in last
    // argument
    UInt evalNewtonSystemRow(UInt x, Mod mod, unsigned row,
                             UInt & eta_row) const;

    // write the coefficient of a std rep. of the same polynomial.
    // mod must be consistent with the current entries of the
    // instance.
    void toStdRep(Mod mod, UInt coeffs[]) const;

    virtual void toSparsePolyVar(unsigned first_var, Mod mod,
                                 SparsePoly & res) const override final;

  private:
    using NewtonPoly::xi_;
    std::vector<UInt> ai_;
  };


  // Multivariate Newton polynomial with spacing = 1.  Defined
  // recursively as a vector of (unique) pointers to NewtonPoly,
  // terminating the recursion to the univariate case
  class NewtonMPoly final : public NewtonPoly {
  public:
    NewtonMPoly(unsigned nvars, std::size_t degree=0, UInt x0=1)
      : NewtonPoly(degree, x0), ai_()
    {
      ai_.reserve(degree+1);
      for (unsigned i=0; i<degree+1; ++i)
        ai_.push_back(new_entry_(nvars-1));
    }

    NewtonMPoly(const NewtonMPoly & oth)
      : NewtonPoly(oth), ai_()
    {
      ai_.reserve(oth.ai_.size());
      for (const auto & p : oth.ai_)
        ai_.push_back((*p).copy());
    }

    NewtonMPoly(NewtonMPoly && oth)
      : NewtonPoly(std::move(oth)), ai_(std::move(oth.ai_)) {}

    void swap(NewtonMPoly & oth)
    {
      xi_.swap(oth.xi_);
      ai_.swap(oth.ai_);
    }

    NewtonMPoly & operator=(const NewtonMPoly & oth)
    {
      NewtonMPoly p(oth);
      p.swap(*this);
      return *this;
    }

    NewtonMPoly & operator=(NewtonMPoly && oth)
    {
      oth.swap(*this);
      return *this;
    }

    virtual void setDegree(std::size_t degree) override final
    {
      unsigned old_degree = getDegree();
      unsigned n = nvars();
      xi_.resize(degree+1);
      ai_.resize(degree+1);
      for (unsigned i=old_degree+1; i<=getDegree(); ++i)
        ai_[i] = new_entry_(n-1);
    }

    virtual std::size_t getDegree() const override final
    {
      return ai_.size()-1;
    }

    virtual std::size_t size() const override final
    {
      return ai_.size();
    }

    virtual unsigned nvars() const override final;

    void push_back(UInt xi, const NewtonPoly & ai)
    {
      xi_.push_back(xi);
      ai_.push_back(ai.copy());
    }

    void push_back(UInt xi, NewtonPoly::Ptr && ai)
    {
      xi_.push_back(xi);
      ai_.push_back(std::move(ai));
    }

    virtual bool is_zero() const override final
    {
      return std::all_of(ai_.begin(),ai_.end(),
                         [](const NewtonPoly::Ptr & p)
                         { return p->is_zero(); });
    }

    virtual NewtonPoly::Ptr copy() const override final
    {
      return NewtonPoly::Ptr(new NewtonMPoly(*this));
    }

    // TODO
    //UInt operator() (UInt x, Mod mod) const;

    virtual const NewtonPoly::Ptr * pcoeffs() const override final
    {
      return ai_.data();
    }

    virtual NewtonPoly::Ptr * pcoeffs() override final
    {
      return ai_.data();
    }

    // reset to 0
    virtual void clear() override final
    {
      setDegree(0);
      ai_[0] = new_entry_(nvars()-1);
    }

    // remove higher degree terms which are zeroes
    virtual void remove_high_zeroes() override final
    {
      while (ai_.size()>1 && (*ai_.back()).is_zero())
        ai_.pop_back();
    }

    virtual bool is_zero_quick() const override final
    {
      return (ai_.size() == 1) && ((*ai_[0]).is_zero_quick());
    }

    void remove_high_zeroes_quick()
    {
      while (ai_.size()>1 && (*ai_.back()).is_zero_quick())
        ai_.pop_back();
    }

    // same as remove_high_zeroes, but also requesting a reduction of
    // the capacity of the internal vector
    virtual void shrink_to_fit() override final
    {
      remove_high_zeroes();
      ai_.shrink_to_fit();
    }

    virtual UInt eval(UInt x, UInt xi[], Mod mod) const override final;

    virtual void setAllX0(UInt x0) override final
    {
      setX0(x0);
      for (auto & p : ai_)
        (*p).setAllX0(x0);
    }

    virtual void toSparsePolyVar(unsigned first_var, Mod mod,
                                 SparsePoly & res) const override final;

  private:

    // creates a new zero entry
    NewtonPoly::Ptr new_entry_(unsigned nvars_minus_1) const
    {
      if (nvars_minus_1 == 1)
        return NewtonPoly::Ptr(new NewtonUPoly(0,xi_[0]));
      else
        return NewtonPoly::Ptr(new NewtonMPoly(nvars_minus_1,0,xi_[0]));
    }

  private:
    using NewtonPoly::xi_;
    std::vector<NewtonPoly::Ptr> ai_;
  };


  inline NewtonPoly::Ptr NewtonPoly::create(std::size_t nvars, UInt x0)
  {
    if (nvars == 1)
      return NewtonPoly::Ptr(new NewtonUPoly(0,x0));
    else
      return NewtonPoly::Ptr(new NewtonMPoly(nvars,0,x0));
  }


  // Horner representation

  const UInt HORNER_END = ~UInt(0);

  UInt horner_eval(const UInt * __restrict horner_data, unsigned nvars,
                   UInt * __restrict workspace,
                   const UInt * __restrict x,
                   const UInt * __restrict xp_shoup,
                   Mod mod);

  std::size_t horner_from_sparse_poly(UInt * horner_data,
                                      const Monomial * m,
                                      std::size_t size,
                                      unsigned nvars,
                                      unsigned first_var=0,
                                      unsigned * map_pos = 0);

  std::size_t horner_required_workspace(const UInt * horner_data);

  std::size_t horner_size(const UInt * __restrict horner_data, unsigned nvars);


  typedef std::unique_ptr<UInt[]> HornerPtr;

  HornerPtr horner_clone(const UInt * __restrict horner_data, unsigned nvars);


  class HornerWorkspacePtr {
  public:
    HornerWorkspacePtr() : ptr_(nullptr), size_(0) {}

    HornerWorkspacePtr(const HornerWorkspacePtr & oth) = delete;

    HornerWorkspacePtr(HornerWorkspacePtr && oth) : ptr_(nullptr), size_(0)
    {
      swap(oth);
    }

    void swap(HornerWorkspacePtr & oth)
    {
      std::swap(ptr_, oth.ptr_);
      std::swap(size_, oth.size_);
    }

    void resize(std::size_t size)
    {
      size_ = size;
      ptr_ = static_cast<UInt *>(crealloc(static_cast<void*>(ptr_),
                                          sizeof(UInt)*size));
    }

    std::size_t size() const
    {
      return size_;
    }

    void ensure_size(std::size_t size)
    {
      if (size > size_)
        resize(size);
    }

    ~HornerWorkspacePtr()
    {
      std::free(ptr_);
    }

    const UInt * get() const
    {
      return ptr_;
    }

    UInt * get()
    {
      return ptr_;
    }

  private:
    UInt * ptr_;
    std::size_t size_;
  };

  HornerPtr hornerptr_from_sparse_poly(const Monomial * m,
                                       std::size_t size,
                                       unsigned nvars,
                                       unsigned first_var=0,
                                       unsigned * map_pos = 0);


  // Utility for generating a list of variable names following a
  // pattern.  E.g. vars_with_pattern("x{}",3) will generate
  // {"x0","x1","x2"}.
  inline std::unique_ptr<std::string[]> vars_with_pattern(const char * patt,
                                                          unsigned nvars,
                                                          unsigned start = 0)
  {
    std::unique_ptr<std::string[]> res(new std::string[nvars]);
    for (unsigned i=0; i<nvars; ++i)
      res[i] = format(patt,start+i);
    return res;
  }

} // namespace fflow


#endif // FFLOW_POLYNOMIAL_HH
