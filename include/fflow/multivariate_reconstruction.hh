#ifndef FFLOW_MULTIVARIATE_RECONSTRUCTION_HH
#define FFLOW_MULTIVARIATE_RECONSTRUCTION_HH

#include <fflow/polynomial.hh>
#include <fflow/rational_function.hh>
#include <fflow/matrix.hh>
#include <fflow/primes.hh>
#include <fflow/univariate_reconstruction.hh>
#include <fflow/multivariate_reconstruction_details.hh>
#include <fflow/function_cache.hh>
#include <flint/flint.h>

namespace fflow {

  // defining interface which computes a polynomial or rational
  // function of x_i
  class RatFun {
  public:
    // must return FAILED if the values of x_i make the calculation
    // impossible/meaningless
    virtual UInt evaluate(const UInt x[], Mod mod) = 0;
    virtual ~RatFun() {}
  };


  // conversions

  class MtoURatFun : public URatFun {
  public:

    explicit MtoURatFun(RatFun & fun) : fun_(fun) {}

    virtual UInt evaluate(const UInt x, Mod mod)
    {
      return fun_.evaluate(&x, mod);
    }

  private:
    RatFun & fun_;
  };

  class UtoMRatFun : public RatFun {
  public:

    explicit UtoMRatFun(URatFun & fun) : fun_(fun) {}

    virtual UInt evaluate(const UInt x[], Mod mod)
    {
      return fun_.evaluate(*x, mod);
    }

  private:
    URatFun & fun_;
  };


  namespace detail {
    struct URatFunFreezed;
    struct URatFunFreezedSample;
  }


  // Reconstruct multivariate polynomials using Newton interpolation
  // method recursively
  class NPolyReconstruction {
  public:
    enum {DEFAULT_N_UCHECKS = 2};
    enum {DEFAULT_N_CHECKS = 2};
    enum {DEFAULT_N_SINGULAR = 10};
    enum {DEFAULT_MAX_DEG = 100};

    explicit NPolyReconstruction(unsigned nvars,
                                 unsigned maxdegree = DEFAULT_MAX_DEG)
      : poly_(NewtonPoly::create(nvars)),
        x0_(new UInt[nvars]()),
        nvars_(nvars), maxdeg_(maxdegree),
        n_checks(DEFAULT_N_CHECKS), n_uchecks(DEFAULT_N_UCHECKS),
        n_singular(DEFAULT_N_SINGULAR), xdegs(nullptr),
        check_over_maxdeg(true)
    {
      Mod mod(BIG_UINT_PRIMES[BIG_UINT_PRIMES_SIZE-1]);
      for (unsigned j=0; j<nvars; ++j)
        x0_[j] = sample_uint(OFFSET_1,j+1234567,mod);
    }

    void reset()
    {
      poly_ = NewtonPoly::create(nvars_);
    }

    void reset(unsigned nvars,
               unsigned maxdegree = DEFAULT_MAX_DEG)
    {
      nvars_ = nvars;
      maxdeg_ = maxdegree;
      x0_.reset(new UInt[nvars_]());
      for (unsigned j=0; j<nvars; ++j)
        x0_[j] = SAMPLING_STRIDE+(j+1)*64;
      reset();
    }

    unsigned getMaxDegree() const
    {
      return maxdeg_;
    }

    void setMaxDegree(unsigned deg)
    {
      maxdeg_ = deg;
    }

    void writeResult(Mod mod, SparsePoly & p)
    {
      (*poly_).toSparsePoly(mod,p);
    }

    void writeResultVar(std::size_t first_var, Mod mod, SparsePoly & p)
    {
      (*poly_).toSparsePolyVar(first_var,mod,p);
    }

    UInt getX0(std::size_t var = 0) const
    {
      return x0_[var];
    }

    UInt getXi(std::size_t i, std::size_t var = 0) const
    {
      return x0_[var] + i;
    }

    void setX0(UInt x0)
    {
      for (unsigned i=0; i<nvars_; ++i)
        x0_[i] = x0;
    }

    void setX0(const UInt x0[])
    {
      for (unsigned i=0; i<nvars_; ++i)
        x0_[i] = x0[i];
    }

    void setX0(const std::initializer_list<UInt> & x0)
    {
      std::size_t i=0;
      for (const auto t : x0)
        x0_[i++] = t;
    }

    Ret reconstruct(RatFun & f, Mod mod);
    void sample(RatFun & f, Mod mod);

    unsigned getTrueDegree();

    const NewtonMPoly & getNewtonMPoly() const
    {
      return static_cast<const NewtonMPoly &>(*poly_);
    }

  private:

    std::size_t maxdeg_var_(std::size_t i)
    {
      if (xdegs)
        return xdegs[i];
      return maxdeg_;
    }

    NewtonPoly::Ptr reconstruct_(detail::URatFunFreezed & f,
                                 unsigned nvars,
                                 unsigned fixed_vars,
                                 unsigned current_deg,
                                 Mod mod);

    void sample_(detail::URatFunFreezedSample & f,
                 unsigned nvars,
                 unsigned fixed_vars,
                 unsigned current_deg,
                 Mod mod);

  private:
    NewtonPoly::Ptr poly_;
    std::unique_ptr<UInt[]> x0_;
    std::size_t nvars_, maxdeg_;

  public:
    std::size_t n_checks, n_uchecks, n_singular;
    const std::size_t * xdegs;
    bool check_over_maxdeg;
  };


  class RatFunReconstruction {
  public:
    enum {DEFAULT_N_UNDCHECKS = 0};
    enum {DEFAULT_N_UCHECKS = 3};
    enum {DEFAULT_N_CHECKS = 2};
    enum {DEFAULT_N_SINGULAR = 3};
    enum {DEFAULT_MAX_DEG = 100};

    static bool VERBOSE;

    explicit RatFunReconstruction(unsigned nvars,
                                  unsigned maxdegree = DEFAULT_MAX_DEG)
      : fun_(nvars,true), npolyrec_(nvars > 1 ? nvars-1 : 1, maxdegree),
        tmp_(nvars), shift_(), xdegs_(), n0_(FAILED), t0_(OFFSET_0_HASH),
        nvars_(nvars), maxdeg_(maxdegree), num_pref_(), den_pref_(),
        n_checks(DEFAULT_N_CHECKS), n_uchecks(DEFAULT_N_UCHECKS),
        n_undchecks(DEFAULT_N_UNDCHECKS),
        n_singular(DEFAULT_N_SINGULAR) {}

    unsigned getMaxDegree() const
    {
      return maxdeg_;
    }

    void setMaxDegree(unsigned deg)
    {
      maxdeg_ = deg;
    }

    UInt getX0(std::size_t var = 0) const
    {
      return npolyrec_.getX0(var);
    }

    UInt getXi(std::size_t i, std::size_t var = 0) const
    {
      return npolyrec_.getXi(i,var);
    }

    void setX0(UInt x0)
    {
      npolyrec_.setX0(x0);
    }

    void setX0(const UInt x0[])
    {
      npolyrec_.setX0(x0);
    }

    void setX0(const std::initializer_list<UInt> & x0)
    {
      npolyrec_.setX0(x0);
    }

    UInt getT0() const
    {
      return t0_;
    }

    UInt getTi(std::size_t i) const
    {
      return t0_ + i;
    }

    void setT0(UInt t0)
    {
      t0_ = t0;
    }

    void setShift(const UInt c[])
    {
      if (!c) {
        shift_.clear();
      } else {
        shift_.resize(nvars_);
        for (unsigned i=0; i<nvars_; ++i)
          shift_[i] = c[i];
      }
    }

    void setShift(const std::initializer_list<UInt> & c)
    {
      shift_.clear();
      shift_.reserve(nvars_);
      for (const auto t : c)
        shift_.push_back(t);
    }

    void clearShift()
    {
      shift_.clear();
    }

    // compute the total degree
    Ret degree(RatFun & f, Mod mod, unsigned & numdeg, unsigned & dendeg)
    {
      if (nvars_ == 1)
        return degree_univariate_(f, mod, numdeg, dendeg);
      else
        return degree_(f, mod, numdeg, dendeg);
    }

    // compute min and max degrees wrt one variable
    Ret var_degree(RatFun & f, unsigned var, Mod mod,
                   std::size_t & numdeg_min, std::size_t & numdeg_max,
                   std::size_t & dendeg_min, std::size_t & dendeg_max);

    Ret reconstructWithDegs(RatFun & f, Mod mod,
                            unsigned numdeg, unsigned dendeg,
                            RatFunVarDegrees && degs)
    {
      if (nvars_ == 1)
        return reconstruct_univariate_(f, mod);

      xdegs_ = std::move(degs);
      get_monomial_prefactors_();

      return reconstruct_(f, mod,
                          numdeg - (num_pref_ ? num_pref_.degree() : 0),
                          dendeg - (den_pref_ ? den_pref_.degree() : 0));
    }

    Ret reconstructWithDegs(RatFun & f, Mod mod,
                            unsigned numdeg, unsigned dendeg,
                            const RatFunVarDegrees & degs)
    {
      RatFunVarDegrees newdegs;
      newdegs.copy(degs, nvars_);
      return reconstructWithDegs(f, mod, numdeg, dendeg, std::move(newdegs));
    }

    Ret reconstruct(RatFun & f, Mod mod, unsigned numdeg, unsigned dendeg)
    {
      if (nvars_ == 1)
        return reconstruct_univariate_(f, mod);

      Ret ret;
      if ((ret = vars_degrees_(f, mod)) != SUCCESS)
        return ret;
      get_monomial_prefactors_();

      return reconstruct_(f, mod,
                          numdeg - (num_pref_ ? num_pref_.degree() : 0),
                          dendeg - (den_pref_ ? den_pref_.degree() : 0));
    }

    Ret reconstruct(RatFun & f, Mod mod)
    {
      if (nvars_ == 1)
        return reconstruct_univariate_(f, mod);

      Ret ret;
      if ((ret = vars_degrees_(f, mod)) != SUCCESS)
        return ret;
      get_monomial_prefactors_();

      unsigned numdeg, dendeg;
      if ((ret = degree(f, mod, numdeg, dendeg)) != SUCCESS)
        return ret;
      return reconstruct_(f, mod, numdeg, dendeg);
    }

    void sample(RatFun & f, Mod mod,
                unsigned numdeg, unsigned dendeg,
                RatFunVarDegrees && degs)
    {
      if (nvars_ == 1) {
        sample_univariate_(f, mod);
        return;
      }

      xdegs_ = std::move(degs);
      get_monomial_prefactors_();

      sample_(f, mod,
              numdeg - (num_pref_ ? num_pref_.degree() : 0),
              dendeg - (den_pref_ ? den_pref_.degree() : 0));
    }

    void sample(RatFun & f, Mod mod,
                unsigned numdeg, unsigned dendeg,
                const RatFunVarDegrees & degs)
    {
      RatFunVarDegrees newdegs;
      newdegs.copy(degs, nvars_);
      sample(f, mod, numdeg, dendeg, std::move(newdegs));
    }

    const SparseRationalFunction & getFunction() const
    {
      return fun_;
    }

    SparseRationalFunction & getFunction()
    {
      return fun_;
    }

  private:
    Ret vars_degrees_(RatFun & f, Mod mod);

    void get_monomial_prefactors_();

    Ret degree_(RatFun & f, Mod mod, unsigned & numdeg, unsigned & dendeg);

    Ret degree_univariate_(RatFun & f, Mod mod,
                           unsigned & numdeg, unsigned & dendeg);

    Ret reconstruct_(RatFun & f, Mod mod, unsigned numdeg, unsigned dendeg);

    Ret reconstruct_univariate_(RatFun & f, Mod mod,
                                unsigned numdeg, unsigned dendeg);

    Ret reconstruct_univariate_(RatFun & f, Mod mod);

    void sample_(RatFun & f, Mod mod, unsigned numdeg, unsigned dendeg);
    void sample_univariate_(RatFun & f, Mod mod);

  private:
    SparseRationalFunction fun_;
    NPolyReconstruction npolyrec_;
    SparsePoly tmp_;
    std::vector<UInt> shift_;
    RatFunVarDegrees xdegs_;
    UInt n0_, t0_;
    std::size_t nvars_, maxdeg_;
    Monomial num_pref_, den_pref_;

  public:
    std::size_t n_checks, n_uchecks, n_undchecks, n_singular;
  };


} // namespace fflow

#endif // FFLOW_MULTIVARIATE_RECONSTRUCTION_HH
