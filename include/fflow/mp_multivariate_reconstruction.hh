#ifndef FFLOW_MP_MULTIVARIATE_RECONSTRUCTION_HH
#define FFLOW_MP_MULTIVARIATE_RECONSTRUCTION_HH

#include <fflow/multivariate_reconstruction.hh>
#include <fflow/mp_functions.hh>

namespace fflow {

  class MPRatFunReconstruction {
  public:

    enum {DEFAULT_N_UCHECKS = RatFunReconstruction::DEFAULT_N_UCHECKS};
    enum {DEFAULT_N_CHECKS = RatFunReconstruction::DEFAULT_N_CHECKS};
    enum {DEFAULT_N_SINGULAR = RatFunReconstruction::DEFAULT_N_SINGULAR};
    enum {DEFAULT_MAX_DEG = RatFunReconstruction::DEFAULT_MAX_DEG};
    enum {DEFAULT_MAX_PRIMES = 5};

    static bool VERBOSE;

    explicit MPRatFunReconstruction(unsigned nvars,
                                    unsigned maxdegree = DEFAULT_MAX_DEG)
      : shift_(), t0_(OFFSET_0_HASH),
        nvars_(nvars), maxdeg_(maxdegree),
        n_checks(DEFAULT_N_CHECKS), n_uchecks(DEFAULT_N_UCHECKS),
        n_singular(DEFAULT_N_SINGULAR),
        start_mod(0), max_primes(DEFAULT_MAX_PRIMES) {}

    unsigned getMaxDegree() const
    {
      return maxdeg_;
    }

    void setMaxDegree(unsigned deg)
    {
      maxdeg_ = deg;
    }

    UInt getT0() const
    {
      return t0_;
    }

    void setT0(UInt t0)
    {
      t0_ = t0;
    }

    void setShift(const UInt c[])
    {
      shift_.resize(nvars_);
      for (unsigned i=0; i<nvars_; ++i)
        shift_[i] = c[i];
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

    unsigned nvars() const
    {
      return nvars_;
    }

    // compute the degree
    Ret degree(RatFun & f, Mod mod, unsigned & numdeg, unsigned & dendeg);

    // reconstruct
    Ret reconstruct(RatFun & f, MPReconstructedRatFun & recf)
    {
      return reconstruct_(f, 0, 0, nullptr, recf);
    }

    Ret reconstructWithDegs(RatFun & f,
                            unsigned numdeg, unsigned dendeg,
                            const RatFunVarDegrees & degs,
                            MPReconstructedRatFun & recf)
    {
      return reconstruct_(f, numdeg, dendeg, &degs, recf);
    }

    // compute min and max degrees wrt one variable
    Ret var_degree(RatFun & f, unsigned var, Mod mod,
                   std::size_t & numdeg_min, std::size_t & numdeg_max,
                   std::size_t & dendeg_min, std::size_t & dendeg_max);

    void sample(RatFun & f,
                unsigned numdeg, unsigned dendeg,
                const RatFunVarDegrees & degs);

    void setup_rat_rec(RatFunReconstruction & rec)
    {
      rec.setT0(t0_);
      if (!shift_.empty())
        rec.setShift(shift_.data());
      rec.setMaxDegree(maxdeg_);
      rec.n_checks = n_checks;
      rec.n_uchecks = n_uchecks;
      rec.n_singular = n_singular;
    }

  private:

    Ret reconstruct_(RatFun & f,
                     unsigned numdeg, unsigned dendeg,
                     const RatFunVarDegrees * degs,
                     MPReconstructedRatFun & recf);

  private:
    std::vector<UInt> shift_;
    UInt t0_;
    std::size_t nvars_, maxdeg_;

  public:
    std::size_t n_checks, n_uchecks, n_singular;
    unsigned start_mod, max_primes, polymethod;
  };

} // namespace fflow


#endif // FFLOW_MP_MULTIVARIATE_RECONSTRUCTION_HH
