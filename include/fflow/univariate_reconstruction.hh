#ifndef FFLOW_UNIVARIATE_RECONSTRUCTION_HH
#define FFLOW_UNIVARIATE_RECONSTRUCTION_HH

#include <fflow/polynomial.hh>
#include <fflow/matrix.hh>
#include <fflow/rational_function.hh>

namespace fflow {

  // Generic interface for univariate rational finctions and
  // polynomials
  class URatFun {
  public:
    // must return FAILED if the value of x makes the calculation
    // impossible/meaningless
    virtual UInt evaluate(UInt x, Mod mod) = 0;
    virtual ~URatFun() {}
  };


  const UInt SAMPLING_STRIDE = 12345678901ULL;


  // Reconstruct univariate polynomials using Newton interpolation
  // method
  class UPolyReconstruction {
  public:
    enum {DEFAULT_N_CHECKS = 3, DEFAULT_N_SINGULAR = 10};

    explicit UPolyReconstruction(unsigned maxdegree = 0)
      : poly_(0), maxdeg_(maxdegree),
        n_checks(DEFAULT_N_CHECKS), n_singular(DEFAULT_N_SINGULAR),
        check_over_maxdeg(true) {}

    UPolyReconstruction(unsigned mindegree, unsigned maxdegree)
      : poly_(mindegree), maxdeg_(maxdegree),
        n_checks(DEFAULT_N_CHECKS), n_singular(DEFAULT_N_SINGULAR),
        check_over_maxdeg(true) {}

    unsigned getMinDegree() const
    {
      return poly_.getDegree();
    }

    unsigned getMaxDegree() const
    {
      return maxdeg_;
    }

    void setMinDegree(unsigned deg)
    {
      if (poly_.getDegree() < deg)
        poly_.setDegree(deg);
    }

    void setMaxDegree(unsigned deg)
    {
      maxdeg_ = deg;
      if (poly_.getDegree() > deg)
        poly_.setDegree(deg);
    }

    void writeResult(Mod mod, UPoly & p)
    {
      UInt * coeffs = p.coeffs();
      poly_.toStdRep(mod,coeffs);
    }

    UInt getX0() const
    {
      return poly_.getX0();
    }

    UInt getXi(std::size_t i) const
    {
      return poly_.getXi(i);
    }

    void setX0(UInt x0)
    {
      poly_.setX0(x0);
    }

    Ret reconstruct(URatFun & f, Mod mod);
    void sample(URatFun & f, Mod mod);

    unsigned getTrueDegree();

    const NewtonUPoly & getNewtonPoly() const
    {
      return poly_;
    }

    NewtonUPoly & getNewtonPoly()
    {
      return poly_;
    }

    void remove_high_zeroes()
    {
      poly_.remove_high_zeroes();
    }

  private:
    NewtonUPoly poly_;
    std::size_t maxdeg_;

  public:
    std::size_t n_checks, n_singular;
    bool check_over_maxdeg;
  };



  // Reconstruct univariate rational functions using Thiele's
  // interpolation method
  class URatFunReconstruction {
  public:
    enum {DEFAULT_N_CHECKS = 3};
    enum {DEFAULT_N_SINGULAR = 10};

    explicit URatFunReconstruction(unsigned maxterms = 0)
      : fun_(), maxterms_(maxterms),
        n_checks(DEFAULT_N_CHECKS), n_singular(DEFAULT_N_SINGULAR) {}

    URatFunReconstruction(unsigned minterms, unsigned maxterms)
      : fun_(minterms), maxterms_(maxterms),
        n_checks(DEFAULT_N_CHECKS), n_singular(DEFAULT_N_SINGULAR) {}

    unsigned getMinTerms() const
    {
      return fun_.size();
    }

    unsigned getMaxTerms() const
    {
      return maxterms_;
    }

    void setMinTerms(unsigned t)
    {
      if (fun_.size() < t)
        fun_.resize(t);
    }

    void setMaxTerms(unsigned t)
    {
      maxterms_ = t;
      if (fun_.size() > t)
        fun_.resize(t);
    }

    void reset()
    {
      fun_.clear();
    }

    void reset(unsigned maxterms)
    {
      maxterms_ = maxterms;
      fun_.clear();
    }

    void reset(unsigned minterms, unsigned maxterms)
    {
      maxterms_ = maxterms;
      fun_.clear(minterms);
    }

    void writeResult(Mod mod, URationalFunction & p)
    {
      fun_.toStdRep(mod,p);
    }

    UInt getX0() const
    {
      return fun_.getX0();
    }

    UInt getXi(std::size_t i) const
    {
      return fun_.getXi(i);
    }

    void setX0(UInt x0)
    {
      fun_.setX0(x0);
    }

    Ret reconstruct(URatFun & f, Mod mod);
    void sample(URatFun & f, Mod mod);

  private:
    ThieleURationalFunction fun_;
    std::size_t maxterms_;

  public:
    std::size_t n_checks, n_singular;
  };



  // Reconstruct univariate the numerator and denominator of a
  // rational functions with a straightforward linear solve method, in
  // the case where den[0] is normalizable to 1 and the degrees of
  // numerator and denominators are known
  class URatFunNumDemReconstruction {
  public:
    enum {DEFAULT_N_CHECKS = 1, DEFAULT_N_SINGULAR = 10};

    URatFunNumDemReconstruction()
      : mat_(), fun_(0,0), t0_(1),
        n_checks(DEFAULT_N_CHECKS), n_singular(DEFAULT_N_SINGULAR) {}

    URatFunNumDemReconstruction(std::size_t numdeg, std::size_t dendeg)
      : mat_(), fun_(numdeg, dendeg), t0_(1),
        n_checks(DEFAULT_N_CHECKS), n_singular(DEFAULT_N_SINGULAR) {}

    void setNumDegree(std::size_t deg)
    {
      fun_.setNumDegree(deg);
    }

    void setDenDegree(std::size_t deg)
    {
      fun_.setDenDegree(deg);
    }

    void setDegree(std::size_t numdeg, std::size_t dendeg)
    {
      fun_.setNumDegree(numdeg);
      fun_.setDenDegree(dendeg);
    }

    Ret reconstruct(URatFun & f, Mod mod);
    void sample(URatFun & f, Mod mod);

    URationalFunction & getFunction()
    {
      return fun_;
    }

    const URationalFunction & getFunction() const
    {
      return fun_;
    }

    void setX0(UInt x0)
    {
      t0_ = x0;
    }

    UInt getX0() const
    {
      return t0_;
    }

  private:
    DynamicMatrix mat_;
    URationalFunction fun_;
    UInt t0_;

  public:
    std::size_t n_checks, n_singular;
  };


  class URatFunNumDenRecoHighDegs {
  public:
    enum {DEFAULT_N_CHECKS = 1, DEFAULT_N_SINGULAR = 10};

    URatFunNumDenRecoHighDegs()
      : mat_(), fun_(0,0), t0_(1),
        n_checks(DEFAULT_N_CHECKS), n_singular(DEFAULT_N_SINGULAR),
        rnum_min(0), rden_min(1) {}

    URatFunNumDenRecoHighDegs(std::size_t numdeg, std::size_t dendeg)
      : mat_(), fun_(numdeg, dendeg), t0_(1),
        n_checks(DEFAULT_N_CHECKS), n_singular(DEFAULT_N_SINGULAR),
        rnum_min(0), rden_min(1) {}

    void setMinUnknownDegs(unsigned num, unsigned den)
    {
      rnum_min = num;
      rden_min = den;
    }

    void setNumDegree(std::size_t deg)
    {
      fun_.setNumDegree(deg);
    }

    void setDenDegree(std::size_t deg)
    {
      fun_.setDenDegree(deg);
    }

    void setDegree(std::size_t numdeg, std::size_t dendeg)
    {
      fun_.setNumDegree(numdeg);
      fun_.setDenDegree(dendeg);
    }

    Ret reconstruct(URatFun & f, Mod mod);
    void sample(URatFun & f, Mod mod);

    URationalFunction & getFunction()
    {
      return fun_;
    }

    const URationalFunction & getFunction() const
    {
      return fun_;
    }

    void setX0(UInt x0)
    {
      t0_ = x0;
    }

    UInt getX0() const
    {
      return t0_;
    }

  private:
    DynamicMatrix mat_;
    URationalFunction fun_;
    UInt t0_;

  public:
    std::size_t n_checks, n_singular;
    unsigned rnum_min, rden_min;
  };


} // namespace fflow


#endif // FFLOW_UNIVARIATE_RECONSTRUCTION_HH
