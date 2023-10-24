#include <iostream>
#include <unordered_map>
#include <memory>

#include <fflow/multivariate_reconstruction.hh>
#include <fflow/univariate_reconstruction.hh>
#include <fflow/function_cache.hh>

namespace fflow {

  /////////////////////////
  // Multivariate Newton //
  /////////////////////////

  namespace detail {

    // univariate function obtained by "freezing" the first n-1
    // variables of a multivariate function
    struct URatFunFreezed : public URatFun {

      URatFunFreezed(RatFun & fun, unsigned nvars)
        : f_(fun), xi_(nvars), ai_(nvars-1), i_(nvars-1) {}

      virtual UInt evaluate(UInt x, Mod mod)
      {
        xi_.back() = x;
        UInt res = f_.evaluate(xi_.data(),mod);
        if (res == FAILED)
          return res;

        unsigned nv = nvars();

        for (unsigned i=0; i<nv-1; ++i) {

          if (i_[i] == 0)
            continue;

          NewtonPoly & ri = (**ai_[i]);

          for (unsigned j=0; j<=i_[i]-1; ++j) {
            NewtonPoly & ai = *(ri.pcoeffs()[j]);
            res = sub_mod(res,ai.eval(xi_[i+1],xi_.data()+i+2,mod),mod);
            res = div_mod(res, sub_mod(xi_[i],ri.getXi(j),mod), mod);
          }

        }

        return res;
      }

      unsigned nvars() const
      {
        return xi_.size();
      }

      RatFun & f_;
      std::vector<UInt> xi_;
      std::vector<NewtonPoly::Ptr *> ai_;
      std::vector<unsigned> i_;
    };


    struct URatFunFreezedSample : public URatFun {

      URatFunFreezedSample(RatFun & fun, unsigned nvars)
        : f_(fun), xi_(nvars) {}

      virtual UInt evaluate(UInt x, Mod mod)
      {
        xi_.back() = x;
        f_.evaluate(xi_.data(),mod);
        return 0;
      }

      unsigned nvars() const
      {
        return xi_.size();
      }

      RatFun & f_;
      std::vector<UInt> xi_;
    };

  } // namespace detail


  Ret NPolyReconstruction::reconstruct(RatFun & f, Mod mod)
  {
    detail::URatFunFreezed fr(f,nvars_);
    (*poly_).setX0(x0_[0]);
    NewtonPoly::Ptr res = reconstruct_(fr,nvars_,0,0,mod);

    if (res.get() != nullptr) {
      poly_ = std::move(res);
      return SUCCESS;
    }

    return FAILED;
  }


  NewtonPoly::Ptr NPolyReconstruction::reconstruct_(detail::URatFunFreezed & f,
                                                    unsigned nvars,
                                                    unsigned fixed_vars,
                                                    unsigned current_deg,
                                                    Mod mod)
  {
    if (!check_over_maxdeg && current_deg > maxdeg_)
      return nullptr;

    if (fixed_vars == nvars-1) {

      // fit univariate polynomial

      std::size_t this_deg = current_deg > maxdeg_ ? 0 : maxdeg_ - current_deg;
      if (xdegs)
        this_deg = std::min(this_deg, maxdeg_var_(fixed_vars));
      UPolyReconstruction poly_rec(this_deg);
      poly_rec.setX0(x0_[fixed_vars]);
      poly_rec.n_checks = n_uchecks;
      poly_rec.n_singular = n_singular;
      if (!check_over_maxdeg)
        poly_rec.check_over_maxdeg = false;

      if (!(poly_rec.reconstruct(f,mod) == SUCCESS))
        return nullptr;

      poly_rec.remove_high_zeroes();
      NewtonPoly * res = new NewtonUPoly(std::move(poly_rec.getNewtonPoly()));
      return NewtonPoly::Ptr(res);

    } else {

      // fit multivariate
      const UInt stride = SAMPLING_STRIDE+(fixed_vars+1)*64;
      unsigned check_count = 0;
      NewtonPoly::Ptr resptr(new NewtonMPoly(nvars-fixed_vars));
      NewtonMPoly & res(static_cast<NewtonMPoly&>(*resptr));

      UInt x = x0_[fixed_vars];
      res.setX0(x0_[fixed_vars]);

      unsigned this_deg = std::min(maxdeg_ - current_deg,
                                   maxdeg_var_(fixed_vars));

      for (unsigned i=0; i<this_deg+n_checks+1; ++i) {

        f.i_[fixed_vars] = i;
        f.ai_[fixed_vars] = & resptr;

        NewtonPoly::Ptr a;
        for (unsigned ev=0; ev<=n_singular; ++ev) {
          f.xi_[fixed_vars] = x;
          a = reconstruct_(f,nvars,fixed_vars+1,current_deg+i,mod);
          if (a != nullptr)
            break;
          else if (ev == n_singular)
            return nullptr;
          x = add_mod(x, stride, mod);
        }

        bool a_is_zero = (*a).is_zero_quick();

        if (a_is_zero) {
          ++check_count;
          if (check_count >= n_checks) {
            res.remove_high_zeroes_quick();
            return resptr;
          }
        } else {
          check_count = 0;
        }

        if (res.size() <= i) {
          res.push_back(x, std::move(a));
        } else {
          res.setXi(i, x);
          res.pcoeffs()[i] = std::move(a);
        }

        if (!check_over_maxdeg && i==this_deg) {
          res.remove_high_zeroes_quick();
          return resptr;
        }

        if (!a_is_zero && res.getDegree()>this_deg)
          return nullptr;

        x = add_mod(x, stride, mod);
      }

    }

    return nullptr;
  }


  void NPolyReconstruction::sample(RatFun & f, Mod mod)
  {
    detail::URatFunFreezedSample fr(f,nvars_);
    (*poly_).setX0(x0_[0]);
    sample_(fr,nvars_,0,0,mod);
  }


  void NPolyReconstruction::sample_(detail::URatFunFreezedSample & f,
                                    unsigned nvars,
                                    unsigned fixed_vars,
                                    unsigned current_deg,
                                    Mod mod)
  {
    if (fixed_vars == nvars-1) {

      // sample univariate polynomial

      std::size_t this_deg = current_deg > maxdeg_ ? 0 : maxdeg_ - current_deg;
      if (xdegs)
        this_deg = std::min(this_deg, maxdeg_var_(fixed_vars));
      UPolyReconstruction poly_rec(this_deg);
      poly_rec.setX0(x0_[fixed_vars]);
      poly_rec.n_checks = n_uchecks;
      poly_rec.n_singular = n_singular;
      if (!check_over_maxdeg)
        poly_rec.check_over_maxdeg = false;

      poly_rec.sample(f,mod);

    } else {

      // sample multivariate

      const UInt stride = SAMPLING_STRIDE+(fixed_vars+1)*64;

      UInt x = x0_[fixed_vars];

      unsigned this_deg = std::min(maxdeg_ - current_deg,
                                   maxdeg_var_(fixed_vars));

      std::size_t n_samples = this_deg + 1 + n_singular
        + (check_over_maxdeg ? n_checks : 0);

      for (unsigned i=0; i<n_samples; ++i) {
        f.xi_[fixed_vars] = x;
        unsigned min_current_deg = std::max(0,
                                            int(current_deg+i)-int(n_singular));
        sample_(f,nvars,fixed_vars+1,min_current_deg,mod);
        x = add_mod(x, stride, mod);
      }

    }
  }


  //////////////////////////////
  // Multivariate Vandermonde //
  //////////////////////////////


  // Generate a list of monomials and count its length
  unsigned VPolyReconstruction::gen_mons_(unsigned curr_deg,
                                          unsigned var, UInt * toadd)
  {
    unsigned maxdeg = std::min(maxdeg_-curr_deg, maxdeg_var_(var));
    unsigned mindeg = (curr_deg >= mindeg_ || var != nvars_-1) ? 0
      : mindeg_ - curr_deg;
    unsigned nmons = 0;

    if (var == nvars_-1) {

      if (toadd)
        for (unsigned j=mindeg; j<=maxdeg; ++j) {
          toadd[var] = j;
          auto mon = Monomial(nvars_);
          mon.degree() = curr_deg + j;
          for (unsigned k=0; k<nvars_; ++k)
            mon.exponent(k) = toadd[k];
          mon.coeff() = 1;
          mons_.emplace_back(std::move(mon));
        }

      return maxdeg >= mindeg ? maxdeg - mindeg + 1 : 0;

    } else {

      for (unsigned j=mindeg; j<=maxdeg; ++j) {
        if (toadd)
          toadd[var] = j;
        nmons += gen_mons_(curr_deg+j, var+1, toadd);
      }

    }

    return nmons;
  }


  void VPolyReconstruction::sample(RatFun & f, Mod mod)
  {
    unsigned npoints = gen_mons_();
    if (npoints == 0)
      return;

    f.evaluate(x0_.get(), mod);
    if (npoints == 1)
      return;

    std::unique_ptr<UInt[]> x(new UInt[nvars_]);
    std::copy(x0_.get(), x0_.get()+nvars_, x.get());

    for (unsigned j=1; j<npoints; ++j) {
      for (unsigned k=0; k<nvars_; ++k)
        x[k] = mul_mod(x0_[k], x[k], mod);
      f.evaluate(x.get(), mod);
    }
  }

  Ret VPolyReconstruction::reconstruct(RatFun & f, Mod mod)
  {
    mons_.clear();

    std::unique_ptr<UInt[]> x(new UInt[nvars_]);
    const unsigned npoints = gen_mons_(0,0,x.get());

    if (npoints == 0)
      return SUCCESS;

    std::unique_ptr<UInt[]> v(new UInt[npoints]);
    std::unique_ptr<UInt[]> y(new UInt[npoints]);

    // evaluate monomials at x0
    const auto op = MonomialOp{nvars_};
    for (unsigned j=0; j<npoints; ++j)
      v[j] = op.eval(mons_[j], x0_.get(), mod);

    y[0] = f.evaluate(x0_.get(), mod);

    if (npoints == 1) {
      if (v[0] == 0)
        return FAILED;
      mons_[0].coeff() = div_mod(y[0],v[0],mod);
    }

    // get all other evaluation points
    std::copy(x0_.get(), x0_.get()+nvars_, x.get());
    for (unsigned j=1; j<npoints; ++j) {
      for (unsigned k=0; k<nvars_; ++k)
        x[k] = mul_mod(x0_[k], x[k], mod);
      y[j] = f.evaluate(x.get(), mod);
      if (y[j] == FAILED)
        return FAILED;
    }

    // p(z) = (z-v[0])*(z-v[1])*...
    std::unique_ptr<UInt[]> p(new UInt[npoints+1]);
    std::fill(p.get(), p.get()+npoints,0);
    p[npoints]=1;
    p[npoints-1] = neg_mod(v[0], mod);
    for (unsigned j=1; j<npoints; ++j)
      for (unsigned k=npoints-j-1; k<npoints; ++k)
        p[k] = sub_mod(p[k], mul_mod(p[k+1], v[j], mod), mod);

    // Build the matrix M such that M.y = c
    // It can be computed using:
    //
    //   p_j(z) = \sum_k M_jk * z^k
    //
    // where
    //
    //   p_j(z) = 1/norm*z*p(z)/(z-v[j])
    //
    // with norm = v[j] * \prod_k (v[j]-v[k]) for k != j
    std::unique_ptr<UInt[]> pj(new UInt[npoints]);
    for (unsigned j=0; j<npoints; ++j) {

      // q(z) = p_j(z)*z
      UInt * q = pj.get();

      const UInt vj = v[j];
      const UInt vjp = precomp_mul_shoup(vj, mod);

      // quotient
      q[npoints-1] = 1;
      for (int k=npoints-2; k>=0; --k)
        q[k] = add_mod(p[k+1], mul_mod_shoup(q[k+1],vj,vjp,mod), mod);

      UInt norm=1;
      for (unsigned k=0; k<npoints; ++k)
        if (k!=j)
          norm = mul_mod(norm, sub_mod(vj, v[k], mod), mod);
      norm = mul_mod(norm, vj, mod);
      if (!norm)
        return FAILED;

      const UInt inorm = mul_inv(norm, mod);
      const UInt inormp = precomp_mul_shoup(inorm, mod);
      for (unsigned k=0; k<npoints; ++k)
        q[k] = mul_mod_shoup(q[k], inorm, inormp, mod);

      UInt tmp = 0;
      for (unsigned k=0; k<npoints; ++k)
        tmp = add_mod(tmp, mul_mod(q[k],y[k],mod), mod);
      mons_[j].coeff() = tmp;

    }

    return SUCCESS;
  }

  void VPolyReconstruction::writeResultVar(std::size_t first_var, Mod,
                                           SparsePoly & p)
  {
    p = SparsePoly(nvars_ + first_var);
    if (mons_.size())
      p.fromMonomials(mons_.data(), mons_.size(), first_var);
  }

  void VPolyReconstruction::writeResult(Mod mod, SparsePoly & p)
  {
    writeResultVar(0, mod, p);
  }



  ////////////////////////
  // Rational functions //
  ////////////////////////

  namespace  {

    // univariate rational function from multivariate one, with all
    // variables fixed but one
    struct UFromMRatFun : public URatFun {

      UFromMRatFun(RatFun & f, std::size_t nvars)
        : xi_(new UInt[nvars]), f_(f), var_(0) {}

      virtual UInt evaluate(UInt x, Mod mod);

      std::unique_ptr<UInt[]> xi_; // assigned externally
      RatFun & f_;
      std::size_t var_;
    };

    UInt UFromMRatFun::evaluate(UInt x, Mod mod)
    {
      xi_[var_] = x;
      return f_.evaluate(xi_.get(), mod);
    }


    // "homogenized" univariate rational function from multivariate
    // one
    struct UHomogenizeMRatFun : public URatFun {

      UHomogenizeMRatFun(RatFun & f, std::size_t nvars)
        : xi_(new UInt[nvars]), mx_(new UInt[nvars]), shift_(nullptr),
          num_pref_(nullptr), den_pref_(nullptr),
          f_(f), n_(nvars) {}

      virtual UInt evaluate(UInt x, Mod mod);

      std::unique_ptr<UInt[]> xi_; // assigned externally
      std::unique_ptr<UInt[]> mx_; // used internally
      const UInt * shift_; // assigned externally
      const Monomial * num_pref_;
      const Monomial * den_pref_;
      RatFun & f_;
      std::size_t n_;
    };

    UInt UHomogenizeMRatFun::evaluate(UInt t, Mod mod)
    {
      UInt res;
      mx_[0] = t;
      for (unsigned i=1; i<n_; ++i)
        mx_[i] = mul_mod(t, xi_[i], mod);

      if (shift_)
        for (unsigned i=0; i<n_; ++i)
          mx_[i] = add_mod(mx_[i], shift_[i], mod);

      res = f_.evaluate(mx_.get(),mod);

      if (res != FAILED && res) {
        MonomialOp mop{n_};
        if (num_pref_)
          res = div_mod(res, mop.eval(*num_pref_,mx_.get(),mod), mod);
        if (den_pref_)
          res = mul_mod(res, mop.eval(*den_pref_,mx_.get(),mod), mod);
      }

      return res;
    }


    struct UHomogenizeMRatFunSample : public URatFun {

      UHomogenizeMRatFunSample(RatFun & f, std::size_t nvars)
        : xi_(new UInt[nvars]), mx_(new UInt[nvars]), shift_(nullptr),
          f_(f), n_(nvars) {}

      virtual UInt evaluate(UInt x, Mod mod);

      std::unique_ptr<UInt[]> xi_; // assigned externally
      std::unique_ptr<UInt[]> mx_; // used internally
      const UInt * shift_; // assigned externally
      RatFun & f_;
      std::size_t n_;
    };

    UInt UHomogenizeMRatFunSample::evaluate(UInt t, Mod mod)
    {
      mx_[0] = t;
      for (unsigned i=1; i<n_; ++i)
        mx_[i] = mul_mod(t, xi_[i], mod);

      if (shift_)
        for (unsigned i=0; i<n_; ++i)
          mx_[i] = add_mod(mx_[i], shift_[i], mod);

      f_.evaluate(mx_.get(), mod);

      return 0;
    }


    typedef FunctionCache<0> IntArrayXiHashTable;


    // evaluates the polynomials defined by the homogenized rational
    // function
    struct HomogenizeMRatFun : public RatFun {

      enum {NUMERATOR, DENOMINATOR};

      enum {BAD_SHIFT = fflow::BAD_SHIFT};

      HomogenizeMRatFun(RatFun & f, std::size_t nvars)
        : hash_(binomial(nvars+5,5),
                IntArrayXiHash<0>{nvars-1},
                IntArrayXiEqual<0>{nvars-1}),
          rat_(), ndrat_(), f_(f,nvars), ww_(),
          hornernum_(nullptr), hornerden_(nullptr), mx_(nullptr), ai_(nullptr),
          n0_(0), numdeg_(0), dendeg_(0),
          i_(0), type_(0),
          n_checks(RatFunReconstruction::DEFAULT_N_UCHECKS),
          n_ndchecks(RatFunReconstruction::DEFAULT_N_UNDCHECKS),
          n_learning(RatFunReconstruction::DEFAULT_N_UCHECKS+1),
          n_singular(RatFunReconstruction::DEFAULT_N_SINGULAR),
          is_learning(true) {}

      Ret reconstructThiele(Mod mod);
      Ret reconstructNumDen(Mod mod);

      virtual UInt evaluate(const UInt x[], Mod mod);

      UInt evaluate_learning(const UInt x[], Mod mod);

      URationalFunction & urat()
      {
        return ndrat_.getFunction();
      }

      const URationalFunction & urat() const
      {
        return ndrat_.getFunction();
      }

      std::size_t numterms() const
      {
        return urat().getNumDegree() + 1;
      }

      std::size_t denterms() const
      {
        return urat().getDenDegree() + 1;
      }

      void set_degree(unsigned numdegree, unsigned dendegree)
      {
        ndrat_.setNumDegree(numdegree);
        ndrat_.setDenDegree(dendegree);
      }

      void init_hornerpoly()
      {
        hornernum_.reset(new HornerPtr[numterms()]);
        hornerden_.reset(new HornerPtr[denterms()]);
      }

      IntArrayXiHashTable hash_;
      URatFunReconstruction rat_;
      URatFunNumDenRecoHighDegs ndrat_;
      UHomogenizeMRatFun f_;
      HornerWorkspacePtr ww_;
      std::unique_ptr<HornerPtr[]> hornernum_, hornerden_;
      std::unique_ptr<UInt[]> mx_, mxp_, ai_;
      UInt n0_;
      std::size_t numdeg_, dendeg_;
      std::size_t i_, type_; // assigned externally
      std::size_t n_checks, n_ndchecks, n_learning, n_singular;
      bool is_learning;
    };

    UInt HomogenizeMRatFun::evaluate_learning(const UInt x[], Mod mod)
    {
      if (mx_.get() == nullptr)
        mx_.reset(new UInt[f_.n_ - 1]);

      for (unsigned i=0; i<f_.n_-1; ++i)
        mx_[i] = x[i];

      if (reconstructThiele(mod) != SUCCESS)
        return FAILED;

      if (urat().den()[0] == 0)
        return BAD_SHIFT;

      numdeg_ = std::max(numdeg_, urat().getNumDegree());
      dendeg_ = std::max(dendeg_, urat().getDenDegree());

      return (type_ == DENOMINATOR) ? urat().den()[i_] : urat().num()[i_];
    }

    UInt HomogenizeMRatFun::evaluate(const UInt x[], Mod mod)
    {
      static std::size_t count = 0;

      if (mx_.get() == nullptr)
        mx_.reset(new UInt[f_.n_ - 1]);

      for (unsigned i=0; i<f_.n_-1; ++i)
        mx_[i] = x[i];

      // look for result
      auto lookup = hash_.find(mx_);

      if (lookup != hash_.end()) {

        if (lookup->second == nullptr)
          return FAILED;

        std::size_t index = i_;
        if (type_ == DENOMINATOR)
          index += numterms();

        return lookup->second[index];

      } else {

        if (RatFunReconstruction::VERBOSE) {
          ++count;
          if (!((count) % 100))
            std::cout << "h rec. count = " << count << std::endl;
        }

        // fit univariate rational function
        if (reconstructNumDen(mod) != SUCCESS) {
          hash_.emplace(std::make_pair(std::move(mx_), nullptr));
          return FAILED;
        }

        // store the result for later lookup

        ai_.reset(new UInt[numterms()+denterms()]);
        UInt * aiptr = ai_.get();

        std::copy(urat().num().coeffs(), urat().num().coeffs() + numterms(),
                  ai_.get());

        std::copy(urat().den().coeffs(), urat().den().coeffs() + denterms(),
                  ai_.get() + numterms());

        hash_.emplace(std::make_pair(std::move(mx_), std::move(ai_)));


        // return

        std::size_t index = i_;
        if (type_ == DENOMINATOR)
          index += numterms();

        return aiptr[index];

      }

      // it should never reach this point
      return FAILED;
    }

    Ret HomogenizeMRatFun::reconstructThiele(Mod mod)
    {
      rat_.reset();

      for (unsigned i=1; i<f_.n_; ++i)
        f_.xi_[i] = mx_[i-1];

      rat_.n_checks = n_checks;
      rat_.n_singular = n_singular;

      if (rat_.reconstruct(f_,mod) != SUCCESS)
        return FAILED;

      rat_.writeResult(mod, urat());

      return SUCCESS;
    }


    Ret HomogenizeMRatFun::reconstructNumDen(Mod mod)
    {
      if (mxp_.get() == nullptr)
        mxp_.reset(new UInt[f_.n_ - 1]);

      for (unsigned i=1; i<f_.n_; ++i)
        f_.xi_[i] = mx_[i-1];

      ndrat_.n_checks = n_ndchecks;
      ndrat_.n_singular = n_singular;
      unsigned nrmin = std::min((type_ == HomogenizeMRatFun::DENOMINATOR) ?
                                i_+1 : i_,
                                numdeg_+1);
      unsigned drmin = ndrat_.rden_min = std::min(i_,dendeg_+1);
      if (drmin == 0)
        drmin = 1;
      ndrat_.rnum_min = nrmin;
      ndrat_.rden_min = drmin;

      {
        precomp_array_mul_shoup(mx_.get(), f_.n_-1, mod, mxp_.get());

        UInt * uc = urat().num().coeffs();
        if (nrmin)
          uc[0] = n0_;
        for (unsigned i=1; i<nrmin; ++i)
          uc[i] = horner_eval(hornernum_[i].get(), f_.n_-1, ww_.get(),
                              mx_.get(), mxp_.get(), mod);

        uc = urat().den().coeffs();
        uc[0] = 1;
        for (unsigned i=1; i<drmin; ++i)
          uc[i] = horner_eval(hornerden_[i].get(), f_.n_-1, ww_.get(),
                              mx_.get(), mxp_.get(), mod);
      }

      if (ndrat_.reconstruct(f_,mod)  != SUCCESS) {
        // // try recovering with Thiele
        // if (reconstructThiele(mod) != SUCCESS)
        return FAILED;
      }

      if (urat().getNumDegree() != numdeg_ || (numdeg_ && !urat().num().lt())) {
        set_degree(numdeg_, dendeg_);
        return FAILED;
      }

      return SUCCESS;
    }



    struct HomogenizeMRatFunSample : public RatFun {

      enum {NUMERATOR, DENOMINATOR};

      enum {BAD_SHIFT = fflow::BAD_SHIFT};

      HomogenizeMRatFunSample(RatFun & f, std::size_t nvars)
        : ndrat_(), f_(f,nvars),
          mx_(nullptr), ai_(nullptr),
          numdeg_(0), dendeg_(0),
          i_(0), type_(0),
          n_ndchecks(RatFunReconstruction::DEFAULT_N_UNDCHECKS),
          n_singular(RatFunReconstruction::DEFAULT_N_SINGULAR) {}

      Ret sampleNumDen(Mod mod);

      virtual UInt evaluate(const UInt x[], Mod mod);

      URationalFunction & urat()
      {
        return ndrat_.getFunction();
      }

      const URationalFunction & urat() const
      {
        return ndrat_.getFunction();
      }

      std::size_t numterms() const
      {
        return urat().getNumDegree() + 1;
      }

      std::size_t denterms() const
      {
        return urat().getDenDegree() + 1;
      }

      void set_degree(unsigned numdegree, unsigned dendegree)
      {
        ndrat_.setNumDegree(numdegree);
        ndrat_.setDenDegree(dendegree);
      }

      URatFunNumDenRecoHighDegs ndrat_;
      UHomogenizeMRatFunSample f_;
      std::unique_ptr<UInt[]> mx_, ai_;
      std::size_t numdeg_, dendeg_;
      std::size_t i_, type_; // assigned externally
      std::size_t n_ndchecks, n_singular;
    };

    Ret HomogenizeMRatFunSample::sampleNumDen(Mod mod)
    {
      for (unsigned i=1; i<f_.n_; ++i)
        f_.xi_[i] = mx_[i-1];

      ndrat_.n_checks = n_ndchecks;
      ndrat_.n_singular = n_singular;
      unsigned nrmin = std::min((type_ == HomogenizeMRatFun::DENOMINATOR) ?
                                i_+1 : i_,
                                numdeg_+1);
      unsigned drmin = ndrat_.rden_min = std::min(i_,dendeg_+1);
      if (drmin == 0)
        drmin = 1;
      ndrat_.rnum_min = nrmin;
      ndrat_.rden_min = drmin;

      ndrat_.sample(f_,mod);

      return 0;
    }

    UInt HomogenizeMRatFunSample::evaluate(const UInt x[], Mod mod)
    {
      if (mx_.get() == nullptr)
        mx_.reset(new UInt[f_.n_ - 1]);

      for (unsigned i=0; i<f_.n_-1; ++i)
        mx_[i] = x[i];

      sampleNumDen(mod);

      return 0;
    }

  } // namespace


  bool RatFunReconstruction::VERBOSE = false;

  Ret RatFunReconstruction::vars_degrees_(RatFun & f, Mod mod)
  {
    xdegs_.resize(nvars_);
    UFromMRatFun uf(f, nvars_);

    unsigned n_learning = 1; // std::max(n_checks,std::size_t(3));
    std::size_t n_fails = 0;
    std::unique_ptr<UInt[]> samples(new UInt[nvars_+n_learning+n_singular]);
    for (unsigned i=0; i<nvars_+n_learning+n_singular; ++i)
      samples[i] = sample_uint(OFFSET_1, i, mod);

    for (unsigned vv=0; vv<nvars_; ++vv) {

      uf.var_ = vv;

      for (unsigned ev=0; ev<n_learning; ++ev) {

        bool okay = false;

        do {

          for (unsigned i=0; i<nvars_; ++i)
            uf.xi_[i] = samples[i+ev+n_fails];

          URatFunReconstruction ndrat(2*maxdeg_ + 1);
          ndrat.n_checks = n_uchecks;
          ndrat.n_singular = n_singular;
          ndrat.setX0(t0_);

          Ret ret = ndrat.reconstruct(uf, mod);

          URationalFunction urat;
          ndrat.writeResult(mod, urat);

          if (ret == SUCCESS) {

            std::size_t num_maxdeg = std::max(xdegs_.num_maxdegs(nvars_)[vv],
                                              urat.getNumDegree());
            std::size_t den_maxdeg = std::max(xdegs_.den_maxdegs(nvars_)[vv],
                                              urat.getDenDegree());
            std::size_t num_mindeg = std::max(xdegs_.num_mindegs(nvars_)[vv],
                                              urat.num().getMinDegree());
            std::size_t den_mindeg = std::max(xdegs_.den_mindegs(nvars_)[vv],
                                              urat.den().getMinDegree());
            if (xdegs_.num_maxdegs(nvars_)[vv] == urat.getNumDegree())
              num_mindeg = std::min(xdegs_.num_mindegs(nvars_)[vv],
                                    urat.num().getMinDegree());
            if (xdegs_.num_maxdegs(nvars_)[vv] == urat.getDenDegree())
              den_mindeg = std::min(xdegs_.den_mindegs(nvars_)[vv],
                                    urat.den().getMinDegree());

            xdegs_.num_maxdegs(nvars_)[vv] = num_maxdeg;
            xdegs_.den_maxdegs(nvars_)[vv] = den_maxdeg;
            xdegs_.num_mindegs(nvars_)[vv] = num_mindeg;
            xdegs_.den_mindegs(nvars_)[vv] = den_mindeg;

            okay = true;

          } else {

            ++n_fails;
            if (n_fails == n_singular)
              return FAILED;

          }

        } while(!okay);
      }
    }

    return SUCCESS;
  }

  // must be called after vars_degrees_vars_degs_
  void RatFunReconstruction::get_monomial_prefactors_()
  {
    // get prefactors from min. degs
    unsigned num_pref_deg = std::accumulate(xdegs_.num_mindegs(nvars_),
                                            xdegs_.num_mindegs(nvars_)+nvars_,
                                            0);
    if (num_pref_deg) {
      num_pref_ = Monomial(nvars_);
      num_pref_.coeff() = 1;
      num_pref_.degree() = num_pref_deg;
      std::copy(xdegs_.num_mindegs(nvars_), xdegs_.num_mindegs(nvars_)+nvars_,
                num_pref_.exponents());
    }
    unsigned den_pref_deg = std::accumulate(xdegs_.den_mindegs(nvars_),
                                            xdegs_.den_mindegs(nvars_)+nvars_,
                                            0);
    if (den_pref_deg) {
      den_pref_ = Monomial(nvars_);
      den_pref_.coeff() = 1;
      den_pref_.degree() = den_pref_deg;
      std::copy(xdegs_.den_mindegs(nvars_), xdegs_.den_mindegs(nvars_)+nvars_,
                den_pref_.exponents());
    }

    // remove prefactors from max. degs
    for (unsigned vv=0; vv<nvars_; ++vv) {
      xdegs_.num_maxdegs(nvars_)[vv] -= xdegs_.num_mindegs(nvars_)[vv];
      xdegs_.den_maxdegs(nvars_)[vv] -= xdegs_.den_mindegs(nvars_)[vv];
    }
  }

  Ret RatFunReconstruction::var_degree(RatFun & f, unsigned vv, Mod mod,
                                       std::size_t & xnumdeg_min,
                                       std::size_t & xnumdeg_max,
                                       std::size_t & xdendeg_min,
                                       std::size_t & xdendeg_max)
  {
    UFromMRatFun uf(f, nvars_);

    xnumdeg_min = 0;
    xnumdeg_max = 0;
    xdendeg_min = 0;
    xdendeg_max = 0;

    unsigned n_learning = 1; // std::max(n_checks,std::size_t(3));
    std::size_t n_fails = 0;
    std::unique_ptr<UInt[]> samples(new UInt[nvars_+n_learning+n_singular]);
    for (unsigned i=0; i<nvars_+n_learning+n_singular; ++i)
      samples[i] = sample_uint(OFFSET_1, i, mod);

    {

      uf.var_ = vv;

      for (unsigned ev=0; ev<n_learning; ++ev) {

        bool okay = false;

        do {

          for (unsigned i=0; i<nvars_; ++i)
            uf.xi_[i] = samples[i+ev+n_fails];

          URatFunReconstruction ndrat(2*maxdeg_ + 1);
          ndrat.n_checks = n_uchecks;
          ndrat.n_singular = n_singular;
          ndrat.setX0(t0_);

          Ret ret = ndrat.reconstruct(uf, mod);

          URationalFunction urat;
          ndrat.writeResult(mod, urat);

          if (ret == SUCCESS) {

            std::size_t num_maxdeg = std::max(xnumdeg_max,
                                              urat.getNumDegree());
            std::size_t den_maxdeg = std::max(xdendeg_max,
                                              urat.getDenDegree());
            std::size_t num_mindeg = std::max(xnumdeg_min,
                                              urat.num().getMinDegree());
            std::size_t den_mindeg = std::max(xdendeg_min,
                                              urat.den().getMinDegree());
            if (xnumdeg_max == urat.getNumDegree())
              num_mindeg = std::min(xnumdeg_min, urat.num().getMinDegree());
            if (xnumdeg_max == urat.getDenDegree())
              den_mindeg = std::min(xdendeg_min, urat.den().getMinDegree());

            xnumdeg_max = num_maxdeg;
            xdendeg_max = den_maxdeg;
            xnumdeg_min = num_mindeg;
            xdendeg_min = den_mindeg;

            okay = true;

          } else {

            ++n_fails;
            if (n_fails >= n_singular)
              return FAILED;

          }

        } while(!okay);
      }
    }

    return SUCCESS;
  }

#define FFLOW_SETUP_HF(hf) \
    hf.n_checks = n_uchecks; \
    hf.n_ndchecks = n_undchecks; \
    hf.n_singular = n_singular; \
    hf.n_learning = std::max({n_uchecks,n_checks,std::size_t(3)}); \
    npolyrec_->n_checks = n_checks; \
    npolyrec_->n_uchecks = n_uchecks;   \
    npolyrec_->n_singular = n_singular; \
    npolyrec_->check_over_maxdeg = false; \
                                       \
    hf.rat_.setX0(t0_); \
    hf.rat_.setMaxTerms(2*maxdeg_ + 2); \
    hf.ndrat_.setX0(t0_); \
    \
    if (num_pref_) \
      hf.f_.num_pref_ = &num_pref_; \
    if (den_pref_)                  \
      hf.f_.den_pref_ = &den_pref_; \
    \
    if (!shift_.empty()) { \
      hf.f_.shift_ = shift_.data(); \
    }

  Ret RatFunReconstruction::degree_(RatFun & f, Mod mod,
                                    unsigned & numdeg, unsigned & dendeg)
  {
    std::unique_ptr<UInt []> x(new UInt[nvars_]);
    UInt res;

    HomogenizeMRatFun hf(f,nvars_);
    FFLOW_SETUP_HF(hf);

    // a few evaluations to get the rank
    const std::size_t n_learning = 1; //hf.n_learning;
    std::size_t n_fails = 0;
    std::unique_ptr<UInt[]> samples(new UInt[nvars_+n_learning+n_singular]);
    for (unsigned i=0; i<nvars_+n_learning+n_singular; ++i)
      samples[i] = sample_uint(OFFSET_1, i, mod);

    for (unsigned ev=0; ev<n_learning; ++ev) {

      bool okay = false;

      do {

        for (unsigned i=0; i<nvars_; ++i)
          x[i] = samples[i+ev+n_fails];

        res = hf.evaluate_learning(x.get(), mod);

        if (res == FAILED) {
          ++n_fails;
          if (n_fails >= n_singular)
            return FAILED;
        } else {
          okay = true;
        }

      } while(!okay);

      if (res == HomogenizeMRatFun::BAD_SHIFT)
        return BAD_SHIFT;

      if (hf.urat().num().getDegree() == hf.numdeg_)
        n0_ = res;
    }
    //hf.set_degree(hf.numdeg_, hf.dendeg_);
    if (VERBOSE) {
      std::cout << "Num. total degree: " << hf.numdeg_ << "\n"
                << "Den. total degree: " << hf.dendeg_ << std::endl;
    }

    numdeg = hf.numdeg_;
    dendeg = hf.dendeg_;

    return SUCCESS;
  }

  Ret RatFunReconstruction::reconstruct_(RatFun & f, Mod mod,
                                         unsigned numdeg, unsigned dendeg)
  {
    std::unique_ptr<UInt []> x(new UInt[nvars_]);

    HomogenizeMRatFun hf(f,nvars_);
    FFLOW_SETUP_HF(hf);

    // rank
    hf.numdeg_ = numdeg;
    hf.dendeg_ = dendeg;
    hf.set_degree(hf.numdeg_, hf.dendeg_);
    hf.init_hornerpoly();

    hf.n0_ = n0_;

    unsigned maxdeg = std::max(hf.numterms(), hf.denterms()) - 1;
    unsigned istart = n0_ == FAILED ? 0 : 1;
    for (unsigned i=istart; i<=maxdeg; ++i) {

      // numerator
      if (i < hf.numterms()) {
        hf.type_ = HomogenizeMRatFun::NUMERATOR;
        npolyrec_->reset();
        npolyrec_->setMinDegree(0);
        npolyrec_->setMaxDegree(i);
        if (xdegs_) {
          npolyrec_->xdegs = xdegs_.num_maxdegs(nvars_)+1;
          if (xdegs_.num_maxdegs(nvars_)[0]<i)
            npolyrec_->setMinDegree(i-xdegs_.num_maxdegs(nvars_)[0]);
        }
        tmp_.clear();

        hf.i_ = i;
        if (npolyrec_->reconstruct(hf,mod) != SUCCESS)
          return FAILED;

        npolyrec_->writeResultVar(1,mod,tmp_);
        if (i == 0) {
          n0_ = tmp_.constTerm();
          hf.n0_ = n0_;
          continue;
        }
        hf.hornernum_[i] = hornerptr_from_sparse_poly(tmp_.data(), tmp_.size(),
                                                      nvars_, 1);
        hf.ww_.ensure_size(horner_required_workspace(hf.hornernum_[i].get()));
        tmp_.homogenize(0,i);

        fun_.numerator().add(std::move(tmp_),mod);
      }

      // denominator
      if (i < hf.denterms()) {
        hf.type_ = HomogenizeMRatFun::DENOMINATOR;
        npolyrec_->reset();
        npolyrec_->setMinDegree(0);
        npolyrec_->setMaxDegree(i);
        if (xdegs_) {
          npolyrec_->xdegs = xdegs_.den_maxdegs(nvars_)+1;
          if (xdegs_.den_maxdegs(nvars_)[0]<i)
            npolyrec_->setMinDegree(i-xdegs_.den_maxdegs(nvars_)[0]);
        }
        tmp_.clear();

        hf.i_ = i;
        if (npolyrec_->reconstruct(hf,mod) != SUCCESS)
          return FAILED;

        npolyrec_->writeResultVar(1,mod,tmp_);
        hf.hornerden_[i] = hornerptr_from_sparse_poly(tmp_.data(), tmp_.size(),
                                                      nvars_, 1);
        hf.ww_.ensure_size(horner_required_workspace(hf.hornerden_[i].get()));
        tmp_.homogenize(0,i);

        fun_.denominator().add(std::move(tmp_),mod);
      }

    }

    if (n0_) {
      tmp_.fromConst(n0_);
      fun_.numerator().add(std::move(tmp_),mod);
    }

    if (!shift_.empty()) {

      std::unique_ptr<UInt[]> minus_shift(new UInt[nvars_]);
      for (unsigned i=0; i<nvars_; ++i)
        minus_shift[i] = neg_mod(shift_[i], mod);

      fun_.numerator().shift(minus_shift.get(), mod);
      fun_.denominator().shift(minus_shift.get(), mod);
      fun_.normalize(mod);

    }

    if (num_pref_)
      fun_.numerator().mul(num_pref_,mod);
    if (den_pref_)
      fun_.denominator().mul(den_pref_,mod);

    return SUCCESS;
  }


  // Fallback to univariate case when nvars_ == 1.

  Ret RatFunReconstruction::reconstruct_univariate_(RatFun & f, Mod mod,
                                                    unsigned numdeg,
                                                    unsigned dendeg)
  {
    URatFunNumDemReconstruction ndrat(numdeg, dendeg);
    ndrat.n_checks = n_uchecks;
    ndrat.n_singular = n_singular;
    ndrat.setX0(t0_);

    MtoURatFun ufun(f);

    Ret ret = ndrat.reconstruct(ufun, mod);
    if (ret == FAILED)
      return FAILED;

    fun_.fromUnivariate(ndrat.getFunction());

    return SUCCESS;
  }

  Ret RatFunReconstruction::reconstruct_univariate_(RatFun & f, Mod mod)
  {
    URatFunReconstruction ndrat(2*maxdeg_ + 1);
    ndrat.n_checks = n_uchecks;
    ndrat.n_singular = n_singular;
    ndrat.setX0(t0_);

    MtoURatFun ufun(f);

    Ret ret = ndrat.reconstruct(ufun, mod);
    if (ret == FAILED)
      return FAILED;

    URationalFunction urat;
    ndrat.writeResult(mod, urat);
    fun_.fromUnivariate(urat);

    return SUCCESS;
  }

  Ret RatFunReconstruction::degree_univariate_(RatFun & f, Mod mod,
                                               unsigned & numdeg,
                                               unsigned & dendeg)
  {
    Ret ret = reconstruct_univariate_(f, mod);
    if (ret != SUCCESS)
      return ret;

    numdeg = fun_.numerator().degree();
    dendeg = fun_.denominator().degree();

    return ret;
  }

  void RatFunReconstruction::sample_univariate_(RatFun & f, Mod mod)
  {
    URatFunReconstruction ndrat(2*maxdeg_ + 1);
    ndrat.n_checks = n_uchecks;
    ndrat.n_singular = n_singular;
    ndrat.setX0(t0_);

    MtoURatFun ufun(f);

    ndrat.sample(ufun, mod);
  }

#define FFLOW_SETUP_HF_SAMPLE(hf)   \
    hf.n_ndchecks = n_undchecks; \
    hf.n_singular = n_singular; \
    npolyrec_->n_checks = n_checks; \
    npolyrec_->n_uchecks = n_uchecks;   \
    npolyrec_->n_singular = n_singular; \
    npolyrec_->check_over_maxdeg = false; \
                                       \
    hf.ndrat_.setX0(t0_); \
    \
    if (!shift_.empty()) { \
      hf.f_.shift_ = shift_.data(); \
    }

  void RatFunReconstruction::sample_(RatFun & f, Mod mod,
                                     unsigned numdeg, unsigned dendeg)
  {
    std::unique_ptr<UInt []> x(new UInt[nvars_]);

    HomogenizeMRatFunSample hf(f,nvars_);
    FFLOW_SETUP_HF_SAMPLE(hf);

    // rank
    hf.numdeg_ = numdeg;
    hf.dendeg_ = dendeg;
    hf.set_degree(hf.numdeg_, hf.dendeg_);

    unsigned maxdeg = std::max(hf.numterms(), hf.denterms()) - 1;
    unsigned istart = n0_ == FAILED ? 0 : 1;
    for (unsigned i=istart; i<=maxdeg; ++i) {

      // numerator
      if (i < hf.numterms()) {
        hf.type_ = HomogenizeMRatFun::NUMERATOR;
        npolyrec_->reset();
        npolyrec_->setMinDegree(0);
        npolyrec_->setMaxDegree(i);
        if (xdegs_) {
          npolyrec_->xdegs = xdegs_.num_maxdegs(nvars_)+1;
          if (xdegs_.num_maxdegs(nvars_)[0]<i)
            npolyrec_->setMinDegree(i-xdegs_.num_maxdegs(nvars_)[0]);
        }

        hf.i_ = i;
        npolyrec_->sample(hf,mod);
      }

      // denominator
      if (i < hf.denterms()) {
        hf.type_ = HomogenizeMRatFun::DENOMINATOR;
        npolyrec_->reset();
        npolyrec_->setMinDegree(0);
        npolyrec_->setMaxDegree(i);
        if (xdegs_) {
          npolyrec_->xdegs = xdegs_.den_maxdegs(nvars_)+1;
          if (xdegs_.den_maxdegs(nvars_)[0]<i)
            npolyrec_->setMinDegree(i-xdegs_.den_maxdegs(nvars_)[0]);
        }

        hf.i_ = i;
        npolyrec_->sample(hf,mod);
      }

    }

  }

} // namespace fflow
