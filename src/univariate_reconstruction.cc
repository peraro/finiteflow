#include <fflow/univariate_reconstruction.hh>
#include <fflow/gcd.hh>

namespace fflow {

  Ret UPolyReconstruction::reconstruct(URatFun & f, Mod mod)
  {
    UInt x, y=0, t0, t1;
    unsigned check_count=0;
    unsigned mindeg = poly_.getDegree();

    for (unsigned i=0; i<=maxdeg_+n_checks; ++i) {

      UInt a;
      UInt * ai = poly_.coeffs();

      x = i ? add_mod(getXi(i-1), SAMPLING_STRIDE, mod) : getX0();

      for (unsigned ev=0; ev<=n_singular; ++ev) {
        if ((y = f.evaluate(x,mod)) != FAILED)
          break;
        else if (ev == n_singular)
          return FAILED;
        x = add_mod(x, SAMPLING_STRIDE, mod);
      }

      t0 = poly_.evalNewtonSystemRow(x, mod, i, t1);
      a = div_mod(sub_mod(y,t0,mod), t1, mod);

      if (a == 0 && i > mindeg) {
        ++check_count;
        if (check_count >= n_checks)
          return SUCCESS;
      } else {
        check_count = 0;
      }

      if (poly_.getDegree() < i) {
        poly_.push_back(x, a);
      } else {
        poly_.setXi(i,x);
        ai[i] = a;
      }

      if (!check_over_maxdeg && i == maxdeg_)
        return SUCCESS;

      if (a !=0 && poly_.getDegree()>maxdeg_)
        return FAILED;

    }

    return FAILED;
  }

  void UPolyReconstruction::sample(URatFun & f, Mod mod)
  {
    UInt x=getX0();
    std::size_t n_samples = maxdeg_ + 1 + n_singular
      + (check_over_maxdeg ? n_checks : 0);
    for (unsigned i=0; i<n_samples; ++i) {
      f.evaluate(x, mod);
      x = add_mod(x, SAMPLING_STRIDE, mod);
    }
  }

  unsigned UPolyReconstruction::getTrueDegree()
  {
    int deg = poly_.getDegree();
    UInt * ai = poly_.coeffs();
    while ((deg > 0) && (ai[deg] == 0))
      --deg;
    return deg;
  }


  Ret URatFunReconstruction::reconstruct(URatFun & f, Mod mod)
  {
    UInt x, y, yrec;
    unsigned check_count = 0;
    //unsigned minterms = fun_.size();

    for (unsigned i=0; i<maxterms_+n_checks; ++i) {

      UInt a = 0;
      UInt * ai = fun_.coeffs();

      x = i ? add_mod(getXi(i-1), SAMPLING_STRIDE, mod) : getX0();

      for (unsigned ev=0; ev<=n_singular; ++ev) {

        bool term_ok = true;

        if ((y = f.evaluate(x,mod)) == FAILED)
          term_ok = false;

        yrec = fun_.safe_eval(x,mod);
        if (yrec == FAILED)
          term_ok = false;

        if (term_ok && yrec == y) {

          check_count = 1;
          while (check_count < n_checks) {

            for (unsigned ev2=0; ev2<=n_singular; ++ev2) {
              x = add_mod(x, SAMPLING_STRIDE, mod);
              if ((y = f.evaluate(x,mod)) != FAILED)
                break;
              else if (ev2==n_singular)
                return FAILED;
            }

            if (fun_.safe_eval(x,mod) == y) {
              ++ check_count;
            } else {
              check_count = 0;
              term_ok = false;
              break;
            }

          }

          if (check_count >= n_checks) {
            fun_.remove_high_zeroes();
            return SUCCESS;
          }

        }

        if (term_ok) {

          a = fun_.evalThieleSystemCoeff(x, y, mod, i);

          if (a == FAILED || a == 0)
            term_ok = false;

        }

        if (!term_ok) {
          if (ev == n_singular)
            return FAILED;
          x = add_mod(x, 1, mod);
        } else {
          break;
        }
      }

      if (fun_.size()>=maxterms_)
        return FAILED;

      if (fun_.size() <= i) {
        fun_.push_back(x,a);
      } else {
        fun_.setXi(i,x);
        ai[i] = a;
      }

    }

    return FAILED;
  }

  void URatFunReconstruction::sample(URatFun & f, Mod mod)
  {
    UInt x = getX0();
    for (unsigned i=0; i<maxterms_+n_checks+n_singular; ++i) {
      f.evaluate(x,mod);
      x = add_mod(x, SAMPLING_STRIDE, mod);
    }
  }



  Ret URatFunNumDemReconstruction::reconstruct(URatFun & f, Mod mod)
  {
    const std::size_t num_coeffs = fun_.getNumDegree() + 1;
    const std::size_t den_coeffs = fun_.getDenDegree();
    const std::size_t n_coeffs = num_coeffs + den_coeffs;
    const std::size_t n_eqs = n_coeffs + n_checks;

    // We have:
    //
    // (n[0] + n[1]*t + n[2]*t^2 + ...)/(1 + d[1]*t + d[2]*t^2 + ...) = f(t)
    //
    // so the equations for (n[0],n[1],...,d[1],d[2],...) are
    //
    // n[0] + n[1]*t + n[2]*t^2 + ... = d[1]*t*f(t) - d[2]*t^2*f(t) + ... + f(t)

    if (! mat_.isDimension(n_eqs,n_coeffs+1))
      mat_.reset(n_eqs,n_coeffs+1);

    UInt x = t0_;
    for (std::size_t i=0; i<n_eqs; ++i) {

      UInt fres=0;

      for (unsigned ev = 0; ev<=n_singular; ++ev) {
        if ((fres = f.evaluate(x,mod)) != FAILED)
          break;
        else if (ev == n_singular)
          return FAILED;
        x = add_mod(x, SAMPLING_STRIDE, mod);
      }

      // fun. value
      mat_(i,n_coeffs) = fres;

      // numerator
      UInt tpow = 1;
      for (std::size_t j=0; j<num_coeffs; ++j) {
        mat_(i,j) = tpow;
        tpow = mul_mod(tpow, x, mod);
      }

      // denominator
      tpow = mul_mod(fres, x, mod);
      for (std::size_t j=1; j<den_coeffs+1; ++j) {
        mat_(i,num_coeffs+j-1) = tpow;
        tpow = mul_mod(tpow, x, mod);
      }

      x = add_mod(x, SAMPLING_STRIDE, mod);
    }

    MatrixView smv = mat_.getMatrixView();
    smv.toReducedRowEcholon(mod);

    if (smv.isZeroDimSystem()) {

      UInt * c = fun_.num().coeffs();
      for (std::size_t j=0; j<num_coeffs; ++j)
        c[j] = mat_(j,n_coeffs);

      c = fun_.den().coeffs();
      c[0] = 1;
      for (std::size_t j=1; j<den_coeffs+1; ++j)
        c[j] = neg_mod(mat_(num_coeffs+j-1,n_coeffs), mod);

      return SUCCESS;
    }

    return FAILED;
  }


  void URatFunNumDemReconstruction::sample(URatFun & f, Mod mod)
  {
    const std::size_t num_coeffs = fun_.getNumDegree() + 1;
    const std::size_t den_coeffs = fun_.getDenDegree();
    const std::size_t n_coeffs = num_coeffs + den_coeffs;
    const std::size_t n_eqs = n_coeffs + n_checks + n_singular;

    UInt x = t0_;
    for (std::size_t i=0; i<n_eqs; ++i) {
      f.evaluate(x,mod);
      x = add_mod(x, SAMPLING_STRIDE, mod);
    }
  }


  Ret URatFunNumDenRecoHighDegs::reconstruct(URatFun & f, Mod mod)
  {
    const std::size_t num_deg = fun_.getNumDegree();
    const std::size_t den_deg = fun_.getDenDegree();

    const std::size_t num_unknowns = fun_.getNumDegree() + 1 - rnum_min;
    const std::size_t den_unknowns = fun_.getDenDegree() + 1 - rden_min;

    const std::size_t n_unknowns = num_unknowns + den_unknowns;
    const std::size_t n_eqs = n_unknowns + n_checks;

    if (! mat_.isDimension(n_eqs,n_unknowns+1))
      mat_.reset(n_eqs,n_unknowns+1);

    UInt x = t0_;
    for (std::size_t i=0; i<n_eqs; ++i) {

      UInt fres=0;

      for (unsigned ev = 0; ev<=n_singular; ++ev) {
        if ((fres = f.evaluate(x,mod)) != FAILED)
          break;
        else if (ev == n_singular)
          return FAILED;
        x = add_mod(x, SAMPLING_STRIDE, mod);
      }

      // numerator
      unsigned col = 0;
      UInt fres_sub = 0;
      UInt tpow = 1;
      UInt * c_known = fun_.num().coeffs();
      col = 0;
      for (std::size_t j=0; j<=num_deg; ++j) {
        if (j<rnum_min)
          addmul_mod(fres_sub, tpow, c_known[j], mod);
        else
          mat_(i,col++) = tpow;
        tpow = mul_mod(tpow, x, mod);
      }

      // denominator
      UInt fres_c = fres;
      tpow = mul_mod(fres, x, mod);
      c_known = fun_.den().coeffs();
      for (std::size_t j=1; j<=den_deg; ++j) {
        if (j<rden_min)
          addmul_mod(fres_c, tpow, c_known[j], mod);
        else
          mat_(i,col++) = tpow;
        tpow = mul_mod(tpow, x, mod);
      }

      mat_(i,n_unknowns) = sub_mod(fres_c, fres_sub, mod);

      x = add_mod(x, SAMPLING_STRIDE, mod);
    }

    MatrixView smv = mat_.getMatrixView();
    smv.toReducedRowEcholon(mod);

    if (smv.isZeroDimSystem()) {

      unsigned var=0;
      UInt * c = fun_.num().coeffs();
      for (std::size_t j=rnum_min; j<=num_deg; ++j)
        c[j] = mat_(var++,n_unknowns);

      c = fun_.den().coeffs();
      c[0] = 1;
      for (std::size_t j=rden_min; j<=den_deg; ++j)
        c[j] = neg_mod(mat_(var++,n_unknowns), mod);

      return SUCCESS;
    }

    return FAILED;
  }

  void URatFunNumDenRecoHighDegs::sample(URatFun & f, Mod mod)
  {
    const std::size_t num_unknowns = fun_.getNumDegree() + 1 - rnum_min;
    const std::size_t den_unknowns = fun_.getDenDegree() + 1 - rden_min;

    const std::size_t n_unknowns = num_unknowns + den_unknowns;
    const std::size_t n_eqs = n_unknowns + n_checks + n_singular;

    UInt x = t0_;
    for (std::size_t i=0; i<n_eqs; ++i) {
      f.evaluate(x,mod);
      x = add_mod(x, SAMPLING_STRIDE, mod);
    }
  }


} // namespace fflow
