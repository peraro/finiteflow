#include <fflow/mp_multivariate_reconstruction.hh>
#include <fflow/mp_gcd.hh>
#include <fflow/function_cache.hh>


namespace fflow {

  Ret MPRatFunReconstruction::degree(RatFun & f, Mod mod,
                                     unsigned & numdeg, unsigned & dendeg)
  {
    RatFunReconstruction rec(nvars_, maxdeg_);
    setup_rat_rec(rec);
    return rec.degree(f, mod, numdeg, dendeg);
  }


  Ret MPRatFunReconstruction::reconstruct_(RatFun & f,
                                           unsigned numdeg, unsigned dendeg,
                                           const RatFunVarDegrees * degs,
                                           MPReconstructedRatFun & recf)
  {
    MPInt totmod, c1, c2;
    SparseRationalFunction tmp;

    std::unique_ptr<UInt[]> x(new UInt[nvars_]);

    for (unsigned mod_i=0; mod_i<max_primes+1; ++mod_i) {

      Mod mod(BIG_UINT_PRIMES[mod_i % BIG_UINT_PRIMES_SIZE]);

      if (mod_i) {

        bool testing_ok = true;

        for (unsigned check_i=0; check_i<n_checks; ++check_i) {

          // evaluate the function
          UInt res;
          for (unsigned ev=0; ev<=n_singular; ++ev) {

            for (unsigned i=0; i<nvars_; ++i)
              x[i] = sample_uint(OFFSET_1, check_i + ev + i, mod);

            if ((res = f.evaluate(x.get(), mod)) != FAILED)
              break;
            else if (ev == n_singular)
              return FAILED;
          }

          if (recf.eval(x.get(), mod) != res) {
            testing_ok = false;
            break;
          }

        }

        if (testing_ok)
          return SUCCESS;

      }

      if (mod_i >= max_primes)
        return MISSING_PRIMES;

      if (mod_i == 0)
        totmod = mod.n();
      else
        chinese_remainder_coeffs(totmod, mod.n(), c1, c2, totmod);

      RatFunReconstruction rec(nvars_, maxdeg_);
      setup_rat_rec(rec);

      Ret ret = 0;
      if (degs == nullptr)
        ret = rec.reconstruct(f, mod);
      else
        ret = rec.reconstructWithDegs(f, mod, numdeg, dendeg, *degs);
      if (ret != SUCCESS)
        return ret;

      tmp = std::move(rec.getFunction());

      if (mod_i) {
        recf.merge(std::move(tmp),c1,c2,totmod);
      } else {
        recf.setMod(totmod);
        recf.copy(std::move(tmp));
      }
    }

    return FAILED;
  }

  void MPRatFunReconstruction::sample(RatFun & f,
                                      unsigned numdeg, unsigned dendeg,
                                      const RatFunVarDegrees & degs)
  {
    std::unique_ptr<UInt[]> x(new UInt[nvars_]);

    for (unsigned mod_i=0; mod_i<max_primes+1; ++mod_i) {

      Mod mod(BIG_UINT_PRIMES[mod_i % BIG_UINT_PRIMES_SIZE]);

      if (mod_i) {

        for (unsigned check_i=0; check_i<n_checks+n_singular; ++check_i) {

          for (unsigned i=0; i<nvars_; ++i)
            x[i] = sample_uint(OFFSET_1, check_i + i, mod);
          f.evaluate(x.get(), mod);

        }

      }

      if (mod_i >= max_primes)
        return;

      RatFunReconstruction rec(nvars_, maxdeg_);
      setup_rat_rec(rec);

      rec.sample(f, mod, numdeg, dendeg, degs);
    }
  }

} // namespace fflow
