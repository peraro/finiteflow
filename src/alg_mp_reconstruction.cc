#include <fstream>
#include <fflow/alg_mp_reconstruction.hh>
#include <fflow/mp_gcd.hh>

namespace fflow {

  namespace  {

    Ret algorithm_get_degree(AlgorithmRatFun & fun,
                             const UInt shift[],
                             const ReconstructionOptions & opt,
                             unsigned & numdeg, unsigned & dendeg)
    {
      unsigned nparsin = fun.get_alg()->nparsin[0];
      RatFunReconstruction funrec(nparsin);
      set_rec_options(opt, shift, funrec);

      Mod mod(BIG_UINT_PRIMES[opt.start_mod]);
      unsigned ndeg, ddeg;
      Ret ret = funrec.degree(fun, mod, ndeg, ddeg);

      numdeg = ndeg;
      dendeg = ddeg;

      return ret;
    }

    Ret algorithm_rec_univariate(AlgorithmRatFun & fun,
                                 const UInt shift[],
                                 const ReconstructionOptions & opt,
                                 MPReconstructedRatFun & res)
    {
      unsigned nparsin = fun.get_alg()->nparsin[0];
      if (nparsin != 1)
        return FAILED;

      MPRatFunReconstruction funrec(nparsin);
      set_mp_rec_options(opt, shift, funrec);

      Ret ret = funrec.reconstruct(fun, res);

      return ret;
    }

    Ret algorithm_rec_univariate_mod(Mod mod,
                                     AlgorithmRatFun & fun,
                                     const UInt shift[],
                                     const ReconstructionOptions & opt,
                                     MPReconstructedRatFun & res)
    {
      unsigned nparsin = fun.get_alg()->nparsin[0];
      if (nparsin != 1)
        return FAILED;

      MPRatFunReconstruction funrec(nparsin);
      set_mp_rec_options(opt, shift, funrec);

      RatFunReconstruction rec(nparsin, RatFunReconstruction::DEFAULT_MAX_DEG);
      funrec.setup_rat_rec(rec);

      Ret ret = rec.reconstruct(fun, mod);

      if (ret == SUCCESS) {
        SparseRationalFunction tmp = std::move(rec.getFunction());
        res.copy_exactly(std::move(tmp));
      }

      return ret;
    }

  } // namespace


  Ret algorithm_get_degrees(Algorithm & alg,
                            AlgorithmData * data,
                            Context * ctxt,
                            const UInt shift[],
                            const ReconstructionOptions & opt,
                            unsigned numdeg[], unsigned dendeg[])
  {
    unsigned nparsout = alg.nparsout;

    AlgorithmRatFun fun;
    fun.set_algorithm(alg);
    fun.set_algorithm_data(data);
    fun.set_context(ctxt);

    for (unsigned idx=0; idx<nparsout; ++idx) {
      fun.idx = idx;
      Ret ret = algorithm_get_degree(fun, shift, opt,
                                     numdeg[idx], dendeg[idx]);
      if (ret != SUCCESS)
        return ret;
    }
    return SUCCESS;
  }


  Ret algorithm_dump_degree_info(const char * filename,
                                 unsigned npars_in,
                                 unsigned npars_out,
                                 const unsigned numdeg[],
                                 const unsigned dendeg[],
                                 const RatFunVarDegrees degs[])
  {
    std::ofstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto dump = [&file](UInt z)
      {
        file.write(reinterpret_cast<char*>(&z), sizeof(UInt));
      };

    dump(npars_in);
    dump(npars_out);

    for (unsigned j=0; j<npars_out; ++j) {
      dump(numdeg[j]);
      dump(dendeg[j]);
      for (unsigned v=0; v<npars_in; ++v) {
        dump(degs[j].num_maxdegs(npars_in)[v]);
        dump(degs[j].num_mindegs(npars_in)[v]);
        dump(degs[j].den_maxdegs(npars_in)[v]);
        dump(degs[j].den_mindegs(npars_in)[v]);
      }
    }

    return SUCCESS;
  }


  Ret algorithm_load_degree_info(const char * filename,
                                 unsigned npars_in,
                                 unsigned npars_out,
                                 unsigned numdeg[],
                                 unsigned dendeg[],
                                 RatFunVarDegrees degs[])
  {
    std::ifstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto load = [&file]() -> UInt
      {
        UInt z;
        file.read(reinterpret_cast<char*>(&z), sizeof(UInt));
        return z;
      };
#define LOAD(z)                                 \
    do {                                        \
      z = load();                               \
      if (file.fail())                          \
        return FAILED;                          \
    } while (0)

    UInt check_npars_in, check_npars_out;

    check_npars_in = load();
    if (file.fail() || check_npars_in != npars_in)
      return FAILED;

    check_npars_out = load();
    if (file.fail() || check_npars_out != npars_out)
      return FAILED;

    for (unsigned j=0; j<npars_out; ++j) {
      LOAD(numdeg[j]);
      LOAD(dendeg[j]);
      degs[j].resize(npars_in);
      for (unsigned v=0; v<npars_in; ++v) {
        LOAD(degs[j].num_maxdegs(npars_in)[v]);
        LOAD(degs[j].num_mindegs(npars_in)[v]);
        LOAD(degs[j].den_maxdegs(npars_in)[v]);
        LOAD(degs[j].den_mindegs(npars_in)[v]);
      }
    }

#undef LOAD

    return SUCCESS;
  }


  Ret algorithm_npars_from_degree_info(const char * filename,
                                       unsigned & npars_in,
                                       unsigned & npars_out)
  {
    std::ifstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto load = [&file]() -> UInt
      {
        UInt z;
        file.read(reinterpret_cast<char*>(&z), sizeof(UInt));
        return z;
      };

    npars_in = load();
    if (file.fail())
      return FAILED;

    npars_out = load();
    if (file.fail())
      return FAILED;

    return SUCCESS;
  }


  Ret algorithm_load_degree_info_for_var(const char * filename,
                                         unsigned npars_in,
                                         unsigned npars_out,
                                         unsigned var,
                                         unsigned & numdeg,
                                         unsigned & dendeg,
                                         RatFunVarDegrees & degs)
  {
    std::ifstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto load = [&file]() -> UInt
      {
        UInt z;
        file.read(reinterpret_cast<char*>(&z), sizeof(UInt));
        return z;
      };
#define LOAD(z)                                 \
    do {                                        \
      z = load();                               \
      if (file.fail())                          \
        return FAILED;                          \
    } while (0)

    UInt check_npars_in, check_npars_out;

    check_npars_in = load();
    if (file.fail() || check_npars_in != npars_in)
      return FAILED;

    check_npars_out = load();
    if (file.fail() || check_npars_out != npars_out)
      return FAILED;

    file.ignore(var*(2+4*npars_in)*sizeof(UInt));
    if (file.fail())
      return FAILED;

    LOAD(numdeg);
    LOAD(dendeg);

    degs.resize(npars_in);

    for (unsigned v=0; v<npars_in; ++v) {
      LOAD(degs.num_maxdegs(npars_in)[v]);
      LOAD(degs.num_mindegs(npars_in)[v]);
      LOAD(degs.den_maxdegs(npars_in)[v]);
      LOAD(degs.den_mindegs(npars_in)[v]);
    }

#undef LOAD

    return SUCCESS;
  }


  Ret algorithm_reconstruct_univariate(const Algorithm & alg,
                                       AlgorithmData * data,
                                       Context * ctxt,
                                       const UInt shift[],
                                       const ReconstructionOptions & opt,
                                       MPReconstructedRatFun res[])
  {
    unsigned nparsout = alg.nparsout;

    AlgorithmRatFun fun;
    fun.set_algorithm(alg);
    fun.set_algorithm_data(data);
    fun.set_context(ctxt);

    for (unsigned idx=0; idx<nparsout; ++idx) {
      fun.idx = idx;
      Ret ret = algorithm_rec_univariate(fun, shift, opt, res[idx]);
      if (ret != SUCCESS)
        return ret;
    }
    return SUCCESS;
  }

  Ret algorithm_reconstruct_univariate_mod(Mod mod,
                                           const Algorithm & alg,
                                           AlgorithmData * data,
                                           Context * ctxt,
                                           const UInt shift[],
                                           const ReconstructionOptions & opt,
                                           MPReconstructedRatFun res[])
  {
    unsigned nparsout = alg.nparsout;

    AlgorithmRatFun fun;
    fun.set_algorithm(alg);
    fun.set_algorithm_data(data);
    fun.set_context(ctxt);

    for (unsigned idx=0; idx<nparsout; ++idx) {
      fun.idx = idx;
      Ret ret = algorithm_rec_univariate_mod(mod, fun, shift, opt, res[idx]);
      if (ret != SUCCESS)
        return ret;
    }
    return SUCCESS;
  }


  namespace  {

    Ret algorithm_get_var_degree(AlgorithmRatFun & fun,
                                 unsigned var,
                                 const UInt shift[],
                                 const ReconstructionOptions & opt,
                                 std::size_t & num_mindeg,
                                 std::size_t & num_maxdeg,
                                 std::size_t & den_mindeg,
                                 std::size_t & den_maxdeg)
    {
      unsigned nparsin = fun.get_alg()->nparsin[0];
      RatFunReconstruction funrec(nparsin);
      set_rec_options(opt, shift, funrec);

      Mod mod(BIG_UINT_PRIMES[opt.start_mod]);
      std::size_t nmindeg, nmaxdeg, dmindeg, dmaxdeg;
      Ret ret = funrec.var_degree(fun, var, mod,
                                  nmindeg, nmaxdeg, dmindeg, dmaxdeg);

      num_mindeg = nmindeg;
      num_maxdeg = nmaxdeg;
      den_mindeg = dmindeg;
      den_maxdeg = dmaxdeg;

      return ret;
    }

  } // namespace


  Ret algorithm_get_var_degrees(Algorithm & alg,
                                AlgorithmData * data,
                                Context * ctxt,
                                unsigned var,
                                const UInt shift[],
                                const ReconstructionOptions & opt,
                                RatFunVarDegrees degs[])
  {
    unsigned nparsin = alg.nparsin[0], nparsout = alg.nparsout;

    AlgorithmRatFun fun;
    fun.set_algorithm(alg);
    fun.set_algorithm_data(data);
    fun.set_context(ctxt);

    for (unsigned idx=0; idx<nparsout; ++idx) {
      fun.idx = idx;
      Ret ret = algorithm_get_var_degree(fun, var,
                                         shift, opt,
                                         degs[idx].num_mindegs(nparsin)[var],
                                         degs[idx].num_maxdegs(nparsin)[var],
                                         degs[idx].den_mindegs(nparsin)[var],
                                         degs[idx].den_maxdegs(nparsin)[var]);
      if (ret != SUCCESS)
        return FAILED;
    }

    return SUCCESS;
  }


  namespace  {

    void algorithm_merge_degrees(const RatFunVarDegrees degs[],
                                 const unsigned numdegs[],
                                 const unsigned dendegs[],
                                 unsigned nparsin, unsigned nparsout,
                                 RatFunVarDegrees & mdegs,
                                 unsigned & mnumdeg,
                                 unsigned & mdendeg)
    {
      unsigned nv = nparsin;
      mdegs.resize(nv);

      for (unsigned idx=0; idx<nparsout; ++idx)
        for (unsigned v=0; v<nv; ++v) {
          mdegs.num_maxdegs(nv)[v] = std::max(mdegs.num_maxdegs(nv)[v],
                                              degs[idx].num_maxdegs(nv)[v]
                                              - degs[idx].num_mindegs(nv)[v]);
          mdegs.den_maxdegs(nv)[v] = std::max(mdegs.den_maxdegs(nv)[v],
                                              degs[idx].den_maxdegs(nv)[v]
                                              - degs[idx].den_mindegs(nv)[v]);
        }

      if (!nparsout) {
        mnumdeg = mdendeg = 0;
      } else {
        mnumdeg = *std::max_element(numdegs, numdegs + nparsout);
        mdendeg = *std::max_element(dendegs, dendegs + nparsout);
      }
    }


    void algorithm_get_sample_points(GenerateSamplePoints & alg,
                                     unsigned nparsin, unsigned,
                                     const UInt shift[],
                                     const ReconstructionOptions & opt,
                                     const RatFunVarDegrees & mdegs,
                                     unsigned mnumdeg, unsigned mdendeg)
    {
      MPRatFunReconstruction funrec(nparsin);
      set_mp_rec_options(opt, shift, funrec);
      funrec.sample(alg, mnumdeg, mdendeg, mdegs);
    }

    void algorithm_verify_sample_points_idx(VerifySamplePoints & alg,
                                            unsigned nparsin, unsigned,
                                            const UInt shift[],
                                            const ReconstructionOptions & opt,
                                            const RatFunVarDegrees & mdegs,
                                            unsigned mnumdeg, unsigned mdendeg,
                                            unsigned idx)
    {
      MPRatFunReconstruction funrec(nparsin);
      set_mp_rec_options(opt, shift, funrec);
      alg.idx = idx;
      funrec.sample(alg, mnumdeg, mdendeg, mdegs);
    }

  } // namespace


  void algorithm_generate_sample_points(GenerateSamplePoints & samples,
                                        unsigned nparsin, unsigned nparsout,
                                        const UInt shift[],
                                        const ReconstructionOptions & opt,
                                        const RatFunVarDegrees degs[],
                                        const unsigned numdegs[],
                                        const unsigned dendegs[])
  {
    RatFunVarDegrees mdegs;
    unsigned nmdeg=0, dmdeg=0;
    if (nparsin > 1)
      algorithm_merge_degrees(degs, numdegs, dendegs, nparsin, nparsout,
                              mdegs, nmdeg, dmdeg);

    algorithm_get_sample_points(samples, nparsin, nparsout,
                                shift, opt, mdegs, nmdeg, dmdeg);
  }


  void algorithm_verify_sample_points(std::unique_ptr<UInt[]> * samples,
                                      unsigned nsamples,
                                      unsigned nparsin, unsigned nparsout,
                                      const std::vector<unsigned> & needed,
                                      const UInt shift[],
                                      const ReconstructionOptions & opt,
                                      const RatFunVarDegrees degs[],
                                      const unsigned numdegs[],
                                      const unsigned dendegs[])
  {
    VerifySamplePoints alg(nparsin);
    UIntCache positions;
    positions.init(nparsin+1);
    positions.reserve_more(nsamples, nsamples*(nparsin+1+1));
    for (unsigned j=0; j<nsamples; ++j) {
      UInt * xin = positions.get_new_inptr();
      std::copy(samples[j].get(), samples[j].get()+nparsin+1, xin);
      UInt * xout = positions.new_entry(xin, 1);
      xout[0] = j;
    }
    alg.set_samples(samples, &positions);
    if (nparsin>1) {
      for (auto idx : needed) {
        algorithm_verify_sample_points_idx(alg, nparsin, nparsout,
                                           shift, opt, degs[idx],
                                           numdegs[idx], dendegs[idx],
                                           idx);
      }
    } else {
      RatFunVarDegrees mdegs;
      for (auto idx : needed) {
        algorithm_verify_sample_points_idx(alg, nparsin, nparsout,
                                           shift, opt, mdegs,
                                           0, 0, idx);
      }
    }
  }


  namespace {

    Ret algorithm_get_rec_fail_status(const std::vector<UInt> & missing_mods,
                                      const ReconstructionOptions & opt)
    {
      if (missing_mods.size() == 0)
        return FAILED;

      if (missing_mods.size() == 1 &&
          missing_mods[0] == BIG_UINT_PRIMES[opt.start_mod])
        return MISSING_SAMPLES;

      return MISSING_PRIMES;
    }

  } // namespace



  Ret algorithm_reconstruct(const UIntCache & cache,
                            unsigned nparsin, unsigned,
                            unsigned idx,
                            const UInt shift[],
                            const ReconstructionOptions & opt,
                            const RatFunVarDegrees degs[],
                            const unsigned numdeg[], const unsigned dendeg[],
                            MPReconstructedRatFun & recf)
  {
    SampledRatFun sampledfun;
    sampledfun.set_function_cache(&cache);
    sampledfun.idx = idx;

    MPRatFunReconstruction funrec(nparsin);
    set_mp_rec_options(opt, shift, funrec);

    Ret ret = 0;

    if (nparsin != 1)
      ret = funrec.reconstructWithDegs(sampledfun,
                                       numdeg[idx], dendeg[idx], degs[idx],
                                       recf);
    else
      ret = funrec.reconstruct(sampledfun, recf);

    if (ret == FAILED)
      ret = algorithm_get_rec_fail_status(sampledfun.missing_mods(), opt);

    return ret;
  }


  Ret algorithm_sparse_reconstruct(const UIntCache & cache,
                                   unsigned nparsin, unsigned nparsout,
                                   unsigned idx,
                                   const UInt shift[],
                                   const ReconstructionOptions & opt,
                                   const RatFunVarDegrees degs[],
                                   const unsigned numdeg[],
                                   const unsigned dendeg[],
                                   MPReconstructedRatFun & recf)
  {
    SparselySampledRatFun sampledfun;
    sampledfun.set_function_cache(&cache, nparsin, nparsout);
    sampledfun.idx = idx;

    MPRatFunReconstruction funrec(nparsin);
    set_mp_rec_options(opt, shift, funrec);

    Ret ret = 0;

    if (nparsin != 1)
      ret = funrec.reconstructWithDegs(sampledfun,
                                       numdeg[idx], dendeg[idx], degs[idx],
                                       recf);
    else
      ret = funrec.reconstruct(sampledfun, recf);

    if (ret == FAILED)
      ret = algorithm_get_rec_fail_status(sampledfun.missing_mods(), opt);

    return ret;
  }


  Ret algorithm_sparse_reconstruct_mod(Mod mod,
                                       const UIntCache & cache,
                                       unsigned nparsin, unsigned nparsout,
                                       unsigned idx,
                                       const UInt shift[],
                                       const ReconstructionOptions & opt,
                                       const RatFunVarDegrees degs[],
                                       const unsigned numdeg[],
                                       const unsigned dendeg[],
                                       MPReconstructedRatFun & recf)
  {
    SparselySampledRatFun sampledfun;
    sampledfun.set_function_cache(&cache, nparsin, nparsout);
    sampledfun.idx = idx;

    MPRatFunReconstruction funrec(nparsin);
    set_mp_rec_options(opt, shift, funrec);

    RatFunReconstruction rec(nparsin, RatFunReconstruction::DEFAULT_MAX_DEG);
    funrec.setup_rat_rec(rec);

    Ret ret = 0;

    if (nparsin != 1)
      ret = rec.reconstructWithDegs(sampledfun, mod,
                                    numdeg[idx], dendeg[idx], degs[idx]);
    else
      ret = rec.reconstruct(sampledfun, mod);

    if (ret == SUCCESS) {
      SparseRationalFunction tmp = std::move(rec.getFunction());
      recf.copy_exactly(std::move(tmp));
    }

    return ret;
  }


  Ret algorithm_reconstruct_from_evals(const std::string files[],
                                       unsigned nfiles,
                                       const char * degs_file,
                                       unsigned nparsin, unsigned nparsout,
                                       unsigned idx,
                                       const UInt shift[],
                                       const ReconstructionOptions & opt,
                                       MPReconstructedRatFun & recf)
  {
    UIntCache cache;
    cache.init(nparsin+1);

    unsigned numdeg, dendeg;
    RatFunVarDegrees vdegs;

    Ret ret = algorithm_load_degree_info_for_var(degs_file, nparsin, nparsout,
                                                 idx, numdeg, dendeg, vdegs);
    if (ret == FAILED)
      return FAILED;

    for (unsigned j=0; j<nfiles; ++j) {
      ret = load_evaluations_for_var(files[j].c_str(), nparsin, nparsout,
                                     cache, idx);
      if (ret == FAILED)
        return FAILED;
    }

    ret = algorithm_reconstruct(cache, nparsin, 1, 0,
                                shift, opt, &vdegs, &numdeg, &dendeg, recf);

    return ret;
  }


  Ret algorithm_reconstruct_numeric(const Algorithm & alg,
                                    AlgorithmData * data,
                                    Context * ctxt,
                                    const ReconstructionOptions & opt,
                                    MPRational res[])
  {
    std::unique_ptr<UInt[]> xout(new UInt[alg.nparsout]);
    std::unique_ptr<MPInt[]> zout(new MPInt[alg.nparsout]);
    MPInt totmod, c1, c2;

    const std::size_t max_primes = opt.max_primes;
    std::size_t prime_i = opt.start_mod;
    unsigned singular = 0;
    const UInt * xin = nullptr;

    for (std::size_t ip=0; ip<max_primes+1; ++ip) {

      Mod mod(BIG_UINT_PRIMES[prime_i]);
      Ret ret = alg.evaluate(ctxt, &xin, mod, data, xout.get());
      if (ret != SUCCESS) {
        if (singular == opt.n_singular)
          return FAILED;
        ++singular;
        continue;
      }

      if (ip == 0) {

        totmod = mod.n();
        for (unsigned j=0; j<alg.nparsout; ++j) {
          if (xout[j]) {
            zout[j] = xout[j];
            res[j] = rat_rec(xout[j], mod.n());
          } else {
            zout[j] = UInt(0);
            res[j] = Int(0);
          }
        }

      } else {

        // check
        MPInt this_mod(mod.n()), xmp;
        bool test_success = true;
        for (unsigned j=0; j<alg.nparsout; ++j) {
          if (!rat_mod(res[j], this_mod, xmp))
            return FAILED;
          if (xout[j] != xmp.to_uint()) {
            test_success = false;
            break;
          }
        }

        if (test_success)
          return SUCCESS;

        if (ip==max_primes)
          return FAILED;

        // merge
        chinese_remainder_coeffs(totmod, mod.n(), c1, c2, totmod);
        for (unsigned j=0; j<alg.nparsout; ++j) {
          chinese_remainder_from_coeffs(zout[j], xout[j], c1, c2, totmod,
                                        zout[j]);
          rat_rec(zout[j], totmod, res[j]);
        }

      }

      prime_i = (prime_i + 1) % BIG_UINT_PRIMES_SIZE;
    }

    return FAILED;
  }


} // namespace fflow
