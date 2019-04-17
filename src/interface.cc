#include <thread>
#include <future>
#include <fflow/interface.hh>
#include <fflow/mp_gcd.hh>

namespace fflow {

  // Algorithms

  void ExternalContext::alg_reserve(unsigned n)
  {
    alg_data.reserve(n);
  }

  void ExternalContext::alg_resize(unsigned n)
  {
    alg_data.resize(n);
  }

  void ExternalContext::init_subcontexts(unsigned n)
  {
    if (subctxt.size() < n) {
      subctxt.resize(n);
      for (unsigned i=0; i<n; ++i)
        if (subctxt[i].get() == nullptr)
          subctxt[i].reset(new ExternalContext());
    }
  }

  void ExternalContext::clear_subcontexts()
  {
    subctxt.clear();
  }

  void ExternalContext::ensure_subcontext_space()
  {
    unsigned idsize = alg_data.size();
    unsigned rfsize = regfun.size();
    for (auto & sc : subctxt) {
      if ((*sc).alg_data.size() < idsize)
        (*sc).alg_resize(idsize);
      if ((*sc).regfun.size() < rfsize)
        (*sc).regfun_resize(rfsize);
    }
    for (auto & sc : subctxt)
      (*sc).ww.ensure_size(ww.size());
  }

  void ExternalContext::alg_share_with_subctxt(unsigned id, unsigned i)
  {
    if ((*subctxt[i]).alg_data[id].alg.get() == nullptr)
      (*subctxt[i]).alg_data[id].alg = algid(id).clone(*subctxt[i]);
  }

  void ExternalContext::alg_share_with_subctxt(unsigned id,
                                               ExternalContext & sctxt)
  {
    if (sctxt.alg_data[id].alg.get() == nullptr)
      sctxt.alg_data[id].alg = algid(id).clone(sctxt);
  }

  unsigned ExternalContext::alg_get_freeid()
  {
    unsigned id = 0;

    if (alg_freeslots.empty()) {

      id = alg_data.size();
      if (id == 0)
        alg_reserve(DEFAULT_ALG_SLOTS);

      alg_data.resize(id+1);

    } else {

      id = alg_freeslots.back();
      alg_freeslots.pop_back();

    }

    return id;
  }


  void ExternalContext::alg_unregister(unsigned id)
  {
    alg_data[id] = AlgData();
    for (auto & t : subctxt)
      if (id < (*t).alg_data.size())
        (*t).alg_data[id] = AlgData();

    // insert id in free-slots, sorted from large to small values
    alg_freeslots.insert(std::upper_bound(alg_freeslots.begin(),
                                          alg_freeslots.end(),
                                          id,
                                          [](unsigned i, unsigned j)
                                          {
                                            return i > j;
                                          }), id);
  }


  void ExternalContext::alg_clear_all()
  {
    alg_data.clear();
    alg_freeslots.clear();
  }


  unsigned ExternalContext::dummy_alg(unsigned nparsin, unsigned nparsout)
  {
    std::unique_ptr<DummyAlgorithm> algptr(new DummyAlgorithm());
    DummyAlgorithm & alg = *algptr;
    alg.nparsin = nparsin;
    alg.nparsout = nparsout;
    unsigned id = alg_get_freeid();
    alg_data[id].alg = std::move(algptr);
    return id;
  }


  void ExternalContext::set_default_shift(unsigned id)
  {
    Mod mod(BIG_UINT_PRIMES[BIG_UINT_PRIMES_SIZE-1]);
    const unsigned nvars = algid(id).nparsin;
    auto & shift = algrecid(id).shift;
    if (nvars)
      shift.reset(new UInt[nvars]);
    for (unsigned i=0; i<nvars; ++i)
      shift[i] = sample_uint(OFFSET_4, i, mod);
  }


  void ExternalContext::alg_make_reconstructible(unsigned id)
  {
    if (alg_data[id].alg_rec.get() == nullptr) {
      alg_data[id].alg_rec.reset(new AlgRecData());
      set_default_shift(id);
    }
  }


  Ret ExternalContext::alg_learn(unsigned id)
  {
    const LearningOptions & opt = alglearnid(id);
    return algid(id).learn_all(Mod(BIG_UINT_PRIMES[opt.prime_no]),
                               opt.n_singular);
  }


  Ret ExternalContext::alg_degrees(unsigned id)
  {
    alg_make_reconstructible(id);

    Algorithm & a = algid(id);
    AlgRecData & arec = algrecid(id);
    AlgDegs & degs = arec.degs;
    const unsigned nparsout = a.nparsout;

    degs.numdeg.reset(new unsigned[nparsout]);
    degs.dendeg.reset(new unsigned[nparsout]);

    return algorithm_get_degrees(a, arec.shift.get(), arec.opt,
                                 degs.numdeg.get(), degs.dendeg.get());
  }


  Ret ExternalContext::alg_var_degrees(unsigned id, unsigned var)
  {
    Algorithm & a = algid(id);
    AlgRecData & arec = algrecid(id);
    AlgDegs & degs = arec.degs;

    RatFunVarDegrees * vdegs = degs.vdegs.get();
    return algorithm_get_var_degrees(a, var, arec.shift.get(), arec.opt,
                                     vdegs);
  }


  Ret ExternalContext::alg_all_var_degrees(unsigned id)
  {
    alg_make_reconstructible(id);

    Algorithm & a = algid(id);
    const unsigned nparsin = a.nparsin;
    const unsigned nparsout = a.nparsout;
    AlgDegs & degs = algrecid(id).degs;

    degs.vdegs.reset(new RatFunVarDegrees[nparsout]);
    for (unsigned i=0; i<nparsout; ++i)
      degs.vdegs[i].resize(a.nparsin);

    for (unsigned var=0; var<nparsin; ++var) {
      Ret ret = alg_var_degrees(id, var);
      if (ret != SUCCESS)
        return ret;
    }
    return SUCCESS;
  }


  Ret ExternalContext::alg_all_degrees(unsigned id)
  {
    Ret ret = alg_degrees(id);
    if (ret != SUCCESS)
      return ret;
    ret = alg_all_var_degrees(id);
    return ret;
  }


  struct ExtAlgParallelTotDegs {

    Ret operator() ()
    {
      par_ctxt->alg_share_with_subctxt(id, *sub_ctxt);

      Algorithm & a = sub_ctxt->algid(id);
      ExternalContext::AlgRecData & arec = par_ctxt->algrecid(id);
      ExternalContext::AlgDegs & degs = arec.degs;

      return algorithm_get_degrees(a, arec.shift.get(), arec.opt,
                                   degs.numdeg.get(), degs.dendeg.get());
    }

    ExternalContext *par_ctxt, *sub_ctxt;
    unsigned id;
  };

  struct ExtAlgParallelVarDegs {

    Ret operator() (unsigned var_beg, unsigned step)
    {
      par_ctxt->alg_share_with_subctxt(id, *sub_ctxt);

      const unsigned nparsin = sub_ctxt->algid(id).nparsin;
      Algorithm & a = sub_ctxt->algid(id);
      ExternalContext::AlgRecData & arec = par_ctxt->algrecid(id);
      ExternalContext::AlgDegs & degs = arec.degs;
      RatFunVarDegrees * vdegs = degs.vdegs.get();

      Ret ret = SUCCESS;
      for (unsigned v=var_beg; v<nparsin; v+=step) {
        ret = algorithm_get_var_degrees(a, v, arec.shift.get(), arec.opt,
                                        vdegs);
        if (ret != SUCCESS)
          return ret;
      }

      return SUCCESS;
    }

    ExternalContext *par_ctxt, *sub_ctxt;
    unsigned id;
  };

  Ret ExternalContext::alg_parallel_all_degrees(unsigned id,
                                                unsigned nthreads)
  {
    if (nthreads == 1)
      return alg_all_degrees(id);
    if (nthreads == 0)
      nthreads = std::thread::hardware_concurrency();

    alg_make_reconstructible(id);

    Algorithm & a = algid(id);
    AlgRecData & arec = algrecid(id);
    AlgDegs & degs = arec.degs;
    const unsigned nparsin = a.nparsin;
    const unsigned nparsout = a.nparsout;

    if (nthreads > nparsin+1)
      nthreads = nparsin+1;

    init_subcontexts(nthreads);
    ensure_subcontext_space();

    std::future<Ret> deg_ret;
    auto td_eval = ExtAlgParallelTotDegs{this, subctxt[0].get(), id};
    {
      degs.numdeg.reset(new unsigned[nparsout]);
      degs.dendeg.reset(new unsigned[nparsout]);
      deg_ret = std::async(std::launch::async, td_eval);
    }

    std::vector<std::future<Ret>> vdeg_ret(nthreads-1);
    std::vector<ExtAlgParallelVarDegs> t_eval(nthreads-1);
    {
      degs.vdegs.reset(new RatFunVarDegrees[nparsout]);
      for (unsigned i=0; i<nparsout; ++i)
        degs.vdegs[i].resize(a.nparsin);
      for (unsigned i=0; i<nthreads-1; ++i) {
        t_eval[i] = ExtAlgParallelVarDegs{this, subctxt[i+1].get(), id};
        vdeg_ret[i] = std::async(std::launch::async, t_eval[i], i, nthreads-1);
      }
    }

    Ret reti = deg_ret.get();
    if (reti != SUCCESS)
      return reti;
    for (unsigned i=0; i<nthreads-1; ++i) {
      reti = vdeg_ret[i].get();
      if (reti != SUCCESS)
        return reti;
    }

    return SUCCESS;
  }


  Ret ExternalContext::alg_dump_degree_info(unsigned id, const char * filename)
  {
    Algorithm & a = algid(id);
    AlgRecData & arec = algrecid(id);
    return algorithm_dump_degree_info(filename, a.nparsin, a.nparsout,
                                      arec.degs.numdeg.get(),
                                      arec.degs.dendeg.get(),
                                      arec.degs.vdegs.get());
  }


  Ret ExternalContext::alg_load_degree_info(unsigned id, const char * filename)
  {
    alg_make_reconstructible(id);

    Algorithm & a = algid(id);
    AlgRecData & arec = algrecid(id);
    AlgDegs & degs = arec.degs;
    const unsigned nparsout = a.nparsout;

    degs.numdeg.reset(new unsigned[nparsout]);
    degs.dendeg.reset(new unsigned[nparsout]);
    degs.vdegs.reset(new RatFunVarDegrees[nparsout]);

    return algorithm_load_degree_info(filename, a.nparsin, a.nparsout,
                                      arec.degs.numdeg.get(),
                                      arec.degs.dendeg.get(),
                                      arec.degs.vdegs.get());
  }


  void ExternalContext::init_cache_(CacheUPtr & cache,
                                    unsigned nparsin, unsigned nparsout,
                                    unsigned nsamples)
  {
    cache.reset(new UIntCache());
    (*cache).init(nparsin+1, nparsout);
    if (nsamples)
      (*cache).reserve_more(nsamples);
  }


  void ExternalContext::alg_init_cache(unsigned id)
  {
    init_cache_(algrecid(id).cache, algid(id).nparsin, algid(id).nparsout, 0);
  }


  void ExternalContext::alg_gen_sample_points(unsigned id,
                                              SamplePointsVector & samples)
  {
    Algorithm & a = algid(id);
    AlgRecData & arec = algrecid(id);
    AlgDegs & degs = arec.degs;

    GenerateSamplePoints sfun(a.nparsin);
    sfun.set_complement(arec.cache.get());

    algorithm_generate_sample_points(sfun, a.nparsin, a.nparsout,
                                     arec.shift.get(), arec.opt,
                                     degs.vdegs.get(),
                                     degs.numdeg.get(), degs.dendeg.get());
    sfun.append_to_vector(samples);
    sort_by_mod(samples.data(), samples.data()+samples.size(), a.nparsin);
  }


  Ret ExternalContext::alg_save_sample_points(unsigned id,
                                              const char * filename,
                                              const std::string cfiles[],
                                              unsigned n_cfiles)
  {
    Algorithm & a = algid(id);
    AlgRecData & arec = algrecid(id);
    AlgDegs & degs = arec.degs;
    SamplePointsVector samples;

    UIntCache cache;
    cache.init(a.nparsin, 1);
    UIntCache * cacheptr = nullptr;
    if (n_cfiles) {
      for (unsigned j=0; j<n_cfiles; ++j) {
        Ret ret = load_evaluations_for_var(cfiles[j].c_str(),
                                           a.nparsin, a.nparsout, cache, 0);
        if (ret == FAILED)
          return FAILED;
      }
      cacheptr = &cache;
    }

    GenerateSamplePoints sfun(a.nparsin);
    sfun.set_complement(cacheptr);
    algorithm_generate_sample_points(sfun, a.nparsin, a.nparsout,
                                     arec.shift.get(), arec.opt,
                                     degs.vdegs.get(),
                                     degs.numdeg.get(), degs.dendeg.get());
    sfun.append_to_vector(samples);
    sort_by_mod(samples.data(), samples.data()+samples.size(), a.nparsin);

    if (dump_samples(filename, a.nparsin, samples) == FAILED)
      return FAILED;

    return SUCCESS;
  }


  void ExternalContext::alg_sample(unsigned id)
  {
    SamplePointsVector samples;
    alg_gen_sample_points(id, samples);

    Algorithm & a = algid(id);
    AlgRecData & arec = algrecid(id);

    if (arec.cache.get() == nullptr)
      alg_init_cache(id);
    UIntCache & cache = *arec.cache;

    evaluate_and_cache_samples(a,
                               samples.data(), samples.data() + samples.size(),
                               cache);
  }


  Ret ExternalContext::alg_save_evaluations(unsigned id,
                                            const char * samples_file,
                                            const char * output_prefix,
                                            std::size_t tot_batches,
                                            std::size_t start_batch,
                                            std::size_t n_batches)
  {
    Algorithm & a = algid(id);
    SamplePointsVector samples;
    UInt tot_samples = samples_file_size(samples_file, a.nparsin);
    if (tot_samples == FAILED)
      return FAILED;

    std::size_t min_bsize = tot_samples / tot_batches;
    std::size_t rem_samples = tot_samples % tot_batches;

    for (std::size_t b = start_batch; b < start_batch + n_batches; ++b) {

      std::size_t bstart = 0;
      if (b < rem_samples)
        bstart = b*(min_bsize+1);
      else
        bstart = rem_samples*(min_bsize+1) + (b-rem_samples)*min_bsize;
      std::size_t bsize = min_bsize + (b < rem_samples);

      Ret ret = load_samples(samples_file, a.nparsin, samples, bstart, bsize);
      if (ret == FAILED)
        return FAILED;

      UIntCache cache;
      cache.init(a.nparsin+1, a.nparsout);
      cache.reserve_more(bsize);
      evaluate_and_cache_samples(a,
                                 samples.data(),
                                 samples.data() + samples.size(),
                                 cache);

      samples.clear();

      std::string output = format("{}.{}.ffloweval", output_prefix, b);
      dump_evaluations(output.c_str(), a.nparsin, a.nparsout, cache);
    }

    return SUCCESS;
  }


  struct ExtAlgParallelEvaluate {

    void operator() (CacheUPtr & tcache)
    {
      mctxt->alg_share_with_subctxt(id, *sctxt);
      Algorithm & alg = sctxt->algid(id);
      evaluate_and_cache_samples(alg, start, end, *tcache);
    }

    ExternalContext *mctxt, *sctxt;
    const std::unique_ptr<UInt[]> *start, *end;
    unsigned id;
  };

  void ExternalContext::alg_parallel_sample(unsigned id, unsigned nthreads)
  {
    if (nthreads == 1) {
      alg_sample(id);
      return;
    }

    if (nthreads == 0)
      nthreads = std::thread::hardware_concurrency();

    Algorithm & a = algid(id);
    AlgRecData & arec = algrecid(id);

    SamplePointsVector samples;
    alg_gen_sample_points(id, samples);

    const unsigned tot_samples = samples.size();
    std::unique_ptr<UInt[]> * start = samples.data();

    nthreads = std::min(nthreads, tot_samples);

    std::vector<CacheUPtr> caches(nthreads);
    std::vector<std::future<void>> ret(nthreads);
    std::vector<ExtAlgParallelEvaluate> t_eval(nthreads);

    init_subcontexts(nthreads);
    ensure_subcontext_space();

    for (unsigned i=0; i<nthreads; ++i) {
      const unsigned n = tot_samples/nthreads + (i < (tot_samples % nthreads));
      ExternalContext::init_cache_(caches[i], a.nparsin, a.nparsout, n);
    }

    for (unsigned i=0; i<nthreads; ++i) {
      const unsigned n = tot_samples/nthreads + (i < (tot_samples % nthreads));
      t_eval[i] = ExtAlgParallelEvaluate{this,
                                         subctxt[i].get(),
                                         start, start+n, id};
      ret[i] = std::async(std::launch::async, t_eval[i], std::ref(caches[i]));
      start += n;
    }

    for (unsigned i=0; i<nthreads; ++i)
      ret[i].get();

    if (arec.cache.get() == nullptr)
      alg_init_cache(id);
    UIntCache  & cache = *arec.cache;

    for (unsigned i=0; i<nthreads; ++i)
      merge_function_caches(caches.data(), nthreads, cache);
  }


  Ret ExternalContext::alg_reconstruct_(ExternalContext & parent,
                                        unsigned id,
                                        unsigned var_beg, unsigned step,
                                        MPReconstructedRatFun res[])
  {
    Algorithm & a = algid(id);
    AlgRecData & arec = parent.algrecid(id);
    AlgDegs & degs = arec.degs;
    UIntCache  & cache = *arec.cache;

    unsigned nparsin = a.nparsin, nparsout = a.nparsout;

    for (unsigned i=var_beg; i<nparsout; i+=step) {
      Ret ret = algorithm_reconstruct(cache, nparsin, nparsout, i,
                                      arec.shift.get(), arec.opt,
                                      degs.vdegs.get(),
                                      degs.numdeg.get(), degs.dendeg.get(),
                                      res[i]);
      if (ret != SUCCESS)
        return ret;
    }

    return SUCCESS;
  }

  Ret ExternalContext::alg_reconstruct(unsigned id,
                                       MPReconstructedRatFun res[])
  {
    return alg_reconstruct_(*this, id, 0, 1, res);
  }

  struct ExtAlgParallelReconstruct {

    Ret operator() (unsigned var_beg, unsigned step,
                    MPReconstructedRatFun res[])
    {
      par_ctxt->alg_share_with_subctxt(id, *this_ctxt);
      return this_ctxt->alg_reconstruct_(*par_ctxt, id, var_beg, step, res);
    }

    ExternalContext *par_ctxt, *this_ctxt;
    unsigned id;
  };

  Ret ExternalContext::alg_parallel_reconstruct(unsigned id,
                                                MPReconstructedRatFun res[],
                                                unsigned nthreads)
  {
    if (nthreads == 1)
      return alg_reconstruct(id, res);

    if (nthreads == 0)
      nthreads = std::thread::hardware_concurrency();

    unsigned nparsout = algid(id).nparsout;
    nthreads = std::min(nthreads, nparsout);

    init_subcontexts(nthreads);
    ensure_subcontext_space();

    std::vector<ExtAlgParallelReconstruct> t_eval(nthreads);
    std::vector<std::future<Ret>> ret(nthreads);
    for (unsigned i=0; i<nthreads; ++i)
      t_eval[i] = ExtAlgParallelReconstruct{this, subctxt[i].get(), id};

    for (unsigned i=0; i<nthreads; ++i)
      ret[i] = std::async(std::launch::async, t_eval[i], i, nthreads, res);

    for (unsigned i=0; i<nthreads; ++i) {
      Ret reti = ret[i].get();
      if (reti != SUCCESS)
          return reti;
    }

    return SUCCESS;
  }


  Ret ExternalContext::alg_reconstruct_univariate(unsigned id,
                                                  MPReconstructedRatFun res[])
  {
    Algorithm & a = algid(id);
    if (a.nparsin != 1)
      return FAILED;
    alg_make_reconstructible(id);
    AlgRecData & arec = algrecid(id);
    return algorithm_reconstruct_univariate(a, arec.shift.get(), arec.opt, res);
  }


  Ret ExternalContext::alg_reconstruct_numeric(unsigned id, MPRational res[])
  {
    Algorithm & a = algid(id);
    if (a.nparsin != 0)
      return FAILED;
    alg_make_reconstructible(id);
    AlgRecData & arec = algrecid(id);
    return algorithm_reconstruct_numeric(a, arec.opt, res);
  }


  UInt ExternalContext::alg_independent_of_var(unsigned id, unsigned var)
  {
    Algorithm & a = algid(id);
    unsigned nparsin = a.nparsin, nparsout = a.nparsout;
    std::unique_ptr<UInt> xdata(new UInt[nparsin+2*nparsout]);

    UInt * xin = xdata.get();
    UInt * xout1 = xdata.get() + nparsin;
    UInt * xout2 = xdata.get() + nparsin + nparsout;
    unsigned n_singular = alg_data[id].alg_learn.n_singular;

    Mod mod(BIG_UINT_PRIMES[0]);

    const unsigned NEEDED_CHECKS = 2;
    unsigned iii = 0;
    unsigned checks = 0;
    unsigned fails = 0;

    while (fails < n_singular) {

      for (unsigned i=0; i<nparsin; ++i)
        if (i != var)
          xin[i] = sample_uint(OFFSET_1, iii++, mod);
      xin[var] = sample_uint(OFFSET_2, 2*checks, mod);

      Ret ret = a.evaluate(xin, mod, xout1);
      if (ret == FAILED) {
        ++fails;
        continue;
      };

      xin[var] = sample_uint(OFFSET_2, 2*checks+1, mod);
      ret = a.evaluate(xin, mod, xout2);
      if (ret == FAILED) {
        ++fails;
        continue;
      };

      if (!std::equal(xout1, xout1+nparsout, xout2))
        return UInt(false);

      ++checks;
      if (checks == NEEDED_CHECKS)
        return UInt(true);
    }

    return FAILED;
  }

  Ret ExternalContext::alg_reconstruct_from_evals(unsigned id,
                                                  const std::string files[],
                                                  unsigned nfiles,
                                                  const char * degs_file,
                                                  unsigned idx,
                                                  MPReconstructedRatFun & recf)
  {
    Algorithm & a = algid(id);
    alg_make_reconstructible(id);
    AlgRecData & arec = algrecid(id);
    return algorithm_reconstruct_from_evals(files, nfiles, degs_file,
                                            a.nparsin, a.nparsout, idx,
                                            arec.shift.get(), arec.opt, recf);
  }


  // Registered functions

  unsigned ExternalContext::regfun_get_freeid()
  {
    unsigned id = 0;

    if (regfun_freeslots.empty()) {

      id = regfun.size();
      regfun_resize(id+1);

    } else {

      id = regfun_freeslots.back();
      regfun_freeslots.pop_back();

    }

    return id;
  }

  void ExternalContext::regfun_unregister(unsigned id)
  {
    regfun[id] = HornerRatFunPtr();
    regfun_map[id].clear();
    regfun_mod[id] = 0;
    regfun_res[id] = 0;

    // insert id in free-slots, sorted from large to small values
    regfun_freeslots.insert(std::upper_bound(regfun_freeslots.begin(),
                                             regfun_freeslots.end(),
                                             id,
                                             [](unsigned i, unsigned j)
                                             {
                                               return i > j;
                                             }), id);
  }

  void ExternalContext::regfun_clear_all()
  {
    regfun.clear();
    regfun_res.clear();
    regfun_map.clear();
    regfun_mod.clear();
    regfun_freeslots.clear();
  }

  void ExternalContext::regfun_resize(unsigned n)
  {
    regfun.resize(n);
    regfun_res.resize(n);
    regfun_map.resize(n);
    regfun_mod.resize(n);
  }

  UInt ExternalContext::regfun_eval(unsigned id, unsigned nvars,
                                    const UInt x[], const UInt xp_shoup[],
                                    Mod mod)
  {
    return regfun[id].eval(nvars, ww.get(), x, xp_shoup, mod);
  }

  UInt ExternalContext::regfun_eval_and_store(unsigned id, unsigned nvars,
                                              const UInt x[],
                                              const UInt xp_shoup[],
                                              Mod mod)
  {
    UInt res = regfun_eval(id, nvars, x, xp_shoup, mod);
    regfun_res[id] = res;
    return res;
  }

  void ExternalContext::regfun_share_with_subctxt(unsigned id,
                                                  unsigned nparsin,
                                                  unsigned i)
  {
    regfun_share_with_subctxt(id, nparsin, *subctxt[i]);
  }

  void ExternalContext::regfun_share_with_subctxt(unsigned id,
                                                  unsigned nparsin,
                                                  ExternalContext & sctxt)
  {
    if (regfun[id].num() != nullptr && sctxt.regfun[id].num() == nullptr)
      horner_mapped_ratfun_clone(regfun_map[id], regfun[id], nparsin,
                                 sctxt.regfun_map[id], sctxt.regfun[id]);
  }


  // External algorithms

  Ret DummyAlgorithm::evaluate(const UInt[], Mod, UInt[])
  {
    return FAILED;
  }

  std::unique_ptr<Algorithm> DummyAlgorithm::clone(ExternalContext &) const
  {
    return nullptr;
  }


  void ExternalDenseSystem::reset_mod(Mod mod)
  {
    this_mod = mod.n();
    if (xp.get() == nullptr)
      xp.reset(new UInt[nparsin]);
    const unsigned n_cols = c.columns();
    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();
    for (unsigned i=0; i<nindepeqs; ++i)
      for (unsigned j=0; j<n_cols; ++j)
        horner_ratfunmap_mprational(cmap(ieq[i],j), mod);
  }

  Ret ExternalDenseSystem::fill_matrix(std::size_t n_rows,
                                       const std::size_t rows[],
                                       const UInt xi[], Mod mod,
                                       MatrixView & m)
  {
    if (mod.n() != this_mod)
      reset_mod(mod);
    precomp_array_mul_shoup(xi, nparsin, mod, xp.get());

    UInt * ww = ctxt->ww.get();
    const unsigned n_cols = c.columns();

    for (unsigned i=0; i<n_rows; ++i) {

      const HornerRatFunPtr * f = c.row(rows[i]);
      const HornerRatFunPtr * fend = f + n_cols;
      UInt * r = m.row(i);
      const DenseLinearSolver::flag_t * info = xinfo();

      for (; f<fend; ++f, ++r, ++info)
        if ((*info) & LSVar::IS_NON_ZERO) {
          UInt res = (*f).eval(nparsin, ww, xi, xp.get(), mod);
          if (res == FAILED)
            return FAILED;
          *r = res;
        }
    }

    return SUCCESS;
  }


  std::unique_ptr<Algorithm>
  ExternalDenseSystem::clone(ExternalContext & ectxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalDenseSystem());
    auto & newalg = static_cast<ExternalDenseSystem &>(*ptr);

    DenseLinearSolver::copy_into(newalg);

    newalg.c.resize(c.rows(), c.columns());
    newalg.cmap.resize(cmap.rows(), cmap.columns());

    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();
    for (unsigned i=0; i<nindepeqs; ++i)
      for (unsigned j=0; j<c.columns(); ++j)
        horner_mapped_ratfun_clone(cmap(ieq[i],j), c(ieq[i],j), nparsin,
                                   newalg.cmap(ieq[i],j), newalg.c(ieq[i],j));

    newalg.ctxt = &ectxt;

    return ptr;
  }



  void ExternalSparseSystem::reset_mod(Mod mod)
  {
    this_mod = mod.n();
    if (xp.get() == nullptr)
      xp.reset(new UInt[nparsin]);
    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();
    for (unsigned ind=0; ind<nindepeqs; ++ind) {
      std::size_t i = ieq[ind];
      std::size_t row_size = rinfo[i].size;
      for (unsigned j=0; j<row_size; ++j)
        horner_ratfunmap_mprational(cmap[i][j], mod);
    }
  }

  Ret ExternalSparseSystem::fill_matrix(std::size_t n_rows,
                                        const std::size_t rows[],
                                        const UInt xi[], Mod mod,
                                        SparseMatrix & m)
  {
    if (mod.n() != this_mod)
      reset_mod(mod);
    precomp_array_mul_shoup(xi, nparsin, mod, xp.get());

    UInt * ww = ctxt->ww.get();

    for (unsigned i=0; i<n_rows; ++i) {

      SparseMatrixRow & r = m.row(i);

      const std::size_t row_size = rinfo[rows[i]].size;
      const unsigned * cols = rinfo[rows[i]].cols.get();
      const HornerRatFunPtr * f = c[rows[i]].get();
      const HornerRatFunPtr * fend = f + row_size;

      r.resize(row_size);
      unsigned j=0;

      for (; f<fend; ++f, ++j) {
        UInt res = (*f).eval(nparsin, ww, xi, xp.get(), mod);
        unsigned col = cols[j];
        if (res == FAILED || res == 0)
          return FAILED;
        r.el(j).col = col;
        r.el(j).val = res;
      }
      r.el(j).col = SparseMatrixRow::END;
    }

    return SUCCESS;
  }


  std::unique_ptr<Algorithm>
  ExternalSparseSystem::clone(ExternalContext & ectxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalSparseSystem());
    auto & newalg = static_cast<ExternalSparseSystem &>(*ptr);

    SparseLinearSolver::copy_into(newalg);

    std::size_t n_rows = rinfo.size();

    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();

    newalg.rinfo.resize(rinfo.size());
    for (unsigned ind=0; ind<nindepeqs; ++ind) {
      std::size_t i = ieq[ind];
      std::size_t size =  newalg.rinfo[i].size = rinfo[i].size;
      newalg.rinfo[i].cols.reset(new unsigned[size]);
      for (unsigned j=0; j<size; ++j)
        newalg.rinfo[i].cols[j] = rinfo[i].cols[j];
    }

    newalg.c.resize(n_rows);
    newalg.cmap.resize(n_rows);

    for (unsigned ind=0; ind<nindepeqs; ++ind) {
      std::size_t i = ieq[ind];
      std::size_t row_size = rinfo[i].size;
      newalg.c[i].reset(new HornerRatFunPtr[row_size]);
      newalg.cmap[i].reset(new MPHornerRatFunMap[row_size]);
      for (unsigned j=0; j<row_size; ++j)
        horner_mapped_ratfun_clone(cmap[i][j], c[i][j], nparsin,
                                   newalg.cmap[i][j], newalg.c[i][j]);
    }

    newalg.ctxt = &ectxt;

    return ptr;
  }


  void ExternalSparseSystem::delete_unneeded_eqs()
  {
    std::size_t n_rows = rinfo.size();
    std::unique_ptr<bool[]> needed(new bool[n_rows]());

    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();

    for (unsigned i=0; i<nindepeqs; ++i)
      needed[ieq[i]] = true;

    for (unsigned i=0; i<n_rows; ++i)
      if (!needed[i]) {
        rinfo[i] = RowInfo();
        c[i].reset(nullptr);
        cmap[i].reset(nullptr);
      }
  }


  Ret ExternalDenseSystemN::fill_matrix(std::size_t n_rows,
                                        const std::size_t rows[],
                                        const UInt[], Mod mod,
                                        MatrixView & m)
  {
    const MPInt mpmod(mod.n());
    MPInt mpres;
    const unsigned n_cols = c.columns();

    for (unsigned i=0; i<n_rows; ++i) {

      const MPRational * f = c.row(rows[i]);
      const MPRational * fend = f + n_cols;
      UInt * r = m.row(i);
      const DenseLinearSolver::flag_t * info = xinfo();

      for (; f<fend; ++f, ++r, ++info)
        if ((*info) & LSVar::IS_NON_ZERO) {
          rat_mod((*f), mpmod, mpres);
          *r = mpres.to_uint();
        }
    }

    return SUCCESS;
  }

  std::unique_ptr<Algorithm>
  ExternalDenseSystemN::clone(ExternalContext & ectxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalDenseSystemN());
    auto & newalg = static_cast<ExternalDenseSystemN &>(*ptr);

    DenseLinearSolver::copy_into(newalg);

    newalg.c.resize(c.rows(), c.columns());

    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();
    for (unsigned i=0; i<nindepeqs; ++i)
      for (unsigned j=0; j<c.columns(); ++j)
        newalg.c(ieq[i], j) = c(ieq[i], j);

    newalg.ctxt = &ectxt;

    return ptr;
  }


  Ret ExternalSparseSystemN::fill_matrix(std::size_t n_rows,
                                         const std::size_t rows[],
                                         const UInt[], Mod mod,
                                         SparseMatrix & m)
  {
    const MPInt mpmod(mod.n());
    MPInt mpres;

    for (unsigned i=0; i<n_rows; ++i) {

      SparseMatrixRow & r = m.row(i);

      const std::size_t row_size = rinfo[rows[i]].size;
      const unsigned * cols = rinfo[rows[i]].cols.get();
      const MPRational * f = c[rows[i]].get();
      const MPRational * fend = f + row_size;

      r.resize(row_size);
      unsigned j=0;

      for (; f<fend; ++f, ++j) {
        rat_mod(*f, mpmod, mpres);
        unsigned col = cols[j];
        r.el(j).col = col;
        r.el(j).val = mpres.to_uint();
      }
      r.el(j).col = SparseMatrixRow::END;
    }

    return SUCCESS;
  }


  std::unique_ptr<Algorithm>
  ExternalSparseSystemN::clone(ExternalContext & ectxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalSparseSystemN());
    auto & newalg = static_cast<ExternalSparseSystemN &>(*ptr);

    SparseLinearSolver::copy_into(newalg);

    std::size_t n_rows = rinfo.size();

    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();

    newalg.rinfo.resize(rinfo.size());
    for (unsigned ind=0; ind<nindepeqs; ++ind) {
      std::size_t i = ieq[ind];
      std::size_t size =  newalg.rinfo[i].size = rinfo[i].size;
      newalg.rinfo[i].cols.reset(new unsigned[size]);
      for (unsigned j=0; j<size; ++j)
        newalg.rinfo[i].cols[j] = rinfo[i].cols[j];
    }

    newalg.c.resize(n_rows);

    for (unsigned ind=0; ind<nindepeqs; ++ind) {
      std::size_t i = ieq[ind];
      std::size_t row_size = rinfo[i].size;
      newalg.c[i].reset(new MPRational[row_size]);
      for (unsigned j=0; j<row_size; ++j)
        newalg.c[i][j] = c[i][j];
    }

    newalg.ctxt = &ectxt;

    return ptr;
  }

  void ExternalSparseSystemN::delete_unneeded_eqs()
  {
    std::size_t n_rows = rinfo.size();
    std::unique_ptr<bool[]> needed(new bool[n_rows]());

    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();

    for (unsigned i=0; i<nindepeqs; ++i)
      needed[ieq[i]] = true;

    for (unsigned i=0; i<n_rows; ++i)
      if (!needed[i]) {
        rinfo[i] = RowInfo();
        c[i].reset(nullptr);
      }
  }


  Ret ExternalAlgCompose::init(ExternalContext * ectxt,
                               unsigned id_in, unsigned id_out)
  {
    ctxt = ectxt;

    Algorithm & algin = ctxt->algid(id_in);
    Algorithm & algout = ctxt->algid(id_out);

    if (algin.nparsout != algout.nparsin)
      return FAILED;

    nparsin = algin.nparsin;
    nparsout = algout.nparsout;
    alg_idin = id_in;
    alg_idout = id_out;
    xoutin.reset(new UInt[algin.nparsout]);

    return SUCCESS;
  }

  UInt ExternalAlgCompose::min_learn_times()
  {
    return ctxt->algid(alg_idout).min_learn_times();
  }

  Ret ExternalAlgCompose::learn(const UInt xin[], Mod mod)
  {
    Ret ret = ctxt->algid(alg_idin).evaluate(xin, mod, xoutin.get());
    if (ret == FAILED)
      return FAILED;
    Algorithm & algout = ctxt->algid(alg_idout);
    ret = algout.learn(xoutin.get(), mod);
    nparsout = algout.nparsout;
    return ret;
  }

  Ret ExternalAlgCompose::evaluate(const UInt xin[], Mod mod, UInt xout[])
  {
    Ret ret = ctxt->algid(alg_idin).evaluate(xin, mod, xoutin.get());
    if (ret == FAILED)
      return FAILED;
    return ctxt->algid(alg_idout).evaluate(xoutin.get(), mod, xout);
  }

  std::unique_ptr<Algorithm>
  ExternalAlgCompose::clone(ExternalContext & sctxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalAlgCompose());
    auto & newalg = static_cast<ExternalAlgCompose &>(*ptr);

    ctxt->alg_share_with_subctxt(alg_idin, sctxt);
    ctxt->alg_share_with_subctxt(alg_idout, sctxt);

    newalg.init(&sctxt, alg_idin, alg_idout);

    return ptr;
  }


  void ExternalAlgMul::new_mod(Mod mod)
  {
    this_mod = mod.n();
    if (xp.get() == nullptr)
      xp.reset(new UInt[nparsin]);
    horner_ratfunmap_mprational(cmap, mod);
  }

  Ret ExternalAlgMul::evaluate(const UInt xin[], Mod mod, UInt xout[])
  {
    if (mod.n() != this_mod)
      new_mod(mod);
    precomp_array_mul_shoup(xin, nparsin, mod, xp.get());

    UInt fval = cfun.eval(nparsin, ctxt->ww.get(), xin, xp.get(), mod);
    if (fval == FAILED)
      return FAILED;

    Ret aret = ctxt->algid(alg_id).evaluate(xin, mod, xout);
    if (aret != SUCCESS)
      return aret;

    UInt fvalp = precomp_mul_shoup(fval, mod);
    for (unsigned i=0; i<nparsout; ++i)
      xout[i] = mul_mod_shoup(xout[i], fval, fvalp, mod);

    return SUCCESS;
  }

  std::unique_ptr<Algorithm>
  ExternalAlgMul::clone(ExternalContext & sctxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalAlgMul());
    auto & newalg = static_cast<ExternalAlgMul &>(*ptr);

    ctxt->alg_share_with_subctxt(alg_id, sctxt);

    newalg.nparsin = nparsin;
    newalg.nparsout = nparsout;
    horner_mapped_ratfun_clone(cmap, cfun, nparsin, newalg.cmap, newalg.cfun);
    newalg.alg_id = alg_id;

    newalg.ctxt = &sctxt;

    return ptr;
  }


  Ret ExternalAlgBindVar::new_mod(Mod mod)
  {
    if (x.get() == nullptr)
      x.reset(new UInt[nparsin+1]);
    MPInt res;
    bool ret = rat_mod(rval, MPInt(mod.n()), res);
    if (!ret)
      return FAILED;
    val = res.to_uint();
    this_mod = mod.n();
    return SUCCESS;
  }

  Ret ExternalAlgBindVar::evaluate(const UInt xin[], Mod mod, UInt xout[])
  {
    if (mod.n() != this_mod) {
      UInt ret = new_mod(mod);
      if (ret == FAILED)
        return FAILED;
    }

    unsigned pos = 0;
    for (unsigned i=0; i<nparsin+1; ++i)
      if (i != var)
        x[i] = xin[pos++];
    x[var] = val;

    return ctxt->algid(alg_id).evaluate(x.get(), mod, xout);
  }

  std::unique_ptr<Algorithm>
  ExternalAlgBindVar::clone(ExternalContext & sctxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalAlgBindVar());
    auto & newalg = static_cast<ExternalAlgBindVar &>(*ptr);

    ctxt->alg_share_with_subctxt(alg_id, sctxt);

    newalg.nparsin = nparsin;
    newalg.nparsout = nparsout;

    newalg.x.reset(new UInt[newalg.nparsin+1]);
    newalg.val = val;
    newalg.this_mod = this_mod;
    newalg.rval = rval;
    newalg.var = var;
    newalg.alg_id = alg_id;

    newalg.ctxt = &sctxt;

    return ptr;
  }


  void ExternalAlgRatFunEval::new_mod(Mod mod)
  {
    if (xp.get() == nullptr)
      xp.reset(new UInt[nparsin]);
    this_mod = mod.n();
    for (unsigned i=0; i<nparsout; ++i)
      horner_ratfunmap_mprational(fmap[i], mod);
  }

  Ret ExternalAlgRatFunEval::evaluate(const UInt xin[], Mod mod, UInt xout[])
  {
    if (mod.n() != this_mod)
      new_mod(mod);

    precomp_array_mul_shoup(xin, nparsin, mod, xp.get());

    UInt * ww = ctxt->ww.get();
    for (unsigned i=0; i<nparsout; ++i) {
      UInt res = f[i].eval(nparsin, ww, xin, xp.get(), mod);
      if (res == FAILED)
        return FAILED;
      xout[i] = res;
    }

    return SUCCESS;
  }

  std::unique_ptr<Algorithm>
  ExternalAlgRatFunEval::clone(ExternalContext & ectxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalAlgRatFunEval());
    auto & newalg = static_cast<ExternalAlgRatFunEval &>(*ptr);

    newalg.nparsin = nparsin;
    newalg.nparsout = nparsout;

    newalg.f.reset(new HornerRatFunPtr[nparsout]);
    newalg.fmap.reset(new MPHornerRatFunMap[nparsout]);
    for (unsigned i=0; i<nparsout; ++i)
      horner_mapped_ratfun_clone(fmap[i], f[i], nparsin,
                                 newalg.fmap[i], newalg.f[i]);

    newalg.ctxt = &ectxt;

    return ptr;
  }


  Ret ExternalAlgSum::evaluate(const UInt xin[], Mod mod, UInt xout[])
  {
    if (partial_xout.get() == nullptr)
      partial_xout.reset(new UInt[nparsout]);

    std::fill(xout, xout + nparsout, 0);
    for (unsigned id : alg_ids) {
      Ret ret = ctxt->algid(id).evaluate(xin, mod, partial_xout.get());
      if (ret != SUCCESS)
        return ret;
      for (unsigned i=0; i<nparsout; ++i)
        xout[i] = add_mod(xout[i], partial_xout[i], mod);
    }

    return SUCCESS;
  }

  std::unique_ptr<Algorithm>
  ExternalAlgSum::clone(ExternalContext & ectxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalAlgSum());
    auto & newalg = static_cast<ExternalAlgSum &>(*ptr);

    for (unsigned id : alg_ids)
      ctxt->alg_share_with_subctxt(id, ectxt);

    newalg.nparsin = nparsin;
    newalg.nparsout = nparsout;
    newalg.alg_ids = alg_ids;
    newalg.ctxt = &ectxt;

    return ptr;
  }


  Ret ExternalAlgMatMul::evaluate(const UInt xin[], Mod mod, UInt xout[])
  {
    Ret ret = 0;

    ret = ctxt->algid(alg_id1).evaluate(xin, mod, out1.get());
    if (ret != SUCCESS)
      return ret;

    ret = ctxt->algid(alg_id2).evaluate(xin, mod, out2.get());
    if (ret != SUCCESS)
      return ret;

    std::fill(xout, xout + nparsout, 0);
    for (unsigned i=0; i<nrows1; ++i)
      for (unsigned k=0; k<ncols1; ++k)
        for (unsigned j=0; j<ncols2; ++j)
          xout[i*ncols2 + j] = add_mod(xout[i*ncols2 + j],
                                       mul_mod(out1[i*ncols1 + k],
                                               out2[k*ncols2 + j],
                                               mod),
                                       mod);

    return SUCCESS;
  }

  std::unique_ptr<Algorithm>
  ExternalAlgMatMul::clone(ExternalContext & ectxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalAlgMatMul());
    auto & newalg = static_cast<ExternalAlgMatMul &>(*ptr);

    newalg.out1.reset(new UInt[nrows1*ncols1]);
    newalg.out2.reset(new UInt[ncols1*ncols2]);

    ctxt->alg_share_with_subctxt(alg_id1, ectxt);
    ctxt->alg_share_with_subctxt(alg_id2, ectxt);

    newalg.nrows1 = nrows1;
    newalg.ncols1 = ncols1;
    newalg.ncols2 = ncols2;
    newalg.alg_id1 = alg_id1;
    newalg.alg_id2 = alg_id2;
    newalg.nparsin = nparsin;
    newalg.nparsout = nparsout;

    newalg.ctxt = &ectxt;

    return ptr;
  }


  void ExternalAlgSparse::init(ExternalContext * ectxt, unsigned id)
  {
    ctxt = ectxt;

    Algorithm & algin = ctxt->algid(id);
    out.reset(new UInt[algin.nparsout]);

    nparsin = algin.nparsin;
    nparsout = UNINITIALIZED;
    alg_id = id;
  }

  UInt ExternalAlgSparse::min_learn_times()
  {
    return 2;
  }

  Ret ExternalAlgSparse::learn(const UInt xin[], Mod mod)
  {
    Algorithm & alg = ctxt->algid(alg_id);
    unsigned nout = alg.nparsout;

    Ret ret = alg.evaluate(xin, mod, out.get());
    if (ret == FAILED)
      return FAILED;

    unsigned nnz = std::count(out.get(), out.get() + nout, UInt(0));

    if (nparsout == UNINITIALIZED) {
      nparsout = nnz;
      nonzero.reset(new unsigned[nnz]);
      unsigned idx=0;
      for (unsigned j=0; j<nout; ++j)
        if (out[j] != 0)
          nonzero[idx++] = j;
    } else if (nparsout == nnz) {
      unsigned idx=0;
      for (unsigned j=0; j<nout; ++j)
        if (out[j] != 0 && nonzero[idx++] != j) {
          nparsout = UNINITIALIZED;
          return FAILED;
        }
    } else {
      nparsout = UNINITIALIZED;
      return FAILED;
    }

    nparsout = nnz;
    return SUCCESS;
  }

  Ret ExternalAlgSparse::evaluate(const UInt xin[], Mod mod, UInt xout[])
  {
    Algorithm & alg = ctxt->algid(alg_id);
    Ret ret = alg.evaluate(xin, mod, out.get());
    if (ret == FAILED)
      return FAILED;
    for (unsigned j=0; j<nparsout; ++j)
      xout[j] = out[nonzero[j]];
    return SUCCESS;
  }

  std::unique_ptr<Algorithm>
  ExternalAlgSparse::clone(ExternalContext & sctxt) const
  {
    std::unique_ptr<Algorithm> ptr(new ExternalAlgSparse());
    auto & newalg = static_cast<ExternalAlgSparse &>(*ptr);

    ctxt->alg_share_with_subctxt(alg_id, sctxt);

    newalg.init(&sctxt, alg_id);

    newalg.nparsout = nparsout;
    newalg.nonzero.reset(new unsigned[nparsout]);
    std::copy(nonzero.get(), nonzero.get()+nparsout, newalg.nonzero.get());

    return ptr;
  }


  void ExternalTreeLevel::reset_ext_spinors_mod(ExternalSpinorsMap & sp,
                                                Mod mod)
  {
    unsigned n = sp.rows();
    for (unsigned i=0; i<n; ++i)
      for (unsigned j=0; j<4; ++j)
        horner_ratfunmap_mprational(sp(i,j), mod);
  }

  Ret ExternalTreeLevel::init_spinors_from_ext(unsigned nvars,
                                               UInt * __restrict ww,
                                               const UInt * __restrict x,
                                               const UInt * __restrict xp,
                                               Mod mod,
                                               ExternalSpinors & sp,
                                               Spinor4D spinor[])
  {
    unsigned n = sp.rows();
    UInt res;

    for (unsigned i=0; i<n; ++i) {

      res = spinor[i].l[0] = sp(i,0).eval(nvars, ww, x, xp, mod);
      if (res == FAILED)
        return FAILED;

      res = spinor[i].l[1] = sp(i,1).eval(nvars, ww, x, xp, mod);
      if (res == FAILED)
        return FAILED;

      res = spinor[i].lt[0] = sp(i,2).eval(nvars, ww, x, xp, mod);
      if (res == FAILED)
        return FAILED;

      res = spinor[i].lt[1] = sp(i,3).eval(nvars, ww, x, xp, mod);
      if (res == FAILED)
        return FAILED;

    }

    return SUCCESS;
  }


  void ExternalUnitarityCut::init(unsigned nlegs, const UInt types[],
                                  unsigned ncoeffs, unsigned nparams,
                                  unsigned nloopexpr, unsigned nloopvars,
                                  unsigned needed_coeffs[],
                                  unsigned needed_size,
                                  unsigned extra_eqs)
  {
    nparsin = nparams;
    UnitarityCutFit::init(ncoeffs, nloopvars, needed_coeffs, needed_size,
                          extra_eqs);
    UnitarityCutFit::init_kinematics(nlegs, types);

    ntau = nloopvars;
    nlvar = nloopexpr ? nloopexpr : nloopvars;

    spexpr.resize(nspinors(), 4);
    spexprmap.resize(nspinors(), 4);
    lvar.reset(new UInt[nlvar]);
    shoup.reset(new UInt[std::max<std::size_t>(nparsin, nlvar)]);
    lvar_shoup.reset(new UInt[nlvar]);
  }

  void ExternalUnitarityCut::reset_mod(Mod mod)
  {
    this_mod = mod.n();

    MPInt tmp, mmod(mod.n());
    for (unsigned i=0; i<nlflav; ++i) {
      rat_mod(lflav_coeff[i], mmod, tmp);
      lflav[i].coeff = tmp.to_uint();
    }

    const flag_t * xinf = xinfo();
    for (auto & cf : cmpmap)
      if ((*xinf++ && LSVar::IS_NON_ZERO))
        map_mp_data(cf, mod);
    for (auto & cf : integrmpmap)
      map_mp_data(cf, mod);
    for (auto & cf : lexprmpmap)
      map_mp_data(cf, mod);

    unsigned loops = kexpr.rows();
    for (unsigned i=0; i<loops; ++i)
      for (unsigned j=0; j<6; ++j)
        for (auto & kmap : kexprmpmap(i,j))
          horner_polymap_mprational(kmap, mod);

    if (weights.get()) {
      unsigned intsize = integr.size();
      for (unsigned i=0; i<intsize; ++i) {
        unsigned id = weights[i];
        if (ctxt->regfun_mod[id] != mod.n()) {
          horner_ratfunmap_mprational(ctxt->regfun_map[id], mod);
          ctxt->regfun_mod[id] = mod.n();
        }
      }
    }

    ExternalTreeLevel::reset_ext_spinors_mod(spexprmap, mod);
    horner_ratfunmap_mprational(prefmap, mod);
  }

  Ret ExternalUnitarityCut::new_xpoint(const UInt x[], Mod mod)
  {
    if (this_mod != mod.n())
      reset_mod(mod);

    UInt * ww = ctxt->ww.get();

    const flag_t * xinf = xinfo();
    precomp_array_mul_shoup(x, nparsin, mod, shoup.get());
    for (auto & cf : cmap)
      if ((*xinf++ && LSVar::IS_NON_ZERO))
        horner_ratfunmap_horner(cf, nparsin, ww, x, shoup.get(), mod);
    for (auto & cf : integrmap)
      horner_ratfunmap_horner(cf, nparsin, ww, x, shoup.get(), mod);
    for (auto & cf : lexprmap)
      horner_ratfunmap_horner(cf, nparsin, ww, x, shoup.get(), mod);

    evaluate_loop_coeffs(x, mod);

    if (weights.get()) {
      unsigned intsize = integr.size();
      for (unsigned i=0; i<intsize; ++i) {
        unsigned id = weights[i];
        UInt ret = ctxt->regfun_eval_and_store(id, nparsin, x, shoup.get(),
                                               mod);
        if (ret == FAILED)
          return FAILED;
      }
    }

    return UnitarityCutFit::new_xpoint(x, mod);
  }

  Ret ExternalUnitarityCut::init_spinors(const UInt x[], Mod mod,
                                         Spinor4D spinor[])
  {
    return ExternalTreeLevel::init_spinors_from_ext(nparsin,
                                                    ctxt->ww.get(),
                                                    x, shoup.get(),
                                                    mod, spexpr,
                                                    spinor);
  }

  void ExternalUnitarityCut::evaluate_loop_coeffs(const UInt x[], Mod mod)
  {
    UInt * ww = ctxt->ww.get();
    unsigned loops = kexpr.rows();
    for (unsigned i=0; i<loops; ++i)
      for (unsigned j=0; j<6; ++j)
        horner_ratfunmap_horner(kexprmap(i,j), nparsin, ww,
                                x, shoup.get(), mod);
  }

  UInt ExternalUnitarityCut::get_pref(const UInt x[], Mod mod)
  {
    return pref.eval(nparsin, ctxt->ww.get(), x, shoup.get(), mod);
  }

  Ret ExternalUnitarityCut::get_cut_momenta(const UInt[],
                                            const UInt tau_in[],
                                            Mod mod,
                                            Momentum6D k[])
  {
    UInt * ww = ctxt->ww.get();

    // get cut solutions
    precomp_array_mul_shoup(tau_in, ntau, mod, shoup.get());
    unsigned loops = kexpr.rows();
    UInt res = 0;
    for (unsigned i=0; i<loops; ++i)
      for (unsigned j=0; j<6; ++j) {
        res = kexpr(i,j).eval(ntau, ww, tau_in, shoup.get(), mod);
        if (res == FAILED)
          return FAILED;
        k[i][j] = res;
      }

    if (lexpr.empty()) {

      std::copy(tau_in, tau_in + ntau, lvar.get());
      precomp_array_mul_shoup(lvar.get(), ntau, mod, shoup.get());

    } else {

      unsigned nlexpr = lexpr.size();
      for (unsigned i=0; i<nlexpr; ++i) {
        UInt tmp = lexpr[i].eval(ntau, ww, tau_in, shoup.get(), mod);
        if (tmp == FAILED)
          return FAILED;
        lvar[i] = tmp;
      }
      precomp_array_mul_shoup(lvar.get(), nlexpr, mod, shoup.get());

    }

    return SUCCESS;
  }

  Ret ExternalUnitarityCut::fill_equation(const UInt x[], const UInt tau[],
                                          Mod mod,
                                          UInt term[])
  {
    UInt rhs = get_integrands(x, tau, mod);
    if (rhs == FAILED)
      return FAILED;

    term[ncoeffs()] = rhs;

    UInt * ww = ctxt->ww.get();

    for (unsigned i=0; i<ncoeffs(); ++i)
      if (xinfo()[i] && LSVar::IS_NON_ZERO) {
        const unsigned * deltavarlv = deltavar[i].get();
        unsigned lsubvarsize = deltavarlv[0];
        for (unsigned j=0; j<lsubvarsize; ++j) {
          lsubvar[j] = lvar[deltavarlv[1+j]];
          lvar_shoup[j] = shoup[deltavarlv[1+j]];
        }
        UInt res = c[i].eval(lsubvarsize, ww,
                             lsubvar.get(), lvar_shoup.get(), mod);
        if (res == FAILED)
          return FAILED;
        term[i] = res;
      }

    return SUCCESS;
  }

  UInt ExternalUnitarityCut::get_integrands(const UInt xi[], const UInt tau[],
                                            Mod mod)
  {
    UInt * ww = ctxt->ww.get();

    UInt cut = UnitarityCutFit::get_unitarity_integrand(xi, tau, mod);
    if (cut == FAILED)
      return FAILED;

    UInt tmp;
    UInt res = 0;
    const unsigned * wid = weights.get();

    const unsigned * deltavarlv = deltavar[ncoeffs()].get();
    unsigned lsubvarsize = deltavarlv[0];
    for (unsigned i=0; i<lsubvarsize; ++i) {
      lsubvar[i] = lvar[deltavarlv[1+i]];
      lvar_shoup[i] = shoup[deltavarlv[1+i]];
    }
    tmp = c[ncoeffs()].eval(lsubvarsize, ww,
                            lsubvar.get(), lvar_shoup.get(), mod);
    if (tmp == FAILED)
      return FAILED;
    res = add_mod(res, tmp, mod);

    const std::unique_ptr<unsigned[]> * integrlv = integrvar.get();
    for (const auto & cf : integr) {
      lsubvarsize = (*integrlv)[0];
      for (unsigned i=0; i<lsubvarsize; ++i) {
        lsubvar[i] = lvar[(*integrlv)[1+i]];
        lvar_shoup[i] = shoup[(*integrlv)[1+i]];
      }
      tmp = cf.eval(lsubvarsize, ww,
                    lsubvar.get(), lvar_shoup.get(), mod);
      if (tmp == FAILED)
        return FAILED;
      if (wid) {
        tmp = mul_mod(tmp, ctxt->regfun_res[*wid], mod);
        ++wid;
      }
      res = add_mod(res, tmp, mod);
      ++integrlv;
    }

    return add_mod(res, cut, mod);
  }

  std::unique_ptr<Algorithm>
  ExternalUnitarityCut::clone(ExternalContext & sctxt) const
  {
    std::unique_ptr<Algorithm> algptr(new ExternalUnitarityCut());
    auto & alg = static_cast<ExternalUnitarityCut&>(*algptr);

    alg.nparsin = nparsin;
    alg.nparsout = nparsout;
    alg.ctxt = &sctxt;

    UnitarityCutFit::copy_into(alg);

    alg.ntau = ntau;
    alg.nlvar = nlvar;
    alg.nlflav = nlflav;
    alg.max_subvar_size = max_subvar_size;
    alg.this_mod = 0;

    alg.spexpr.resize(nspinors(), 4);
    alg.spexprmap.resize(nspinors(), 4);
    alg.lvar.reset(new UInt[nlvar]);
    alg.shoup.reset(new UInt[std::max<std::size_t>(nparsin, nlvar)]);
    alg.lvar_shoup.reset(new UInt[nlvar]);

    if (weights.get()) {
      unsigned wsize = integr.size();
      alg.weights.reset(new unsigned[integr.size()]);
      std::copy(weights.get(), weights.get()+wsize, alg.weights.get());
      for (unsigned j=0; j<wsize; ++j)
        ctxt->regfun_share_with_subctxt(weights[j], nparsin, sctxt);
    }

    alg.lflav.reset(new LoopFlavours[nlflav]);
    std::copy(lflav.get(), lflav.get()+nlflav, alg.lflav.get());
    alg.lflav_coeff.reset(new MPRational[nlflav]);
    std::copy(lflav_coeff.get(), lflav_coeff.get()+nlflav,
              alg.lflav_coeff.get());

    if (max_subvar_size)
      alg.lsubvar.reset(new UInt[max_subvar_size]);

    alg.integrvar.reset(new std::unique_ptr<unsigned[]>[integr.size()]);
    for (unsigned j=0; j<integr.size(); ++j) {
      unsigned intvsize = integrvar[j][0]+1;
      alg.integrvar[j].reset(new unsigned[intvsize]);
      std::copy(integrvar[j].get(), integrvar[j].get() + intvsize,
                alg.integrvar[j].get());
    }

    alg.deltavar.reset(new std::unique_ptr<unsigned[]>[ncoeffs()+1]);
    for (unsigned j=0; j<ncoeffs()+1; ++j) {
      unsigned dvsize = deltavar[j][0]+1;
      alg.deltavar[j].reset(new unsigned[dvsize]);
      std::copy(deltavar[j].get(), deltavar[j].get() + dvsize,
                alg.deltavar[j].get());
    }

    unsigned loops = kexpr.rows();
    alg.kexprmpmap.resize(loops,6);
    alg.kexprmap.resize(loops,6);
    alg.kexpr.resize(loops,6);
    for (unsigned i=0; i<loops; ++i)
      for (unsigned j=0; j<6; ++j)
        horner_horner_ratfun_clone(kexprmpmap(i,j).data(),
                                   kexprmap(i,j),
                                   kexpr(i,j),
                                   ntau, nparsin,
                                   alg.kexprmpmap(i,j),
                                   alg.kexprmap(i,j),
                                   alg.kexpr(i,j));

    unsigned csize = c.size();
    alg.cmpmap.resize(csize);
    alg.cmap.resize(csize);
    alg.c.resize(csize);
    for (unsigned i=0; i<csize; ++i)
      horner_horner_ratfun_clone(cmpmap[i].data(), cmap[i], c[i],
                                 deltavar[i][0], nparsin,
                                 alg.cmpmap[i], alg.cmap[i], alg.c[i]);

    unsigned integrsize = integr.size();
    alg.integrmpmap.resize(integrsize);
    alg.integrmap.resize(integrsize);
    alg.integr.resize(integrsize);
    for (unsigned i=0; i<integrsize; ++i)
      horner_horner_ratfun_clone(integrmpmap[i].data(), integrmap[i], integr[i],
                                 integrvar[i][0], nparsin,
                                 alg.integrmpmap[i], alg.integrmap[i],
                                 alg.integr[i]);

    unsigned lexprsize = lexpr.size();
    alg.lexprmpmap.resize(lexprsize);
    alg.lexprmap.resize(lexprsize);
    alg.lexpr.resize(lexprsize);
    for (unsigned i=0; i<lexprsize; ++i)
      horner_horner_ratfun_clone(lexprmpmap[i].data(), lexprmap[i], lexpr[i],
                                 ntau, nparsin,
                                 alg.lexprmpmap[i], alg.lexprmap[i],
                                 alg.lexpr[i]);

    unsigned n = spexpr.rows();
    for (unsigned i=0; i<n; ++i)
      for (unsigned j=0; j<4; ++j)
        horner_mapped_ratfun_clone(spexprmap(i,j), spexpr(i,j), nparsin,
                                   alg.spexprmap(i,j), alg.spexpr(i,j));

    horner_mapped_ratfun_clone(prefmap, pref, nparsin, alg.prefmap, alg.pref);

    // note: cutinfo is not copied, but just read from parent
    alg.set_cut(cutinfo.get(), alg.lflav.get(), alg.nlflav);

    return algptr;
  }

} // namespace ampf
