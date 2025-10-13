#include <fflow/subgraph_reconstruct.hh>

namespace fflow {

  void SubgraphRec::set_default_shift()
  {
    Mod mod(BIG_UINT_PRIMES[BIG_UINT_PRIMES_SIZE-1]);
    const unsigned nvars = n_rec_vars();
    if (nvars)
      shift_.reset(new UInt[nvars]);
    for (unsigned i=0; i<nvars; ++i)
      shift_[i] = sample_uint(OFFSET_4, nparsin[0] + i, mod);
  }

  Ret SubgraphRec::init(Session & session, unsigned graphid,
                        SubgraphRecData & data,
                        unsigned n_parsin, unsigned n_recvars,
                        bool shiftvars)
  {
    if (n_recvars < 2)
      return FAILED;

    Ret ret =  SubGraph::init_subgraph(session, graphid, data,
                                       &n_parsin, 1,
                                       n_recvars);
    if (ret != SUCCESS)
      return ret;

    if (shiftvars)
      set_default_shift();
    else
      shift_.reset();

    return SUCCESS;
  }


  struct SubgraphToRecAlgData_ : public AlgorithmData {

    SubgraphToRecAlgData_(UInt * in_varsin,
                          SubGraphData * in_subg_data)
      : varsin(in_varsin), subg_data(in_subg_data) {}

    UInt * varsin;
    SubGraphData * subg_data;
  };

  struct SubgraphToRecAlg_ : public Algorithm {

    SubgraphToRecAlg_(Context * in_ctxt, const Graph * in_graph,
                      const MutAlgInput * in_xin, unsigned in_rec_vars)
      : ctxt(in_ctxt), graph(in_graph), xin(in_xin)
    {
      nparsin.resize(1);
      nparsin[0] = in_rec_vars;
      nparsout = in_graph->nparsout;
    }

    virtual Ret evaluate(Context * ctxt,
                         AlgInput tauin[], Mod mod, AlgorithmData * datain,
                         UInt xout[]) const
    {
      auto & data = *static_cast<SubgraphToRecAlgData_*>(datain);
      const UInt * tau = tauin[0];
      std::copy(tau, tau+nparsin[0], data.varsin);
      return graph->evaluate(ctxt, xin, mod, data.subg_data->data(), xout);
    }

    Context * ctxt;
    const Graph * graph;
    const MutAlgInput * xin;
  };


  unsigned SubgraphRec::get_start_mod_(Mod mod)
  {
    unsigned i=0;
    for (i=0; i<BIG_UINT_PRIMES_SIZE; ++i)
      if (BIG_UINT_PRIMES[i] == mod.n())
        return i;
    return 0;
  }

  unsigned SubgraphRec::cache_size_(unsigned nparsin, unsigned nparsout,
                                    const std::unique_ptr<UInt[]> * start,
                                    const std::unique_ptr<UInt[]> * end)
  {
    unsigned flags_size = bit_array_u64size(nparsout);
    unsigned tot = (nparsin+1+flags_size)*(end-start);
    for (; start!=end; ++start) {
      ConstBitArrayView arrv((*start).get()+nparsin+1);
      tot += bit_array_nonzeroes(arrv, nparsout);
    }
    return tot;
  }

  void SubgraphRec::gen_sample_points_(SamplePointsVector & samples,
                                       unsigned nin, unsigned nout,
                                       const ReconstructionOptions & opt) const
  {
    GenerateSamplePoints sfun(nin);
    algorithm_generate_sample_points(sfun, nin, nout,
                                     shift_.get(), opt,
                                     degs_.vdegs.get(),
                                     degs_.numdeg.get(), degs_.dendeg.get());
    unsigned flags_size = bit_array_u64size(nout);
    sfun.append_to_vector(samples, flags_size);
  }

  void SubgraphRec::remove_unneeded_points_(SamplePointsVector & samples,
                                            unsigned nin, unsigned nout,
                                            const ReconstructionOptions & opt) const
  {
    const UInt this_mod = BIG_UINT_PRIMES[opt.start_mod];
    auto nend = std::remove_if(samples.begin(), samples.end(),
                               [nin,nout,this_mod]
                               (const std::unique_ptr<UInt[]> & el)
                               {
                                 if (el[nin] != this_mod)
                                   return true;
                                 ConstBitArrayView av(el.get()+nin+1);
                                 return bit_array_nonzeroes(av, nout) == 0;
                               });
    samples.erase(nend, samples.end());
  }

  void SubgraphRec::verify_sample_points_(SamplePointsVector & samples,
                                          unsigned nin, unsigned nout,
                                          const ReconstructionOptions & opt)
  {
    // verify
    std::vector<unsigned> to_rec(nout);
    std::iota(to_rec.begin(), to_rec.end(), 0);
    algorithm_verify_sample_points(samples.data(), samples.size(),
                                   nin, nout,
                                   to_rec,
                                   shift_.get(), opt,
                                   degs_.vdegs.get(),
                                   degs_.numdeg.get(), degs_.dendeg.get());
    nsamples_ = samples.size();

    // remove unneeded samples and samples with different primes
    remove_unneeded_points_(samples, nin, nout, opt);
    nfiltered_samples_ = samples.size();
  }

  Ret SubgraphRec::filter_points_(SamplePointsVector & samples,
                                  unsigned nin,unsigned nout,
                                  const ReconstructionOptions & opt) const
  {
    if (samples.size() != nsamples_)
      return FAILED;

    // verify
    std::vector<unsigned> to_rec(nout);
    std::iota(to_rec.begin(), to_rec.end(), 0);
    algorithm_verify_sample_points(samples.data(), samples.size(),
                                   nin, nout,
                                   to_rec,
                                   shift_.get(), opt,
                                   degs_.vdegs.get(),
                                   degs_.numdeg.get(), degs_.dendeg.get());

    remove_unneeded_points_(samples, nin, nout, opt);

    if (samples.size() != nfiltered_samples_)
      return FAILED;

    return SUCCESS;
  }


  Ret SubgraphRec::learnFromGraph(Context * ctxt,
                                  const MutAlgInput xin[],  Mod mod,
                                  const Graph * graph,
                                  SubGraphData * data)
  {
    unsigned nrec = n_rec_vars();
    auto a = SubgraphToRecAlg_(ctxt,graph,xin,nrec);
    auto sdata = SubgraphToRecAlgData_(xin[0], data);

    unsigned nout = subgraph()->nparsout;
    ReconstructionOptions opt = opt_;
    opt.max_primes = 1;
    opt.start_mod = get_start_mod_(mod);


    // compute total degrees

    degs_.numdeg.reset(new unsigned[nout]);
    degs_.dendeg.reset(new unsigned[nout]);

    Ret ret = algorithm_get_degrees(a, &sdata, ctxt,
                                    shift_.get(), opt,
                                    degs_.numdeg.get(), degs_.dendeg.get());
    if (ret != SUCCESS)
      return ret;


    // compute partial degrees

    degs_.vdegs.reset(new RatFunVarDegrees[nout]);
    for (unsigned i=0; i<nout; ++i)
      degs_.vdegs[i].resize(nrec);

    for (unsigned var=0; var<nrec; ++var) {
      RatFunVarDegrees * vdegs = degs_.vdegs.get();
      ret = algorithm_get_var_degrees(a, &sdata, ctxt,
                                      var, shift_.get(), opt,
                                      vdegs);
      if (ret != SUCCESS)
        return ret;
    }


    // do a full reconstruction, saving relevant info

    SamplePointsVector samples;
    gen_sample_points_(samples, nrec, nout, opt);
    verify_sample_points_(samples, nrec, nout, opt);

    UIntCache cache;
    cache.init(nrec+1);
    alloc_size_ = cache_size_(nrec, a.nparsout,
                              samples.data(), samples.data()+samples.size());
    cache.reserve_more(samples.size(), alloc_size_);

    evaluate_and_sparse_cache_samples(a, &sdata, ctxt,
                                      samples.data(),
                                      samples.data() + samples.size(),
                                      cache);

    res_.resize(nout);

    nparsout = 0;
    for (unsigned i=0; i<nout; ++i) {
      ret = algorithm_sparse_reconstruct_mod(mod,
                                             cache, nrec, nout,
                                             i,
                                             shift_.get(), opt,
                                             degs_.vdegs.get(),
                                             degs_.numdeg.get(),
                                             degs_.dendeg.get(),
                                             res_[i]);
      if (ret != SUCCESS) {
        nparsout = 0;
        return ret;
      }
      nparsout += res_[i].numerator().size() + res_[i].denominator().size();
    }

    return SUCCESS;
  }


  Ret SubgraphRec::evaluateFromGraph(Context * ctxt,
                                     const MutAlgInput xin[],  Mod mod,
                                     const Graph * graph, SubGraphData * data,
                                     UInt * xout) const
  {
    unsigned nrec = n_rec_vars();
    auto a = SubgraphToRecAlg_(ctxt,graph,xin,nrec);
    auto sdata = SubgraphToRecAlgData_(xin[0], data);

    unsigned nout = subgraph()->nparsout;
    ReconstructionOptions opt = opt_;
    opt.max_primes = 1;
    opt.start_mod = get_start_mod_(mod);

    SamplePointsVector samples;
    gen_sample_points_(samples, nrec, nout, opt);
    Ret ret = filter_points_(samples, nrec, nout, opt);
    if (ret != SUCCESS)
      return ret;

    UIntCache cache;
    cache.init(nrec+1);
    cache.reserve_more(samples.size(), alloc_size_);

    evaluate_and_sparse_cache_samples(a, &sdata, ctxt,
                                      samples.data(),
                                      samples.data() + samples.size(),
                                      cache);

    std::vector<SparseRationalFunction> res(nout);

    unsigned nout_count = 0;
    for (unsigned i=0; i<nout; ++i) {
      ret = algorithm_sparse_reconstruct_mod(mod,
                                             cache, nrec, nout,
                                             i,
                                             shift_.get(), opt,
                                             degs_.vdegs.get(),
                                             degs_.numdeg.get(),
                                             degs_.dendeg.get(),
                                             res[i]);
      if (ret != SUCCESS)
        return ret;
      nout_count += res[i].numerator().size() + res[i].denominator().size();
    }

    if (nout_count != nparsout)
      return FAILED;

    unsigned idx = 0;
    for (unsigned i=0; i<nout; ++i) {
      for (const auto & mon : res[i].numerator())
        xout[idx++] = mon.coeff();
      for (const auto & mon : res[i].denominator())
        xout[idx++] = mon.coeff();
    }

    return SUCCESS;
  }


  AlgorithmData::Ptr
  SubgraphRec::clone_data(const AlgorithmData * datain) const
  {
    auto * data = static_cast<const SubgraphRecData*>(datain);
    std::unique_ptr<SubgraphRecData> ptr(new SubgraphRecData());
    copy_data_into(*data, *ptr);
    return std::move(ptr);
  }



  // Univariate

  void ExpAlgorithmRatFun::new_xpoint(UInt x[])
  {
    xin_ = x;
    cache_.set_first_value_entries(INVALID_ID);
  }

  Ret ExpAlgorithmRatFun::evaluate(const UInt x, Mod mod)
  {
    if (!x_)
      x_ = cache_.get_new_inptr();

    xin_[0] = *x_ = x;
    UInt * xout = nullptr;
    bool found = cache_.find(x_, &xout);

    if (!found || xout[0] == INVALID_ID) {

      if (!found)
        xout = cache_.new_entry(x_,nparsout_);

      const UInt * xin = xin_;
      Ret ret = alg_->evaluate(ctxt_, &xin, mod, alg_data_, xout);
      if (ret == FAILED)
        xout[0] = FAILED;
      if (!found)
        x_ = nullptr;

    }

    if (xout[0] == FAILED)
      return FAILED;

    return mul_mod(xout[idx], power(x, pref_exp_, mod), mod);
  }

  Ret BaseSubgraphUniRec::init(Session & session, unsigned graphid,
                               BaseSubgraphUniRecData & data,
                               const unsigned param_nodes[],
                               unsigned n_params_nodes)
  {
    Ret ret = init_subgraph(session, graphid, data,
                            param_nodes, n_params_nodes, 1);
    if (ret != SUCCESS)
      return ret;

    const Algorithm & alg = *subgraph();
    const unsigned nparsout = alg.nparsout;
    rfdata_.reset(new ExpRatFunData[nparsout]);
    if (laurent_expand_)
      order_.reset(new int[nparsout]);
    data.init_data_(alg);

    return SUCCESS;
  }

  void BaseSubgraphUniRecData::init_data_(const Algorithm & alg,
                                          unsigned initial_sample_points)
  {
    ufun_.set_algorithm(alg, initial_sample_points);
    rec_.n_checks = 0;
    rec_.n_singular = 0;
  }

  AlgorithmData::Ptr
  BaseSubgraphUniRec::clone_data(const AlgorithmData * datain) const
  {
    const auto & data = *static_cast<const BaseSubgraphUniRecData*>(datain);
    std::unique_ptr<BaseSubgraphUniRecData> newdata(new BaseSubgraphUniRecData());
    copy_data_into(data, *newdata);
    newdata->init_data_(*subgraph(), data.ufun_.cache_size());
    return std::move(newdata);
  }

  void BaseSubgraphUniRec::new_xpoint_(UInt xin[], AlgorithmData * datain) const
  {
    auto & data = *static_cast<BaseSubgraphUniRecData*>(datain);
    data.ufun_.new_xpoint(xin);
  }

  Ret BaseSubgraphUniRec::learnFromGraph(Context * ctxt,
                                         const MutAlgInput xin[],  Mod mod,
                                         const Graph * graph,
                                         SubGraphData * datain)
  {
    auto & data = *static_cast<BaseSubgraphUniRecData*>(datain);

    unsigned nout = graph->nparsout;

    URationalFunction & recfun = data.rec_.getFunction();
    ExpAlgorithmRatFun & ufun = data.ufun_;

    ufun.set_algorithm_data(data.data());
    ufun.set_context(ctxt);
    ufun.new_xpoint(xin[0]);

    URatFunReconstruction thiele(0, 2*max_degree+1);
    thiele.setX0(sample_uint(OFFSET_3, 0, mod));

    unsigned tot_xout = 0;
    const bool first = xout_size_ == ~unsigned(0);

    if (!laurent_expand_)
      ratrec_.reset(new SparseRationalFunction[nout]);

    for (unsigned idx=0; idx<nout; ++idx) {

      thiele.reset();
      recfun.clear();
      ufun.idx = idx;

      Ret ret = thiele.reconstruct(ufun, mod);
      if (ret != SUCCESS)
        return ret;
      thiele.writeResult(mod, recfun);

      auto & rfdata = rfdata_[idx];
      int order = laurent_expand_ ? order_[idx] : 0;
      int pref = laurent_expansion_learn(recfun); // needed in both cases
      if (laurent_expand_ && recfun.getNumDegree()==0 && recfun.num(0) == 0)
        pref = order_[idx]+1;

      if (first) {
        rfdata.num_deg_ = recfun.getNumDegree();
        rfdata.den_deg_ = recfun.getDenDegree();
        rfdata.pref_ = pref;
      } else {
        if (rfdata.num_deg_ != recfun.getNumDegree()
            || rfdata.den_deg_ != recfun.getDenDegree()
            || rfdata.pref_ != pref)
          return FAILED;
      }

      if (laurent_expand_) {

        tot_xout += pref > order ? 0 : order-pref+1;

      } else {

        auto & rrec = ratrec_[idx];
        rrec = SparseRationalFunction(1);
        rrec.fromUnivariate(recfun);
        if (pref > 0)
          for (auto & mon : rrec.numerator()) {
            mon.exponent(0) += pref;
            mon.degree() += pref;
          }
        else if (pref < 0)
          for (auto & mon : rrec.denominator()) {
            mon.exponent(0) -= pref;
            mon.degree() -= pref;
          }
        tot_xout += rrec.numerator().size() + rrec.denominator().size();

      }

    }

    if (first)
      xout_size_ = nparsout = tot_xout;
    recfun.clear();

    return SUCCESS;
  }

  Ret BaseSubgraphUniRec::evaluateFromGraph(Context * ctxt,
                                            const MutAlgInput xin[],  Mod mod,
                                            const Graph * graph,
                                            SubGraphData * datain,
                                            UInt * xout_in) const
  {
    auto & data = *static_cast<BaseSubgraphUniRecData*>(datain);
    UInt * xout = xout_in;

    unsigned nout = graph->nparsout;

    URatFunNumDenRecoHighDegs & rec = data.rec_;
    ExpAlgorithmRatFun & ufun = data.ufun_;

    rec.setX0(sample_uint(OFFSET_3, 0, mod));
    ufun.set_algorithm_data(data.data());
    ufun.set_context(ctxt);
    ufun.new_xpoint(xin[0]);

    for (unsigned idx=0; idx<nout; ++idx) {

      auto & rfdata = rfdata_[idx];
      int order = laurent_expand_ ? order_[idx] : 0;
      int pref = rfdata.pref_;
      unsigned num_deg = rfdata.num_deg_;
      unsigned den_deg = rfdata.den_deg_;

      if (laurent_expand_ && pref > order)
        continue;

      ufun.idx = idx;

      // Note: the reconstruction assumes den[0]=1, therefore, as a
      // workaround, if we have a singular function in x=0 we multiply
      // it by a power of x such that it becomes regular.

      rec.setNumDegree(num_deg);
      rec.setDenDegree(den_deg);
      if (pref > 0) {
        rec.rnum_min = pref;
        rec.rden_min = 1;
        UInt * coeffs = rec.getFunction().num().coeffs();
        std::fill(coeffs, coeffs+pref, 0);
        rec.getFunction().den().coeffs()[0] = 1;
        ufun.set_pref_exp(0);
      } else {
        unsigned den_pref = -pref;
        rec.rnum_min = 0;
        rec.rden_min = 1;
        rec.getFunction().den().coeffs()[0] = 1;
        ufun.set_pref_exp(den_pref);
        rec.setDenDegree(den_deg-den_pref);
      }

      Ret ret = rec.reconstruct(ufun, mod);
      if (ret != SUCCESS)
        return ret;

      if (laurent_expand_) {

        int eff_pref = pref > 0 ? pref : 0;
        int eff_order = pref > 0 ? order : order-pref;
        unsigned this_xout_size = order-pref+1;
        laurent_expansion(rec.getFunction(), eff_order, mod, eff_pref, xout);
        xout += this_xout_size;

      } else {

        const auto & rrec = ratrec_[idx];
        const auto & recfun = rec.getFunction();

        const unsigned num_pref = pref > 0 ? pref : 0;
        const unsigned num_deg = recfun.getNumDegree()+num_pref;
        if (num_deg != rrec.numerator().degree())
          return FAILED;

        const unsigned den_pref = pref < 0 ? -pref : 0;
        const unsigned den_deg = recfun.getDenDegree()+den_pref;
        if (den_deg != rrec.denominator().degree())
          return FAILED;

        for (const auto & mon : rrec.numerator()) {
          if (xout - xout_in >= nparsout)
            return FAILED;
          *xout = recfun.num()[mon.degree()];
          ++xout;
        }

        for (const auto & mon : rrec.denominator()) {
          if (xout - xout_in >= nparsout)
            return FAILED;
          *xout = recfun.den()[mon.degree()];
          ++xout;
        }

      }

    }

    return SUCCESS;
  }


  void BaseSubgraphUniRec::prefactor_exponent(int pref[]) const
  {
    unsigned nout = subgraph()->nparsout;
    for (unsigned idx=0; idx<nout; ++idx)
      pref[idx] = rfdata_[idx].pref_;
  }


} // namespace fflow
