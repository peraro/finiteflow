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


} // namespace fflow
