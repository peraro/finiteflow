#include <fflow/alg_laurent.hh>

namespace fflow {

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

  Ret LaurentExpansion::init(Session & session, unsigned graphid,
                             LaurentExpansionData & data,
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
    order_.reset(new int[nparsout]);
    data.init_data_(alg);

    return SUCCESS;
  }

  void LaurentExpansionData::init_data_(const Algorithm & alg,
                                        unsigned initial_sample_points)
  {
    ufun_.set_algorithm(alg, initial_sample_points);
    rec_.n_checks = 0;
    rec_.n_singular = 0;
  }

  AlgorithmData::Ptr
  LaurentExpansion::clone_data(const AlgorithmData * datain) const
  {
    const auto & data = *static_cast<const LaurentExpansionData*>(datain);
    std::unique_ptr<LaurentExpansionData> newdata(new LaurentExpansionData());
    copy_data_into(data, *newdata);
    newdata->init_data_(*subgraph(), data.ufun_.cache_size());
    return std::move(newdata);
  }

  void LaurentExpansion::new_xpoint_(UInt xin[], AlgorithmData * datain) const
  {
    auto & data = *static_cast<LaurentExpansionData*>(datain);
    data.ufun_.new_xpoint(xin);
  }

  Ret LaurentExpansion::learnFromGraph(Context * ctxt,
                                       const MutAlgInput xin[],  Mod mod,
                                       const Graph * graph,
                                       SubGraphData * datain)
  {
    auto & data = *static_cast<LaurentExpansionData*>(datain);

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

    for (unsigned idx=0; idx<nout; ++idx) {

      thiele.reset();
      recfun.clear();
      ufun.idx = idx;

      Ret ret = thiele.reconstruct(ufun, mod);
      if (ret != SUCCESS)
        return ret;
      thiele.writeResult(mod, recfun);

      auto & rfdata = rfdata_[idx];
      int order = order_[idx];
      int pref = laurent_expansion_learn(recfun);
      if (recfun.getNumDegree()==0 && recfun.num(0) == 0)
        pref = order_[idx]+1;

      if (first) {
        rfdata.num_deg_ = recfun.getNumDegree();
        rfdata.den_deg_ = recfun.getDenDegree();
        rfdata.pref_ = pref;
        tot_xout += pref > order ? 0 : order-pref+1;
      } else {
        if (rfdata.num_deg_ != recfun.getNumDegree()
            || rfdata.den_deg_ != recfun.getDenDegree()
            || rfdata.pref_ != pref)
          return FAILED;
      }

    }

    if (first)
      xout_size_ = nparsout = tot_xout;
    recfun.clear();

    return SUCCESS;
  }

  Ret LaurentExpansion::evaluateFromGraph(Context * ctxt,
                                          const MutAlgInput xin[],  Mod mod,
                                          const Graph * graph,
                                          SubGraphData * datain,
                                          UInt * xout) const
  {
    auto & data = *static_cast<LaurentExpansionData*>(datain);

    unsigned nout = graph->nparsout;

    URatFunNumDenRecoHighDegs & rec = data.rec_;
    ExpAlgorithmRatFun & ufun = data.ufun_;

    rec.setX0(sample_uint(OFFSET_3, 0, mod));
    ufun.set_algorithm_data(data.data());
    ufun.set_context(ctxt);
    ufun.new_xpoint(xin[0]);

    for (unsigned idx=0; idx<nout; ++idx) {

      auto & rfdata = rfdata_[idx];
      int order = order_[idx];
      int pref = rfdata.pref_;
      unsigned num_deg = rfdata.num_deg_;
      unsigned den_deg = rfdata.den_deg_;

      if (pref > order)
        continue;

      ufun.idx = idx;

      // Note: the reconstruction assumes den[0]=1, therefore, as a
      // workaround, if we have a singular function in x=0 we multiply
      // it by a power of x such that it becomes regular.

      int eff_pref = pref > 0 ? pref : 0;
      int eff_order = pref > 0 ? order : order-pref;
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

      unsigned this_xout_size = order-pref+1;
      laurent_expansion(rec.getFunction(), eff_order, mod, eff_pref, xout);
      xout += this_xout_size;

    }

    return SUCCESS;
  }


  void LaurentExpansion::prefactor_exponent(int pref[]) const
  {
    unsigned nout = subgraph()->nparsout;
    for (unsigned idx=0; idx<nout; ++idx)
      pref[idx] = rfdata_[idx].pref_;
  }

} // namespace fflow
