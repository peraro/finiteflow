#include <fflow/subgraph_fit.hh>

namespace fflow {


  namespace detail {

    Ret SubgraphSubFit::new_xpoint(Context *, AlgInput *, Mod,
                                   AlgorithmData *) const
    {
      return SUCCESS;
    }

    Ret SubgraphSubFit::fill_equation(Context *,
                                      AlgInput *, const UInt tau[], Mod mod,
                                      AlgorithmData * datain,
                                      UInt res[]) const
    {
      auto & data = *static_cast<SubgraphSubFitData*>(datain);

      const Graph * graph = data.graph;

      for (unsigned j=0; j<n_sample_vars(); ++j)
        data.xin[j] = tau[j];

      AlgInput sub_xin = data.xin;
      return graph->evaluate(data.ctxt, &sub_xin, mod, data.data, res);
    }

  } // namespace detail



  Ret SubgraphFit::init(Session & session, unsigned graphid,
                        SubgraphFitData & data,
                        const unsigned * param_nodes,
                        unsigned n_params_nodes,
                        std::size_t ncoeffs, std::size_t nsamplevars,
                        unsigned * needed_coeffs, unsigned needed_size,
                        unsigned extra_eqs)
  {
    impl_.init(ncoeffs, nsamplevars, needed_coeffs, needed_size,
               data.impl_, extra_eqs);

    Ret ret =  SubGraph::init_subgraph(session, graphid, data,
                                       param_nodes, n_params_nodes,
                                       nsamplevars);

    if (subgraph()->nparsout != ncoeffs + 1)
      return FAILED;

    return ret;
  }

  AlgorithmData::Ptr SubgraphFit::clone_data(const AlgorithmData * datain) const
  {
    std::unique_ptr<SubgraphFitData> algptr(new SubgraphFitData());
    auto & data = *static_cast<const SubgraphFitData*>(datain);
    impl_.copy_data(&data.impl_, &algptr.get()->impl_);
    copy_data_into(data, *algptr);
    return std::move(algptr);
  }

  Ret SubgraphFit::learnFromGraph(Context * ctxt,
                                  const MutAlgInput xin[],  Mod mod,
                                  const Graph * graph,
                                  SubGraphData * datain)
  {
    auto & data = static_cast<SubgraphFitData*>(datain)->impl_;
    data.ctxt = ctxt;
    data.xin = xin[0];
    data.graph = graph;
    data.data = datain->data();

    Ret ret = impl_.learn(ctxt, xin, mod, &data);
    nparsout = impl_.nparsout;

    return ret;
  }

  Ret SubgraphFit::evaluateFromGraph(Context * ctxt,
                                     const MutAlgInput xin[],  Mod mod,
                                     const Graph * graph,
                                     SubGraphData * datain,
                                     UInt * xout) const
  {
    auto & data = static_cast<SubgraphFitData*>(datain)->impl_;
    data.ctxt = ctxt;
    data.xin = xin[0];
    data.graph = graph;
    data.data = datain->data();

    return impl_.evaluate(ctxt, xin, mod, &data, xout);
  }


  // Multi-fit

  namespace detail {

    Ret SubgraphSubMultiFit::new_xpoint(Context *, AlgInput *, Mod,
                                        AlgorithmData *) const
    {
      return SUCCESS;
    }

    Ret SubgraphSubMultiFit::fill_equation(Context *,
                                           AlgInput *, const UInt tau[],
                                           Mod mod,
                                           AlgorithmData * datain,
                                           UInt res[]) const
    {
      auto & data = *static_cast<SubgraphSubMultiFitData*>(datain);

      const Graph * graph = data.graph;

      UInt * cached_xout = nullptr;
      bool found = data.cache->find(tau, &cached_xout);

      if (!found || cached_xout[0] == INVALID_ID) {

        for (unsigned j=0; j<n_sample_vars(); ++j)
          data.xin[j] = tau[j];

        AlgInput sub_xin = data.xin;
        Ret ret = graph->evaluate(data.ctxt, &sub_xin, mod, data.data,
                                  data.gout);
        if (ret != SUCCESS)
          return ret;

        if (!found) {

          UInt * cached_xin = data.cache->get_new_inptr();
          for (unsigned j=0; j<n_sample_vars(); ++j)
            cached_xin[j] = tau[j];
          cached_xout = data.cache->new_entry(cached_xin, graph->nparsout);

        }

        for (unsigned j=0; j<graph->nparsout; ++j)
          cached_xout[j] = data.gout[j];

      }

      for (unsigned j=0; j<data.take_size; ++j)
        res[j] = cached_xout[data.take[j]];

      return SUCCESS;
    }

  } // namespace detail


  Ret SubgraphMultiFit::init(Session & session, unsigned graphid,
                             SubgraphMultiFitData & data,
                             const unsigned * param_nodes,
                             unsigned n_params_nodes,
                             std::vector<std::vector<unsigned>> && take,
                             std::size_t nsamplevars,
                             std::vector<std::vector<unsigned>> & needed,
                             unsigned extra_eqs)
  {
    impl_.resize(take.size());
    data.cache_.init(nsamplevars);

    {
      Ret ret =  SubGraph::init_subgraph(session, graphid, data,
                                         param_nodes, n_params_nodes,
                                         nsamplevars);
      if (ret != SUCCESS)
        return ret;
    }

    unsigned gnout = subgraph()->nparsout;
    data.impl_.resize(take.size());
    for (unsigned j=0; j<take.size(); ++j) {

      impl_[j].init(take[j].size()-1, nsamplevars,
                    needed[j].data(), needed[j].size(),
                    data.impl_[j], extra_eqs);

      bool bounds_ok = std::all_of(take[j].begin(), take[j].end(),
                                   [gnout](unsigned j)
                                   {
                                     return j < gnout;
                                   });
      if (!bounds_ok)
        return FAILED;
    }

    data.gout_.reset(new UInt[gnout]);
    take_ = std::move(take);

    return SUCCESS;
  }

  AlgorithmData::Ptr
  SubgraphMultiFit::clone_data(const AlgorithmData * datain) const
  {
    std::unique_ptr<SubgraphMultiFitData> algptr(new SubgraphMultiFitData());
    auto & data = *static_cast<const SubgraphMultiFitData*>(datain);
    unsigned dsize = impl_.size();
    algptr->impl_.resize(dsize);
    algptr->cache_.init(n_sample_vars());
    copy_data_into(data, *algptr);
    for (unsigned j=0; j<dsize; ++j)
      impl_[j].copy_data(&data.impl_[j], &algptr.get()->impl_[j]);
    algptr->gout_.reset(new UInt[subgraph()->nparsout]);
    return std::move(algptr);
  }

  Ret SubgraphMultiFit::learnFromGraph(Context * ctxt,
                                       const MutAlgInput xin[],  Mod mod,
                                       const Graph * graph,
                                       SubGraphData * datain)
  {
    auto & data = static_cast<SubgraphMultiFitData*>(datain)->impl_;
    auto & cache = static_cast<SubgraphMultiFitData*>(datain)->cache_;
    UInt * gout = static_cast<SubgraphMultiFitData*>(datain)->gout_.get();
    unsigned take_size = take_.size();

    cache.set_first_value_entries(INVALID_ID);

    nparsout = 0;

    for (unsigned j=0; j<take_size; ++j) {
      data[j].ctxt = ctxt;
      data[j].xin = xin[0];
      data[j].graph = graph;
      data[j].data = datain->data();
      data[j].cache = &cache;
      data[j].gout = gout;
      data[j].take = take_[j].data();
      data[j].take_size = take_[j].size();

      Ret ret = impl_[j].learn(ctxt, xin, mod, &data[j]);
      if (ret != SUCCESS)
        return ret;
      if (impl_[j].nparsout != ~unsigned(0))
        nparsout += impl_[j].nparsout;
    }

    return SUCCESS;
  }

  Ret SubgraphMultiFit::evaluateFromGraph(Context * ctxt,
                                          const MutAlgInput xin[],  Mod mod,
                                          const Graph * graph,
                                          SubGraphData * datain,
                                          UInt * xout) const
  {
    auto & data = static_cast<SubgraphMultiFitData*>(datain)->impl_;
    auto & cache = static_cast<SubgraphMultiFitData*>(datain)->cache_;
    UInt * gout = static_cast<SubgraphMultiFitData*>(datain)->gout_.get();
    unsigned take_size = take_.size();

    cache.set_first_value_entries(INVALID_ID);

    for (unsigned j=0; j<take_size; ++j) {
      data[j].ctxt = ctxt;
      data[j].xin = xin[0];
      data[j].graph = graph;
      data[j].data = datain->data();
      data[j].cache = &cache;
      data[j].gout = gout;
      data[j].take = take_[j].data();
      data[j].take_size = take_[j].size();

      Ret ret = impl_[j].evaluate(ctxt, xin, mod, &data[j], xout);
      if (ret != SUCCESS)
        return ret;

      xout += impl_[j].nparsout;
    }

    return SUCCESS;
  }


} // namespace fflow
