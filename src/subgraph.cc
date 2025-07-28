#include <fflow/subgraph.hh>

namespace fflow {

  Ret SubGraph::init_subgraph(Session & session, unsigned graphid,
                              SubGraphData & data,
                              const unsigned npars[], unsigned npars_size,
                              unsigned n_aux_vars)
  {
    Graph * g = session.graph(graphid);
    if (!g || !g->has_learned()) {
      logerr("Subgraph does not exist or cannot be evaluated");
      return FAILED;
    }

    nparsin.resize(npars_size);
    unsigned tot_parsin = g->nparsin.size() ? g->nparsin[0] : 0;

    if (!npars_size && tot_parsin != n_aux_vars) {
      logerr("Inputs of subgraph node do not have the expected length");
      return FAILED;
    }

    for (unsigned j=0; j<npars_size; ++j) {
      unsigned nout = npars[j];
      nparsin[j] = nout;
      if (tot_parsin != nparsin[j] + n_aux_vars) {
        logerr("Inputs of subgraph node do not have the expected length");
        return FAILED;
      }
    }

    aux_vars_ = n_aux_vars;
    nparsout = g->nparsout;
    dec_subgraph_count_();
    graph_.reset(g);
    g->set_immutable();
    g->subgraph_count_ += 1;
    if (g->out_node() && g->out_node()->algorithm())
      g->out_node()->algorithm()->set_immutable();

    data.data_ = g->clone_data(session.ctxt_.graph_data(graphid));
    set_xin_(data);

    return SUCCESS;
  }

  void SubGraph::set_xin_(SubGraphData & data) const
  {
    unsigned tot_parsin = 0;
    unsigned npars_size = nparsin.size();
    for (unsigned j=0; j<npars_size; ++j)
      tot_parsin += nparsin[j] + aux_vars_;
    if (!npars_size)
      tot_parsin = aux_vars_;
    data.xin_.reset(new MutAlgInput[std::max(npars_size,unsigned(1))]);
    data.xinval_.reset(new UInt[tot_parsin]);
    MutAlgInput * xiptr = data.xin_.get();
    UInt * xivptr = data.xinval_.get();
    xiptr[0] = xivptr;
    for (unsigned j=0; j<npars_size; ++j) {
      xiptr[j] = xivptr;
      xivptr += nparsin[j] + aux_vars_;
    }
  }

  void SubGraph::chain_xin_(AlgInput xin[], MutAlgInput sub_xin[]) const
  {
    const unsigned nps = nparsin.size();
    unsigned n_aux = aux_vars_;
    for (unsigned i=0; i<nps; ++i) {
      const unsigned nps_i = nparsin[i];
      for (unsigned j=0; j<nps_i; ++j)
        sub_xin[i][j+n_aux] = xin[i][j];
    }
  }

  Ret SubGraph::learn(Context * ctxt, Input * xin, Mod mod,
                      AlgorithmData * datain)
  {
    auto & data = *static_cast<SubGraphData*>(datain);

    MutAlgInput * sub_xin = data.xin_.get();
    chain_xin_(xin, sub_xin);

    return learnFromGraph(ctxt, sub_xin, mod, graph_.get(), &data);
  }

  Ret SubGraph::learnFromGraph(Context *, const MutAlgInput *, Mod,
                               const Graph *, SubGraphData *)
  {
    return SUCCESS;
  }

  Ret SubGraph::evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * datain,
                         UInt xout[]) const
  {
    auto & data = *static_cast<SubGraphData*>(datain);

    MutAlgInput * sub_xin = data.xin_.get();
    chain_xin_(xin, sub_xin);

    return evaluateFromGraph(ctxt, sub_xin, mod, graph_.get(), &data,
                             xout);
  }

  void SubGraph::copy_data_into(const SubGraphData & datain,
                                SubGraphData & dataout) const
  {
    dataout.data_ = graph_->clone_data(datain.data_.get());
    set_xin_(dataout);
  }

  void SubGraph::dec_subgraph_count_()
  {
    if (graph_.get()) {
      graph_->subgraph_count_ -= 1;
      if (graph_->subgraph_count_ == 0)
        graph_->set_mutable();
    }
  }

  SubGraph::~SubGraph()
  {
    dec_subgraph_count_();
  }

  AlgorithmData::Ptr SubGraph::clone_data(const AlgorithmData * datain) const
  {
    auto * data = static_cast<const SubGraphData*>(datain);
    std::unique_ptr<SubGraphData> ptr(new SubGraphData());
    copy_data_into(*data, *ptr);
    return std::move(ptr);
  }


  Ret SimpleSubGraph::init(Session & session, unsigned graphid,
                           SubGraphData & data,
                           const unsigned * param_nodes,
                           unsigned n_params_nodes)
  {
    return SubGraph::init_subgraph(session, graphid, data,
                                   param_nodes, n_params_nodes);
  }

  Ret SimpleSubGraph::evaluateFromGraph(Context * ctxt,
                                        const MutAlgInput * xin, Mod mod,
                                        const Graph * graph,
                                        SubGraphData * data,
                                        UInt * xout) const
  {
    return graph->evaluate(ctxt, xin, mod, data->data(), xout);
  }


  Ret MemoizedSubGraph::init(Session & session, unsigned graphid,
                             MemoizedSubGraphData & data,
                             const unsigned * param_nodes,
                             unsigned n_params_nodes)
  {
    Ret ret = SubGraph::init_subgraph(session, graphid, data,
                                      param_nodes, n_params_nodes);
    return ret;
  }

  AlgorithmData::Ptr
  MemoizedSubGraph::clone_data(const AlgorithmData * datain) const
  {
    auto * data = static_cast<const MemoizedSubGraphData*>(datain);
    std::unique_ptr<MemoizedSubGraphData> ptr(new MemoizedSubGraphData());
    SubGraph::copy_data_into(*data, *ptr);
    return std::move(ptr);
  }

  Ret MemoizedSubGraph::evaluateFromGraph(Context * ctxt,
                                          const MutAlgInput * xin, Mod mod,
                                          const Graph * graph,
                                          SubGraphData * datain,
                                          UInt * xout) const
  {
    auto & data = *static_cast<MemoizedSubGraphData*>(datain);
    const unsigned nin = graph->nparsin.size() ? graph->nparsin[0] : 0;
    const unsigned nout = nparsout;

    if (!data.xin_.get()) {

      data.xin_.reset(new UInt[nin]);
      data.xout_.reset(new UInt[nout]);

    } else if (mod.n() == data.mod_ &&
               std::equal(xin[0], xin[0]+nin, data.xin_.get())) {

      std::copy(data.xout_.get(), data.xout_.get()+nout, xout);
      return SUCCESS;

    }

    Ret ret = graph->evaluate(ctxt, xin, mod, data.data(), xout);
    data.mod_ = mod.n();
    std::copy(xin[0], xin[0]+nin, data.xin_.get());
    std::copy(xout, xout+nout, data.xout_.get());
    return ret;
  }


  Ret SubGraphMap::init(Session & session, unsigned graphid,
                        SubGraphData & data,
                        const unsigned * param_nodes,
                        unsigned n_params_nodes)
  {
    Ret ret = SubGraph::init_subgraph(session, graphid, data,
                                      param_nodes, n_params_nodes);
    nparsout = nparsin.size() ? nparsin.size()*subgraph()->nparsout : 0;
    return ret;
  }

  Ret SubGraphMap::evaluateFromGraph(Context * ctxt,
                                     const MutAlgInput * xin, Mod mod,
                                     const Graph * graph,
                                     SubGraphData * data,
                                     UInt * xout) const
  {
    unsigned nin = nparsin.size();
    unsigned nout = graph->nparsout;
    for (unsigned j=0; j<nin; ++j) {
      Ret ret = graph->evaluate(ctxt, xin+j, mod, data->data(), xout);
      if (ret != SUCCESS)
        return ret;
      xout += nout;
    }
    return SUCCESS;
  }

} // namespace fflow
