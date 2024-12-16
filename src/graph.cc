#include <thread>
#include <chrono>
#include <fflow/graph.hh>

namespace fflow {

  Node::~Node()
  {
    if (graph_ && graph_->session_)
      graph_->session_->delete_node_data_(graph_->id_, id_);
  }

  const Session * Node::session() const
  {
    if (!graph_)
      return nullptr;
    return graph_->session_;
  }

  Session * Node::session()
  {
    if (!graph_)
      return nullptr;
    return graph_->session_;
  }


  void GraphData::add_data_(unsigned id, std::nullptr_t)
  {
    add_data_(id, AlgorithmData::Ptr(nullptr));
  }

  void GraphData::add_data_(unsigned id, AlgorithmData::Ptr && data)
  {
    if (algdata_.size() < id+1) {
      algdata_.resize(id+1);
      xout_.resize(id+1);
    }
    algdata_[id] = std::move(data);
  }

  void GraphData::set_nparsout_(unsigned id, unsigned nvars)
  {
    if (nvars)
      xout_[id].resize(nvars);
  }

  void GraphData::delete_data_(unsigned id)
  {
    if (id < algdata_.size())
      algdata_[id].reset(nullptr);
  }


  Graph::~Graph()
  {
  }

  GraphData & Graph::graph_data_()
  {
    return *(session_->ctxt_.graph_data_[id_]);
  }

  unsigned Graph::new_node(Algorithm * alg,
                           std::unique_ptr<AlgorithmData> && algdata,
                           const unsigned * inputs)
  {
    AlgorithmPtr algptr(alg);
    return new_node(std::move(algptr), std::move(algdata), inputs);
  }

  unsigned Graph::new_node(AlgorithmPtr && alg,
                           std::unique_ptr<AlgorithmData> && algdata,
                           const unsigned * inputs)
  {
    if (!is_mutable())
      return ALG_NO_ID;

    // check number inputs is consistent
    unsigned ninputs = (*alg).nparsin.size();
    for (unsigned j=0; j<ninputs; ++j)
      if (!node_exists(inputs[j])
          || nodes_[inputs[j]]->alg_->nparsout != alg->nparsin[j]) {
        return ALG_NO_ID;
      }

    // allocate new node
    unsigned id = get_free_node_id_();
    Node * nodeptr = new Node();
    nodes_[id].reset(nodeptr);
    nodeptr->id_ = id;
    nodeptr->graph_ = this;

    // set inputs, immutability and compute depth
    nodeptr->inputs_.resize(ninputs);
    unsigned depth = 0;
    for (unsigned j=0; j<ninputs; ++j) {
      Node & node = *nodes_[inputs[j]].get();
      nodeptr->inputs_[j] = node.id_;
      node.alg_->set_immutable();
      depth = std::max(depth, node.depth_);
    }
    nodes_[id]->depth_ = depth+1;

    // set the algorithm and algorithm data
    nodes_[id]->alg_ = std::move(alg);
    if (session_) {
      graph_data_().add_data_(id, std::move(algdata));
      auto & ctxt_xin = session_->ctxt_.xin_;
      if (ctxt_xin.size() < ninputs)
        ctxt_xin.resize(ninputs);
    }

    // if learning is not needed, we already set the output
    if (nodeptr->alg_->min_learn_times() == 0) {
      nodeptr->alg_->set_learned();
      if (session_)
        graph_data_().set_nparsout_(id, nodeptr->alg_->nparsout);
    }

    return id;
  }

  unsigned Graph::delete_node(unsigned id)
  {
    return delete_node_(id, false);
  }

  unsigned Graph::delete_node_(unsigned id, bool force)
  {
    Node * n = node(id);
    if (!n || (!force && !n->alg_->is_mutable()))
      return ALG_NO_ID;
    if (id == out_node_)
      out_node_ = ALG_NO_ID;
    nodes_[id].reset();
    free_node_slots_.push_back(id);
    return id;
  }

  unsigned Graph::get_free_node_id_()
  {
    unsigned id=0;
    if (!free_node_slots_.empty()) {

      id = free_node_slots_.back();
      free_node_slots_.pop_back();

    } else {

      if (input_vars_not_set())
        set_input_vars(0);
      id = nodes_.size();
      nodes_.push_back(nullptr);

    }

    return id;
  }

  unsigned Graph::edges(std::vector<unsigned> & out) const
  {
    return edges_(out, false);
  }

  unsigned Graph::marked_edges(std::vector<unsigned> & out) const
  {
    return edges_(out, true);
  }

  unsigned Graph::edges_(std::vector<unsigned> & out,
                         bool marked_only) const
  {
    out.clear();

    if (marked_only && out_node_ == ALG_NO_ID)
      return out_node_;

    for (auto & n : nodes_) {
      if (!n.get())
        continue;
      unsigned nodeid = n->id();
      if (marked_only && !mark_nodes_[nodeid])
        continue;
      for (auto inputid : n->inputs_) {
        out.push_back(inputid);
        out.push_back(nodeid);
      }
    }
    return out_node_;
  }

  unsigned Graph::nodes(std::vector<unsigned> & out) const
  {
    return get_nodes_(out, false);
  }

  unsigned Graph::marked_nodes(std::vector<unsigned> & out) const
  {
    return get_nodes_(out, true);
  }

  unsigned Graph::get_nodes_(std::vector<unsigned> & out,
                             bool marked_only) const
  {
    out.clear();

    if (marked_only && out_node_ == ALG_NO_ID)
      return out_node_;

    for (auto & n : nodes_) {
      if (!n.get())
        continue;
      unsigned nodeid = n->id();
      if (marked_only && !mark_nodes_[nodeid])
        continue;
      out.push_back(nodeid);
    }
    return out_node_;
  }

  unsigned Graph::set_input_vars(unsigned nvars)
  {
    if (!is_mutable())
      return ALG_NO_ID;

    if (nodes_.empty()) {
      nodes_.push_back(std::unique_ptr<Node>(new Node()));
      nodes_[0]->alg_.reset(new InputAlgorithm());
    } else if (!(nodes_[0]->alg_->is_mutable())) {
      return ALG_NO_ID;
    }

    if (!(nodes_[0]->alg_->is_mutable()))
      return ALG_NO_ID;

    nodes_[0]->alg_->nparsout = nvars;
    nodes_[0]->depth_ = 0;
    nodes_[0]->id_ = 0;

    nparsin.resize(1);
    nparsin[0] = nvars;

    if (session_) {
      GraphData & graphdata = graph_data_();
      graphdata.add_data_(0, nullptr);
      graphdata.set_nparsout_(0, nvars);
    }

    return 0;
  }

  void Graph::mark_node_and_deps_(unsigned id)
  {
    mark_node_(id);
    for (auto dep : nodes_[id]->inputs_) {
      if (!node_is_marked_(dep))
        mark_node_and_deps_(dep);
    }
  }

  unsigned Graph::set_output_node(unsigned id)
  {
    if (!is_mutable())
      return ALG_NO_ID;

    if (!node_exists(id))
      return ALG_NO_ID;

    if (id != out_node_) {

      reset_node_marks_();
      mark_node_and_deps_(id);
      comput_seq_.clear();
      for (unsigned j=0; j<nodes_.size(); ++j)
        if (node_is_marked_(j)) {
          if (j != id && !nodes_[j]->alg_->has_learned())
            return ALG_NO_ID;
          comput_seq_.push_back(j);
        }
      const auto * nodes = nodes_.data();
      std::sort(comput_seq_.begin(), comput_seq_.end(),
                [nodes] (unsigned i, unsigned j)
                {
                  return nodes[i]->depth_ < nodes[j]->depth_;
                });

      out_node_ = id;
      set_learned(nodes[id]->alg_->has_learned());
      nparsout = nodes[id]->alg_->nparsout;

      rec_data_.reset(nullptr);
    }

    return id;
  }

  void Graph::invalidate_reconstruction_cache()
  {
    rec_data_.reset(nullptr);
  }

  Ret Graph::evaluate_comput_seq_(Context * ctxt,
                                  Input xinin[], Mod mod,
                                  AlgorithmData * datain,
                                  bool skip_out_node) const
  {
    GraphData & data = *static_cast<GraphData*>(datain);

    const UInt * xin = xinin[0];
    std::copy(xin, xin+nparsin[0], data.xout_[0].get());

    UInt ** ctxt_xin = ctxt->xin_.data();

    for (auto algid : comput_seq_) {

      if (skip_out_node && algid == out_node_)
        continue;

      Node & node = *nodes_[algid];
      const Algorithm * alg = node.alg_.get();
      if (alg->is_input_alg())
        continue;

      unsigned algnin = node.inputs_.size();
      for (unsigned j=0; j<algnin; ++j)
        ctxt_xin[j] = data.xout_[node.inputs_[j]].get();

      Ret ret = alg->evaluate(ctxt, ctxt_xin, mod,
                              data.algdata_[algid].get(),
                              data.xout_[algid].get());
      if (ret != SUCCESS)
        return ret;
    }

    return SUCCESS;
  }

  Ret Graph::evaluate(Context * ctxt,
                      Input xinin[], Mod mod, AlgorithmData * datain,
                      UInt xout[]) const
  {
    Ret ret = evaluate_comput_seq_(ctxt, xinin, mod, datain, false);
    if (ret != SUCCESS)
      return ret;

    GraphData & data = *static_cast<GraphData*>(datain);
    const UInt * alg_xout = data.xout_[out_node_].get();
    std::copy(alg_xout, alg_xout+nparsout, xout);

    return SUCCESS;
  }

  UInt Graph::min_learn_times()
  {
    if (out_node_ == 0 || out_node_ == ALG_NO_ID)
      return 0;
    return nodes_[out_node_]->alg_->min_learn_times();
  }

  Ret Graph::learn(Context * ctxt,
                   Input xinin[], Mod mod, AlgorithmData * datain)
  {
    Ret ret = evaluate_comput_seq_(ctxt, xinin, mod, datain, true);
    if (ret != SUCCESS)
      return ret;

    GraphData & data = *static_cast<GraphData*>(datain);
    UInt ** ctxt_xin = ctxt->xin_.data();
    const unsigned algid = out_node_;
    Node & node = *nodes_[algid];
    Algorithm * alg = node.alg_.get();

    if (!alg->is_input_alg()) {
      unsigned algnin = node.inputs_.size();
      for (unsigned j=0; j<algnin; ++j)
        ctxt_xin[j] = data.xout_[node.inputs_[j]].get();
      ret = alg->learn(ctxt, ctxt_xin, mod, data.algdata_[algid].get());
      if (ret != SUCCESS)
        return ret;
    }

    nparsout = nodes_[out_node_]->alg_->nparsout;
    return SUCCESS;
  }

  AlgorithmData::Ptr Graph::clone_data(const AlgorithmData * datain) const
  {
    std::unique_ptr<GraphData> newdata(new GraphData());
    sync_graph_data(*static_cast<const GraphData*>(datain), *newdata, true);
    return std::move(newdata);
  }

  void Graph::sync_graph_data(const GraphData & adatain,
                              GraphData & adataout, bool all_nodes) const
  {
    adataout.algdata_.resize(adatain.algdata_.size());
    adataout.xout_.resize(adatain.xout_.size());

    unsigned idx = 0;
    for (const auto & node : nodes_) {
      if (!node.get()) {
        ++idx;
        continue;
      }
      Algorithm & alg = *node->algorithm();
      if (all_nodes || mark_nodes_[idx])
        if (adataout.algdata_[idx].get() == nullptr) {
          adataout.algdata_[idx] = alg.clone_data(adatain.algdata_[idx].get());
          adataout.xout_[idx].resize(alg.nparsout);
        }
      ++idx;
    }
  }

  Graph::AlgRecData & Graph::rec_data()
  {
    return *rec_data_;
  }

  Graph::AlgDegs & Graph::degs_data()
  {
    return rec_data_->degs;
  }

  const Graph::AlgRecData & Graph::rec_data() const
  {
    return *rec_data_;
  }

  const Graph::AlgDegs & Graph::degs_data() const
  {
    return rec_data_->degs;
  }


  bool Context::has_graph_data(unsigned id) const
  {
    return (graph_data_.size() >= id + 1) && (graph_data_[id].get() != nullptr);
  }

  void Context::allocate_graph_data_(unsigned id)
  {
    if (graph_data_.size() < id + 1)
      graph_data_.resize(id+1);
    graph_data_[id].reset(new GraphData());
  }

  void Context::delete_graph_data_(unsigned id)
  {
    if (id < graph_data_.size())
      graph_data_[id].reset(new GraphData());
  }

  void Context::delete_node_data_(unsigned graphid, unsigned id)
  {
    if (graphid < graph_data_.size() && graph_data_[graphid].get())
      graph_data_[graphid]->delete_data_(id);
  }


  Graph * Session::graph(unsigned id)
  {
    if (!graph_exists(id))
      return nullptr;
    return graphs_[id].get();
  }

  const Graph * Session::graph(unsigned id) const
  {
    if (!graph_exists(id))
      return nullptr;
    return graphs_[id].get();
  }

  bool Session::graph_exists(unsigned id) const
  {
    return graphs_.size() > id && graphs_[id].get() != nullptr;
  }

  bool Session::graph_can_be_evaluated(unsigned id) const
  {
    const Graph * g = graph(id);
    if (!g)
      return false;
    return g->out_node_id() != ALG_NO_ID;
  }

  Node * Session::node(unsigned graphid, unsigned nodeid)
  {
    Graph * g = graph(graphid);
    if (!g)
      return nullptr;
    return g->node(nodeid);
  }

  const Node * Session::node(unsigned graphid, unsigned nodeid) const
  {
    const Graph * g = graph(graphid);
    if (!g)
      return nullptr;
    return g->node(nodeid);
  }

  const Algorithm * Session::algorithm(unsigned graphid, unsigned nodeid) const
  {
    const Node * n = node(graphid, nodeid);
    if (!n)
      return nullptr;
    return n->algorithm();
  }

  Algorithm * Session::algorithm(unsigned graphid, unsigned nodeid)
  {
    Node * n = node(graphid, nodeid);
    if (!n)
      return nullptr;
    return n->algorithm();
  }

  bool Session::node_exists(unsigned graphid, unsigned nodeid) const
  {
    const Graph * g = graph(graphid);
    if (!g)
      return false;
    return g->node_exists(nodeid);
  }

  AlgorithmData * Session::alg_data(unsigned graphid, unsigned nodeid)
  {
    if (!node_exists(graphid, nodeid))
      return nullptr;
    return ctxt_.graph_data(graphid)->algdata_[nodeid].get();
  }

  AlgorithmData * Session::subctxt_alg_data(unsigned subctxt,
                                            unsigned graphid, unsigned nodeid)
  {
    if (!node_exists(graphid, nodeid))
      return nullptr;
    if (subctxt >= sub_ctxt_.size())
      return nullptr;
    if (!sub_ctxt_[subctxt].has_graph_data(graphid))
      return nullptr;
    GraphData * gdata = sub_ctxt_[subctxt].graph_data_[graphid].get();
    if (nodeid >= gdata->algdata_.size())
      return nullptr;
    return gdata->algdata_[nodeid].get();
  }

  void Session::invalidate_subctxt_alg_data(unsigned graphid, unsigned nodeid)
  {
    for (auto & sctxt : sub_ctxt_)
      sctxt.delete_node_data_(graphid, nodeid);
  }

  void Session::invalidate_reconstruction_cache(unsigned graphid)
  {
    Graph * g = graph(graphid);
    if (g)
      g->invalidate_reconstruction_cache();
  }

  unsigned Session::new_graph()
  {
    unsigned id = get_free_graph_id_();
    graphs_[id].reset(new Graph());
    Graph & graph = *graphs_[id];

    graph.id_ = id;
    graph.session_ = this;
    ctxt_.allocate_graph_data_(id);

    return id;
  }

  void Session::delete_graph(unsigned graphid)
  {
    if (!graph_exists(graphid))
      return;

    graphs_[graphid]->session_ = nullptr;
    graphs_[graphid].reset();
    delete_graph_data_(graphid);

    free_graph_slots_.push_back(graphid);
  }

  unsigned Session::delete_node(unsigned graphid, unsigned nodeid)
  {
    Graph * g = graph(graphid);
    if (!g)
      return ALG_NO_ID;
    return g->delete_node(nodeid);
  }

  unsigned Session::set_output_node(unsigned graphid, unsigned nodeid)
  {
    if (!graph_exists(graphid))
      return ALG_NO_ID;
    return graphs_[graphid]->set_output_node(nodeid);
  }

  unsigned Session::set_node_mutable(unsigned graphid, unsigned nodeid)
  {
    if (!node_exists(graphid, nodeid))
      return ALG_NO_ID;

    Graph * g = graph(graphid);

    if (!g->is_mutable())
      return ALG_NO_ID;

    for (const auto & node : g->nodes_)
      if (node.get())
        for (auto ins : node->inputs_)
          if (ins == nodeid)
            return ALG_NO_ID;

    algorithm(graphid, nodeid)->set_mutable();

    return nodeid;
  }

  unsigned Session::prune_graph(unsigned graphid)
  {
    if (!graph_exists(graphid))
      return ALG_NO_ID;

    Graph * graph = graphs_[graphid].get();
    if (!graph->is_mutable())
      return ALG_NO_ID;

    for (auto & node : graph->nodes_) {
      if (!graph->mark_nodes_[node->id()])
        graph->delete_node_(node->id(), true);
    }

    // Now the graph is mutable and the output node is not used by
    // other nodes, so we might as well make it mutable
    graph->out_node()->algorithm()->set_mutable();

    return graphid;
  }

  void Session::make_reconstructible(unsigned graphid)
  {
    if (graph_exists(graphid))
      make_reconstructible_(graphid);
  }

  unsigned Session::get_free_graph_id_()
  {
    unsigned id=0;
    if (!free_graph_slots_.empty()) {

      id = free_graph_slots_.back();
      free_graph_slots_.pop_back();

    } else {

      id = graphs_.size();
      graphs_.push_back(GraphPtr(nullptr));

    }

    return id;
  }

  void Session::delete_node_data_(unsigned graphid, unsigned id)
  {
    ctxt_.delete_node_data_(graphid, id);
    for (auto & sctxt : sub_ctxt_)
      sctxt.delete_node_data_(graphid, id);
  }

  void Session::delete_graph_data_(unsigned graphid)
  {
    ctxt_.delete_graph_data_(graphid);
    for (auto & sctxt : sub_ctxt_)
      sctxt.delete_graph_data_(graphid);
  }

  AlgorithmData & Session::algdata_(unsigned graphid, unsigned nodeid)
  {
    return algdata_(ctxt_, graphid, nodeid);
  }

  AlgorithmData & Session::algdata_(Context & ctxt,
                                    unsigned graphid, unsigned nodeid)
  {
    return *ctxt.graph_data_[graphid]->algdata_[nodeid];
  }

  void Session::init_subcontexts_(unsigned n)
  {
    if (sub_ctxt_.size() < n)
      sub_ctxt_.resize(n);
  }

  void Session::make_reconstructible_(unsigned graphid)
  {
    if (!graphs_[graphid]->rec_data_.get()) {
      graphs_[graphid]->rec_data_.reset(new Graph::AlgRecData());
      set_default_shift_(graphid);
    }
  }

  void Session::set_default_shift_(unsigned id)
  {
    Mod mod(BIG_UINT_PRIMES[BIG_UINT_PRIMES_SIZE-1]);
    const unsigned nvars = graphs_[id]->nparsin[0];
    auto & shift = graphs_[id]->rec_data().shift;
    if (nvars)
      shift.reset(new UInt[nvars]);
    for (unsigned i=0; i<nvars; ++i)
      shift[i] = sample_uint(OFFSET_4, i, mod);
  }

  Ret Session::learn(unsigned graphid)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;

    Graph & graph = *graphs_[graphid];

    Node * node = graph.out_node();
    if (!node)
      return FAILED;
    const LearningOptions & opt = node->learn_opt;

    Ret ret = graph.learn_all(&ctxt_, Mod(BIG_UINT_PRIMES[opt.prime_no]),
                              ctxt_.graph_data(graphid),
                              opt.n_singular);
    node->alg_->set_learned(graph.has_learned());

    if (ret == SUCCESS) {
      auto * gdata = ctxt_.graph_data(graphid);
      gdata->set_nparsout_(graph.out_node_id(), node->alg_->nparsout);
    }

    return ret;
  }

  Ret Session::degrees(unsigned graphid, const ReconstructionOptions & opt)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return FAILED;

    make_reconstructible_(graphid);

    if (!a.has_learned())
      return FAILED;

    Graph::AlgRecData & arec = *a.rec_data_;;
    Graph::AlgDegs & degs = arec.degs;
    const unsigned nparsout = a.nparsout;

    degs.numdeg.reset(new unsigned[nparsout]);
    degs.dendeg.reset(new unsigned[nparsout]);

    return algorithm_get_degrees(a, ctxt_.graph_data(graphid), &ctxt_,
                                 arec.shift.get(), opt,
                                 degs.numdeg.get(), degs.dendeg.get());
  }

  Ret Session::var_degrees_(unsigned graphid, unsigned var,
                            const ReconstructionOptions & opt)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return FAILED;

    Graph::AlgRecData & arec = *a.rec_data_;;
    Graph::AlgDegs & degs = arec.degs;

    RatFunVarDegrees * vdegs = degs.vdegs.get();
    return algorithm_get_var_degrees(a, ctxt_.graph_data(graphid), &ctxt_,
                                     var, arec.shift.get(), opt,
                                     vdegs);
  }

  Ret Session::all_var_degrees(unsigned graphid,
                               const ReconstructionOptions & opt)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return FAILED;

    make_reconstructible_(graphid);

    if (!a.has_learned())
      return FAILED;

    const unsigned nparsin = a.nparsin[0];
    const unsigned nparsout = a.nparsout;
    Graph::AlgDegs & degs = a.rec_data_->degs;

    degs.vdegs.reset(new RatFunVarDegrees[nparsout]);
    for (unsigned i=0; i<nparsout; ++i)
      degs.vdegs[i].resize(nparsin);

    for (unsigned var=0; var<nparsin; ++var) {
      Ret ret = var_degrees_(graphid, var, opt);
      if (ret != SUCCESS)
        return ret;
    }
    return SUCCESS;
  }

  Ret Session::all_degrees(unsigned graphid, const ReconstructionOptions & opt)
  {
    Ret ret = degrees(graphid, opt);
    if (ret != SUCCESS)
      return ret;
    ret = all_var_degrees(graphid, opt);
    return ret;
  }

  void Session::share_with_subctxt_(unsigned graphid, Context & sctxt)
  {
    if (!graph_exists(graphid))
      return;

    const Graph & g = *graphs_[graphid];
    if (!sctxt.has_graph_data(graphid))
      sctxt.allocate_graph_data_(graphid);
    if (sctxt.graph_data_.size() < graphid+1)
      sctxt.graph_data_.resize(graphid+1);
    if (sctxt.xin_.size() < ctxt_.xin_.size())
      sctxt.xin_.resize(ctxt_.xin_.size());
    sctxt.ww.ensure_size(ctxt_.ww.size());
    g.sync_graph_data(*ctxt_.graph_data_[graphid].get(),
                      *sctxt.graph_data_[graphid].get());
  }

  void Session::share_with_subctxt_(unsigned graphid, unsigned i)
  {
    share_with_subctxt_(graphid, sub_ctxt_[i]);
  }

  void Session::alloc_threads_(unsigned n)
  {
#if FFLOW_THREAD_POOL
    pool_.alloc_threads(n);
#else
    (void)(n);
#endif
  }

  template<typename F, typename... Ts>
  auto Session::enqueue_(F&& f, Ts&&... args)
    -> std::future<typename std::result_of<F(Ts...)>::type>
  {
#if FFLOW_THREAD_POOL
    return pool_.enqueue(std::forward<F>(f), std::forward<Ts>(args)...);
#else
    return std::async(std::launch::async,
                      std::forward<F>(f), std::forward<Ts>(args)...);
#endif
  }


  struct GraphParallelTotDegs {

    Ret operator() ()
    {
      session->share_with_subctxt_(graphid, *sctxt);

      Graph & a = *session->graphs_[graphid];
      Graph::AlgRecData & arec = a.rec_data();
      Graph::AlgDegs & degs = arec.degs;

      return algorithm_get_degrees(a, sctxt->graph_data(graphid), sctxt,
                                   arec.shift.get(), *opt,
                                   degs.numdeg.get(), degs.dendeg.get());
    }

    Session * session;
    Context * sctxt;
    const ReconstructionOptions * opt;
    unsigned graphid;
  };

  struct GraphParallelVarDegs {

    Ret operator() (unsigned var_beg, unsigned step)
    {
      session->share_with_subctxt_(graphid, *sctxt);

      Graph & a = *session->graphs_[graphid];
      Graph::AlgRecData & arec = a.rec_data();
      Graph::AlgDegs & degs = arec.degs;

      const unsigned nparsin = a.nparsin[0];
      RatFunVarDegrees * vdegs = degs.vdegs.get();

      Ret ret = SUCCESS;
      for (unsigned v=var_beg; v<nparsin; v+=step) {
        ret = algorithm_get_var_degrees(a, sctxt->graph_data(graphid), sctxt,
                                        v, arec.shift.get(), *opt,
                                        vdegs);
        if (ret != SUCCESS)
          return ret;
      }

      return SUCCESS;
    }

    Session * session;
    Context * sctxt;
    const ReconstructionOptions * opt;
    unsigned graphid;
  };


  Ret Session::parallel_all_degrees(unsigned graphid, unsigned nthreads,
                                    const ReconstructionOptions & opt)
  {
    if (nthreads == 1)
      return all_degrees(graphid, opt);
    if (nthreads == 0)
      nthreads = std::thread::hardware_concurrency();

    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return FAILED;

    make_reconstructible_(graphid);

    Graph::AlgRecData & arec = a.rec_data();
    Graph::AlgDegs & degs = arec.degs;
    const unsigned nparsin = a.nparsin[0];
    const unsigned nparsout = a.nparsout;

    if (nthreads > nparsin+1)
      nthreads = nparsin+1;

    init_subcontexts_(nthreads);
    alloc_threads_(nthreads);

    std::future<Ret> deg_ret;
    auto td_eval = GraphParallelTotDegs{this, &sub_ctxt_[0], &opt, graphid};
    {
      degs.numdeg.reset(new unsigned[nparsout]);
      degs.dendeg.reset(new unsigned[nparsout]);
      deg_ret = enqueue_(td_eval);
    }

    std::vector<std::future<Ret>> vdeg_ret(nthreads-1);
    std::vector<GraphParallelVarDegs> t_eval(nthreads-1);
    {
      degs.vdegs.reset(new RatFunVarDegrees[nparsout]);
      for (unsigned i=0; i<nparsout; ++i)
        degs.vdegs[i].resize(nparsin);
      for (unsigned i=0; i<nthreads-1; ++i) {
        t_eval[i] = GraphParallelVarDegs{this, &sub_ctxt_[i+1], &opt, graphid};
        vdeg_ret[i] = enqueue_(t_eval[i], i, nthreads-1);
      }
    }

    Ret retv = SUCCESS;
    Ret reti = deg_ret.get();
    if (reti != SUCCESS && retv == SUCCESS)
      retv = reti;
    for (unsigned i=0; i<nthreads-1; ++i) {
      reti = vdeg_ret[i].get();
      if (reti != SUCCESS && retv == SUCCESS)
        retv = reti;
    }

    return retv;
  }

  void Session::init_cache_(CacheUPtr & cache,
                            unsigned nparsin, unsigned,
                            unsigned nsamples, unsigned allocsize)
  {
    cache.reset(new UIntCache());
    (*cache).init(nparsin+1);
    if (nsamples)
      (*cache).reserve_more(nsamples, allocsize);
  }

  unsigned Session::compute_cache_size_(unsigned nparsin, unsigned nparsout,
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

  void Session::init_graph_cache_(unsigned graphid)
  {
    Graph & g = *graphs_[graphid];
    init_cache_(g.rec_data().cache, g.nparsin[0], g.nparsout, 0, 0);
  }

  void Session::gen_sample_points_(unsigned graphid,
                                   SamplePointsVector & samples,
                                   const ReconstructionOptions & opt)
  {
    Graph & a = *graphs_[graphid];

    Graph::AlgRecData & arec = a.rec_data();
    Graph::AlgDegs & degs = arec.degs;

    GenerateSamplePoints sfun(a.nparsin[0]);
    sfun.set_complement(arec.cache.get());

    unsigned nparsin = a.nparsin[0];
    unsigned nparsout = a.nparsout;

    algorithm_generate_sample_points(sfun, nparsin, nparsout,
                                     arec.shift.get(), opt,
                                     degs.vdegs.get(),
                                     degs.numdeg.get(), degs.dendeg.get());
    unsigned flags_size = bit_array_u64size(nparsout);
    sfun.append_to_vector(samples, flags_size);
    sort_by_mod(samples.data(), samples.data()+samples.size(), nparsin);

    std::vector<unsigned> to_rec;
    set_up_to_rec_(graphid, to_rec);
    algorithm_verify_sample_points(samples.data(), samples.size(),
                                   nparsin, nparsout,
                                   to_rec,
                                   arec.shift.get(), opt,
                                   degs.vdegs.get(),
                                   degs.numdeg.get(), degs.dendeg.get());
    auto nend = std::remove_if(samples.begin(), samples.end(),
                               [nparsin,nparsout]
                               (const std::unique_ptr<UInt[]> & el)
                               {
                                 ConstBitArrayView av(el.get()+nparsin+1);
                                 return bit_array_nonzeroes(av, nparsout) == 0;
                               });
    samples.erase(nend, samples.end());
    std::vector<unsigned> ns(nparsout);
    for (auto & s : samples) {
      for (unsigned j=0; j<nparsout; ++j)
        if (ConstBitArrayView(s.get()+nparsin+1).get(j))
          ns[j] += 1;
    }

    if (opt.dbginfo) {
      dbgprint(format("Generated {} sample points", samples.size()));
      std::unique_ptr<UInt[]> xi(new UInt[a.nparsin[0]]);
      std::unique_ptr<UInt[]> xo(new UInt[a.nparsout]);
      Mod mod(BIG_UINT_PRIMES[0]);
      {
        // do a warm-up evaluation first, since first evaluation is
        // often slower, so we only time the second one
        for (unsigned i=0; i<a.nparsin[0]; ++i)
          xi[i] = sample_uint(OFFSET_1, 40 + i, mod);
        const UInt * xiptr = xi.get();
        a.evaluate(&ctxt_, &xiptr, mod, ctxt_.graph_data(graphid), xo.get());
      }
      for (unsigned i=0; i<a.nparsin[0]; ++i)
        xi[i] = sample_uint(OFFSET_1, 20 + i, mod);
      const UInt * xiptr = xi.get();
      auto start = std::chrono::system_clock::now();
      a.evaluate(&ctxt_, &xiptr, mod, ctxt_.graph_data(graphid), xo.get());
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end-start;
      dbgprint(format("Approximate time per evaluation (single core): {} sec.",
                      elapsed_seconds.count()));
    }
  }

  void Session::sample(unsigned graphid, const ReconstructionOptions & opt,
                       SamplePointsGenerator * samplegen)
  {
    if (!graph_can_be_evaluated(graphid))
      return;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return;

    if (!a.rec_data_.get()) {
      if (a.nparsin[0] == 1)
        make_reconstructible_(a.id_);
      else
        return;
    }

    SamplePointsVector samples;
    if (samplegen == nullptr) {
      gen_sample_points_(graphid, samples, opt);
    } else {
      Ret dret = samplegen->load_samples(a.nparsin[0], a.nparsout, samples);
      if (dret != SUCCESS)
        return;
    }

    Graph::AlgRecData & arec = a.rec_data();

    if (arec.cache.get() == nullptr)
      init_graph_cache_(graphid);
    UIntCache & cache = *arec.cache;

    {
      const unsigned n = samples.size();
      const unsigned alloc_size = compute_cache_size_(a.nparsin[0], a.nparsout,
                                                      samples.data(),
                                                      samples.data()+n);
      cache.reserve_more(n, alloc_size);
    }

    evaluate_and_sparse_cache_samples(a, ctxt_.graph_data(graphid), &ctxt_,
                                      samples.data(),
                                      samples.data() + samples.size(),
                                      cache);
  }

  struct GraphParallelEvaluate {

    void operator() (CacheUPtr & tcache)
    {
      session->share_with_subctxt_(graphid, *ctxt);
      Algorithm & alg = *session->graphs_[graphid];
      evaluate_and_sparse_cache_samples(alg, ctxt->graph_data(graphid), ctxt,
                                        start, end, *tcache);
    }

    Session * session;
    Context * ctxt;
    unsigned graphid;
    const std::unique_ptr<UInt[]> *start, *end;
  };


  void Session::parallel_sample(unsigned graphid, unsigned nthreads,
                                const ReconstructionOptions & opt,
                                SamplePointsGenerator * samplegen)
  {
    if (nthreads == 1) {
      sample(graphid, opt, samplegen);
      return;
    }

    if (nthreads == 0)
      nthreads = std::thread::hardware_concurrency();

    if (!graph_can_be_evaluated(graphid))
      return;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return;

    if (!a.rec_data_.get()) {
      if (a.nparsin[0] == 1)
        make_reconstructible_(a.id_);
      else
        return;
    }

    Graph::AlgRecData & arec = a.rec_data();

    SamplePointsVector samples;
    if (samplegen == nullptr) {
      gen_sample_points_(graphid, samples, opt);
    } else {
      Ret dret = samplegen->load_samples(a.nparsin[0], a.nparsout, samples);
      if (dret != SUCCESS)
        return;
    }

    const unsigned tot_samples = samples.size();
    std::unique_ptr<UInt[]> * start = samples.data();

    nthreads = std::min(nthreads, tot_samples);

    std::vector<CacheUPtr> caches(nthreads);
    std::vector<std::future<void>> ret(nthreads);
    std::vector<GraphParallelEvaluate> t_eval(nthreads);

    init_subcontexts_(nthreads);

    alloc_threads_(nthreads);
    const std::unique_ptr<UInt[]> * istart = start;
    for (unsigned i=0; i<nthreads; ++i) {
      const unsigned n = tot_samples/nthreads + (i < (tot_samples % nthreads));
      const unsigned alloc_size = compute_cache_size_(a.nparsin[0], a.nparsout,
                                                      istart, istart+n);
      init_cache_(caches[i], a.nparsin[0], a.nparsout, n, alloc_size);
      istart += n;
    }

    for (unsigned i=0; i<nthreads; ++i) {
      const unsigned n = tot_samples/nthreads + (i < (tot_samples % nthreads));
      t_eval[i] = GraphParallelEvaluate{this, &sub_ctxt_[i], graphid,
                                        start, start+n};
      ret[i] = enqueue_(t_eval[i], std::ref(caches[i]));
      start += n;
    }

    for (unsigned i=0; i<nthreads; ++i)
      ret[i].get();

    if (arec.cache.get() == nullptr)
      init_graph_cache_(graphid);

    UIntCache  & cache = *arec.cache;
    merge_function_caches(caches.data(), nthreads, cache);
  }

  Ret Session::dump_sample_points(unsigned graphid,
                                  const ReconstructionOptions & opt,
                                  const char * filename)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return FAILED;

    if (!a.rec_data_.get()) {
      if (a.nparsin[0] == 1)
        make_reconstructible_(a.id_);
      else
        return FAILED;
    }

    SamplePointsVector samples;
    gen_sample_points_(graphid, samples, opt);

    return dump_samples(filename,
                        a.nparsin[0] + bit_array_u64size(a.nparsout),
                        samples);
  }


  Ret Session::evaluate_list(unsigned graphid,
                             const SamplePointsVector & input,
                             SamplePointsVector & res,
                             unsigned nthreads)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0] || a.nparsout<1)
      return FAILED;

    if (!a.rec_data_.get())
      make_reconstructible_(a.id_);

    unsigned nparsin = a.nparsin[0];
    unsigned nparsout = a.nparsout;
    unsigned npoints = input.size();
    unsigned flagsize = bit_array_u64size(nparsout);
    unsigned xsize = nparsin + 1 + flagsize;
    std::unique_ptr<UInt[]> x;

    res.resize(0);
    res.reserve(input.size());

    SamplePointsVector points;

    for (unsigned j=0; j<npoints; ++j) {
      x.reset(new UInt[xsize]);
      std::copy(input[j].get(), input[j].get()+nparsin+1, x.get());
      std::fill(x.get()+nparsin+1, x.get()+xsize, ~UInt(0));
      points.push_back(std::move(x));
    }
    sort_by_mod(points.data(), points.data()+points.size(), nparsin);

    SamplePointsFromVector vec;
    vec.setVector(std::move(points));

    ReconstructionOptions opt;
    parallel_sample(graphid, nthreads, opt, &vec);

    Graph::AlgRecData & arec = a.rec_data();
    UIntCache & cache = *arec.cache;

    res.clear();
    res.reserve(input.size());

    for (unsigned j=0; j<npoints; ++j) {
      x.reset(new UInt[nparsout]);
      UInt * xout = nullptr;
      if (cache.find(input[j].get(), &xout))
        std::copy(xout+flagsize, xout+flagsize+nparsout, x.get());
      else
        x[0] = FAILED;
      res.push_back(std::move(x));
    }

    return SUCCESS;
  }


  Ret Session::dump_evaluations(unsigned graphid,
                                const char * filename)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return FAILED;

    if (!a.rec_data_.get())
      return FAILED;

    Graph::AlgRecData & arec = a.rec_data();

    if (arec.cache.get() == nullptr)
      return FAILED;
    const UIntCache & cache = *arec.cache;

    const unsigned flag_size = bit_array_u64size(a.nparsout);
    return ::fflow::dump_evaluations(filename,
                                     a.nparsin[0], a.nparsout + flag_size,
                                     cache);
  }

  Ret Session::load_evaluations(unsigned graphid,
                                const char ** filename,
                                unsigned n_files)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return FAILED;

    if (!a.rec_data_.get()) {
      if (a.nparsin[0] == 1)
        make_reconstructible_(a.id_);
      else
        return FAILED;
    }

    Graph::AlgRecData & arec = a.rec_data();

    if (arec.cache.get() == nullptr)
      init_graph_cache_(graphid);
    UIntCache & cache = *arec.cache;

    const unsigned flag_size = bit_array_u64size(a.nparsout);

    for (unsigned j=0; j<n_files; ++j) {
      Ret ret = ::fflow::load_evaluations(filename[j],
                                          a.nparsin[0], a.nparsout + flag_size,
                                          cache);
      if (ret != SUCCESS)
        return ret;
    }

    return SUCCESS;
  }

  unsigned Session::new_dummy_graph(unsigned nparsin, unsigned nparsout)
  {
    unsigned graphid = new_graph();
    Graph & a = *graphs_[graphid];
    unsigned in = a.set_input_vars(nparsin);

    std::unique_ptr<NoAlgorithm> algptr(new NoAlgorithm());
    algptr->nparsin.resize(1);
    algptr->nparsin[0] = nparsin;
    algptr->nparsout = nparsout;

    unsigned id = a.new_node(std::move(algptr), nullptr, &in);
    a.set_output_node(id);

    return graphid;
  }

  Ret Session::dump_degrees(unsigned graphid, const char * file) const
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    const Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return FAILED;

    if (!graphs_[graphid]->rec_data_.get())
      return FAILED;

    Graph::AlgRecData & arec = *a.rec_data_;;
    Graph::AlgDegs & degs = arec.degs;
    RatFunVarDegrees * vdegs = degs.vdegs.get();

    if (!degs.numdeg.get() || !degs.dendeg.get() || !vdegs)
      return FAILED;

    return algorithm_dump_degree_info(file, a.nparsin[0], a.nparsout,
                                      degs.numdeg.get(),
                                      degs.dendeg.get(),
                                      vdegs);
  }

  Ret Session::load_degrees(unsigned graphid, const char * file)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (!a.nparsin[0])
      return FAILED;

    make_reconstructible_(graphid);

    if (!a.has_learned())
      return FAILED;

    const unsigned nparsin = a.nparsin[0];
    const unsigned nparsout = a.nparsout;
    Graph::AlgDegs & degs = a.rec_data_->degs;

    degs.numdeg.reset(new unsigned[nparsout]);
    degs.dendeg.reset(new unsigned[nparsout]);

    degs.vdegs.reset(new RatFunVarDegrees[nparsout]);
    for (unsigned i=0; i<nparsout; ++i)
      degs.vdegs[i].resize(nparsin);

    return algorithm_load_degree_info(file, a.nparsin[0], a.nparsout,
                                      degs.numdeg.get(),
                                      degs.dendeg.get(),
                                      degs.vdegs.get());
  }

  Ret Session::reconstruct_(unsigned graphid,
                            unsigned var_beg, unsigned step,
                            const std::vector<unsigned> & to_rec,
                            const ReconstructionOptions & opt)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];
    if (!a.nparsin[0])
      return FAILED;
    if (!a.rec_data_.get())
      return MISSING_SAMPLES;
    Graph::AlgRecData & arec = a.rec_data();
    Graph::AlgDegs & degs = arec.degs;
    UIntCache  & cache = *arec.cache;
    MPReconstructedRatFun * res = arec.recratfun_.get();

    unsigned nparsin = a.nparsin[0], nparsout = a.nparsout;

    Ret final_ret = SUCCESS;
    for (unsigned i=var_beg; i<to_rec.size(); i+=step) {
      Ret ret = algorithm_sparse_reconstruct(cache, nparsin, nparsout,
                                             to_rec[i],
                                             arec.shift.get(), opt,
                                             degs.vdegs.get(),
                                             degs.numdeg.get(),
                                             degs.dendeg.get(),
                                             res[to_rec[i]]);
      if (ret != SUCCESS) {
        if (final_ret == SUCCESS)
          final_ret = ret;
      } else {
        arec.completed_[to_rec[i]] = true;
      }
    }

    return final_ret;
  }

  Ret Session::reconstruct_mod_(unsigned graphid,
                                unsigned var_beg, unsigned step,
                                const std::vector<unsigned> & to_rec,
                                const ReconstructionOptions & opt,
                                Mod mod)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];
    if (!a.nparsin[0])
      return FAILED;
    if (!a.rec_data_.get())
      return MISSING_SAMPLES;
    Graph::AlgRecData & arec = a.rec_data();
    Graph::AlgDegs & degs = arec.degs;
    UIntCache  & cache = *arec.cache;
    MPReconstructedRatFun * res = arec.recratfun_.get();

    unsigned nparsin = a.nparsin[0], nparsout = a.nparsout;

    Ret final_ret = SUCCESS;
    for (unsigned i=var_beg; i<to_rec.size(); i+=step) {
      Ret ret = algorithm_sparse_reconstruct_mod(mod,
                                                 cache, nparsin, nparsout,
                                                 to_rec[i],
                                                 arec.shift.get(), opt,
                                                 degs.vdegs.get(),
                                                 degs.numdeg.get(),
                                                 degs.dendeg.get(),
                                                 res[to_rec[i]]);
      if (ret != SUCCESS) {
        if (final_ret == SUCCESS)
          final_ret = ret;
      }
    }

    return final_ret;
  }

  void Session::set_up_to_rec_(unsigned id, std::vector<unsigned> & to_rec)
  {
    if (!graph_can_be_evaluated(id))
      return;
    Graph & a = *graphs_[id];
    if (!a.rec_data_.get())
      return;
    Graph::AlgRecData & arec = a.rec_data();

    if (!arec.completed_.get()) {
      arec.recratfun_.reset(new MPReconstructedRatFun[a.nparsout]);
      arec.completed_.reset(new bool[a.nparsout]());
      to_rec.resize(a.nparsout);
      std::iota(to_rec.begin(), to_rec.end(), 0);
    } else {
      to_rec.resize(0);
      for (unsigned j=0; j<a.nparsout; ++j)
        if (!arec.completed_[j])
          to_rec.push_back(j);
    }
  }

  void Session::move_rec_(unsigned id, MPReconstructedRatFun res[])
  {
    Graph & a = *graphs_[id];
    Graph::AlgRecData & arec = a.rec_data();
    MPReconstructedRatFun * recf = arec.recratfun_.get();
    for (unsigned j=0; j<a.nparsout; ++j)
      res[j].copy(std::move(recf[j]));
    arec.completed_.reset();
    arec.recratfun_.reset();
  }

  Ret Session::reconstruct(unsigned id, MPReconstructedRatFun res[],
                           const ReconstructionOptions & opt)
  {
    std::vector<unsigned> to_rec;
    set_up_to_rec_(id, to_rec);
    Ret ret = reconstruct_(id, 0, 1, to_rec, opt);
    if (ret == SUCCESS)
      move_rec_(id, res);
    return ret;
  }

  Ret Session::reconstruct_mod(unsigned id, MPReconstructedRatFun res[],
                               const ReconstructionOptions & opt)
  {
    Mod mod(prime_no(opt.start_mod));
    std::vector<unsigned> to_rec;
    set_up_to_rec_(id, to_rec);
    Ret ret = reconstruct_mod_(id, 0, 1, to_rec, opt, mod);
    if (ret == SUCCESS)
      move_rec_(id, res);
    return ret;
  }

  struct GraphParallelReconstruct {

    Ret operator() (unsigned var_beg, unsigned step,
                    const std::vector<unsigned> & to_rec)
    {
      return session->reconstruct_(graphid, var_beg, step, to_rec, *opt);
    }

    Session * session;
    const ReconstructionOptions * opt;
    unsigned graphid;
  };

  struct GraphParallelReconstructMod {

    Ret operator() (unsigned var_beg, unsigned step,
                    const std::vector<unsigned> & to_rec)
    {
      return session->reconstruct_mod_(graphid, var_beg, step,
                                       to_rec, *opt, mod);
    }

    Mod mod;
    Session * session;
    const ReconstructionOptions * opt;
    unsigned graphid;
  };

  Ret Session::parallel_reconstruct(unsigned graphid,
                                    MPReconstructedRatFun res[],
                                    unsigned nthreads,
                                    const ReconstructionOptions & opt)
  {
    if (nthreads == 1)
      return reconstruct(graphid, res, opt);

    if (nthreads == 0)
      nthreads = std::thread::hardware_concurrency();

    if (!graph_can_be_evaluated(graphid))
      return FAILED;

    std::vector<unsigned> to_rec;
    set_up_to_rec_(graphid, to_rec);

    nthreads = std::min(nthreads, unsigned(to_rec.size()));
    alloc_threads_(nthreads);

    std::vector<GraphParallelReconstruct> t_eval(nthreads);
    std::vector<std::future<Ret>> ret(nthreads);
    for (unsigned i=0; i<nthreads; ++i)
      t_eval[i] = GraphParallelReconstruct{this, &opt, graphid};

    for (unsigned i=0; i<nthreads; ++i)
      ret[i] = enqueue_(t_eval[i], i, nthreads, to_rec);

    Ret retv = SUCCESS;
    for (unsigned i=0; i<nthreads; ++i) {
      Ret reti = ret[i].get();
      if (reti != SUCCESS && retv == SUCCESS)
        retv = reti;
    }

    if (retv == SUCCESS)
      move_rec_(graphid, res);

    return retv;
  }

  Ret Session::parallel_reconstruct_mod(unsigned graphid,
                                        MPReconstructedRatFun res[],
                                        unsigned nthreads,
                                        const ReconstructionOptions & opt)
  {
    if (nthreads == 1)
      return reconstruct_mod(graphid, res, opt);

    Mod mod(prime_no(opt.start_mod));

    if (nthreads == 0)
      nthreads = std::thread::hardware_concurrency();

    if (!graph_can_be_evaluated(graphid))
      return FAILED;

    std::vector<unsigned> to_rec;
    set_up_to_rec_(graphid, to_rec);

    nthreads = std::min(nthreads, unsigned(to_rec.size()));
    alloc_threads_(nthreads);

    std::vector<GraphParallelReconstructMod> t_eval(nthreads);
    std::vector<std::future<Ret>> ret(nthreads);
    for (unsigned i=0; i<nthreads; ++i)
      t_eval[i] = GraphParallelReconstructMod{mod, this, &opt, graphid};

    for (unsigned i=0; i<nthreads; ++i)
      ret[i] = enqueue_(t_eval[i], i, nthreads, to_rec);

    Ret retv = SUCCESS;
    for (unsigned i=0; i<nthreads; ++i) {
      Ret reti = ret[i].get();
      if (reti != SUCCESS && retv == SUCCESS)
        retv = reti;
    }

    if (retv == SUCCESS)
      move_rec_(graphid, res);

    return retv;
  }

  Ret Session::reconstruct_univariate(unsigned graphid,
                                      MPReconstructedRatFun res[],
                                      const ReconstructionOptions & opt)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];
    if (a.nparsin[0] != 1)
      return FAILED;
    make_reconstructible_(graphid);
    Graph::AlgRecData & arec = a.rec_data();

    return algorithm_reconstruct_univariate(a,
                                            ctxt_.graph_data(graphid), &ctxt_,
                                            arec.shift.get(), opt, res);
  }

  Ret Session::reconstruct_univariate_mod(unsigned graphid,
                                          MPReconstructedRatFun res[],
                                          const ReconstructionOptions & opt)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];
    if (a.nparsin[0] != 1)
      return FAILED;
    make_reconstructible_(graphid);
    Graph::AlgRecData & arec = a.rec_data();

    Mod mod(prime_no(opt.start_mod));
    return algorithm_reconstruct_univariate_mod(mod, a,
                                                ctxt_.graph_data(graphid),
                                                &ctxt_,
                                                arec.shift.get(), opt, res);
  }

  Ret Session::reconstruct_numeric(unsigned graphid, MPRational res[],
                                   const ReconstructionOptions & opt)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];
    if (a.nparsin[0] != 0)
      return FAILED;
    make_reconstructible_(graphid);
    return algorithm_reconstruct_numeric(a, ctxt_.graph_data(graphid), &ctxt_,
                                         opt, res);
  }

  Ret Session::full_reconstruction(unsigned graphid,
                                   MPReconstructedRatFun res[],
                                   unsigned nthreads,
                                   ReconstructionOptions opt,
                                   unsigned min_primes)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (a.nparsin[0] < 1)
      return FAILED;

    if (a.nparsin[0] == 1)
      return reconstruct_univariate(graphid, res, opt);

    const unsigned max_primes = opt.max_primes;

    Ret ret = parallel_all_degrees(graphid, nthreads, opt);
    if (ret != SUCCESS)
      return ret;
    ret = MISSING_PRIMES;

    while (ret == MISSING_PRIMES && min_primes<=max_primes) {
      opt.max_primes = min_primes;

      parallel_sample(graphid, nthreads, opt);

      ret = parallel_reconstruct(graphid, res, nthreads, opt);
      if (ret != SUCCESS && ret != MISSING_PRIMES)
        return ret;

      ++min_primes;
    }

    return ret;
  }

  Ret Session::full_reconstruction_mod(unsigned graphid,
                                       MPReconstructedRatFun res[],
                                       unsigned nthreads,
                                       ReconstructionOptions opt)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (a.nparsin[0] < 1)
      return FAILED;

    if (a.nparsin[0] == 1)
      return reconstruct_univariate(graphid, res, opt);

    Ret ret = parallel_all_degrees(graphid, nthreads, opt);
    if (ret != SUCCESS)
      return ret;

    {
      opt.max_primes = 1;

      parallel_sample(graphid, nthreads, opt);

      ret = parallel_reconstruct_mod(graphid, res, nthreads, opt);
      if (ret != SUCCESS && ret != MISSING_PRIMES)
        return ret;
    }

    return ret;
  }

  UInt Session::independent_of_var(unsigned graphid, unsigned var)
  {
    Graph * a = graph(graphid);
    if (!a)
      return FAILED;

    unsigned nparsin = a->nparsin[0], nparsout = a->nparsout;

    if (var >= nparsin)
      return FAILED;

    std::unique_ptr<UInt> xdata(new UInt[nparsin+2*nparsout]);

    UInt * xin = xdata.get();
    UInt * xout1 = xdata.get() + nparsin;
    UInt * xout2 = xdata.get() + nparsin + nparsout;

    Node * n = a->out_node();
    if (!n)
      return FAILED;

    unsigned n_singular = n->learn_opt.n_singular;

    Mod mod(BIG_UINT_PRIMES[0]);

    const unsigned NEEDED_CHECKS = 2;
    unsigned iii = 0;
    unsigned checks = 0;
    unsigned fails = 0;

    while (fails <= n_singular) {

      for (unsigned i=0; i<nparsin; ++i)
        if (i != var)
          xin[i] = sample_uint(OFFSET_1, iii++, mod);
      xin[var] = sample_uint(OFFSET_2, 2*checks, mod);

      Ret ret = a->evaluate(&ctxt_, &xin, mod, ctxt_.graph_data(graphid),
                            xout1);
      if (ret == FAILED) {
        ++fails;
        continue;
      };

      xin[var] = sample_uint(OFFSET_2, 2*checks+1, mod);
      ret = a->evaluate(&ctxt_, &xin, mod, ctxt_.graph_data(graphid),
                        xout2);
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

  UInt Session::count_sample_points(unsigned graphid,
                                    const ReconstructionOptions & opt,
                                    unsigned count[],
                                    unsigned nthreads,
                                    const char * file)
  {
    if (!graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *graphs_[graphid];

    if (a.nparsin.size() == 0)
      return FAILED;

    unsigned nparsin = a.nparsin[0];
    unsigned nparsout = a.nparsout;

    if (nparsin == 0 || nparsin == 1)
      return FAILED;

    if (file == nullptr) {
      Ret ret = parallel_all_degrees(graphid, nthreads, opt);
      if (ret != SUCCESS)
        return FAILED;
    }

    SamplePointsVector samples;
    if (file == nullptr) {
      gen_sample_points_(graphid, samples, opt);
    } else {
      Ret dret = load_samples(file,
                              nparsin + bit_array_u64size(a.nparsout),
                              samples);
      if (dret != SUCCESS)
        return FAILED;
    }

    UInt nsamples = samples.size();

    if (count) {
      std::fill(count, count+nparsout, 0);
      for (unsigned j=0; j<nsamples; ++j) {
        ConstBitArrayView av(samples[j].get()+nparsin+1);
        for (unsigned k=0; k<nparsout; ++k)
          if (av.get(k))
            ++count[k];
      }
    }

    return nsamples;
  }

  Session::~Session()
  {
    graphs_.clear();
  }

  unsigned Session::default_nthreads()
  {
    return std::thread::hardware_concurrency();
  }


} // namespace fflow
