#include <fflow/cached_subgraph.hh>

namespace fflow {

  Ret CachedSubGraph::init(Session & session, unsigned graphid,
                           CachedSubGraphData & data,
                           const unsigned npars[], unsigned npars_size)
  {
    Ret ret = SubGraph::init_subgraph(session, graphid, data,
                                      npars, npars_size);
    if (ret != SUCCESS)
      return ret;

    const Graph * graph = subgraph();
    const unsigned nin = graph->nparsin.size() ? graph->nparsin[0] : 0;
    cache_.init(nin+1);

    return SUCCESS;
  }


  AlgorithmData::Ptr
  CachedSubGraph::clone_data(const AlgorithmData * datain) const
  {
    auto * data = static_cast<const CachedSubGraphData*>(datain);
    std::unique_ptr<CachedSubGraphData> ptr(new CachedSubGraphData());
    SubGraph::copy_data_into(*data, *ptr);
    return std::move(ptr);
  }


  Ret CachedSubGraph::evaluateFromGraph(Context * ctxt,
                                        const MutAlgInput * xin, Mod mod,
                                        const Graph * graph,
                                        SubGraphData * datain,
                                        UInt * xout) const
  {
    auto & data = *static_cast<CachedSubGraphData*>(datain);
    const unsigned nin = graph->nparsin.size() ? graph->nparsin[0] : 0;
    const unsigned nout = nparsout;

    copy_input_params(xin[0], nin, mod, data.xin_);

    UInt * xo;
    bool lookup = cache_.find(data.xin_.get(), &xo);

    if (!lookup) {

      UIntCache * local_cache = data.cache_.get();

      if (local_cache == nullptr) {

        local_cache = new UIntCache();
        local_cache->init(nin+1);
        if (default_subcache_size_)
          local_cache->reserve_more(default_subcache_size_,
                                    default_subcache_size_*(nin + nout + 1));
        data.cache_.reset(local_cache);

      } else {

        lookup = local_cache->find(data.xin_.get(), &xo);

      }

    }

    if (lookup) {

      std::copy(xo, xo+nout, xout);

    } else {

      Ret ret = subgraph()->evaluate(ctxt, xin, mod, data.data(), xout);
      if (ret != SUCCESS)
        return ret;

      UIntCache * local_cache = data.cache_.get();
      UInt * xi = local_cache->get_new_inptr();
      std::copy(data.xin_.get(), data.xin_.get()+nin+1, xi);
      xo = local_cache->new_entry(xi, nout);
      std::copy(xout, xout+nout, xo);

    }

    return SUCCESS;
  }

  void CachedSubGraph::merge_caches(std::unique_ptr<UIntCache> caches[],
                                    std::size_t n)
  {
    cache_.merge_caches(caches, n);
  }

  void CachedSubGraph::set_default_subcache_size(unsigned n)
  {
    default_subcache_size_ = n;
  }



  Ret CachedFromSubGraph::init(const Session * session,
                               unsigned graphid, unsigned nodeid)
  {
    session_ = session;
    graphid_ = graphid;
    nodeid_ = nodeid;

    const Algorithm * alg = session->algorithm(graphid, nodeid);
    if (!alg)
      return FAILED;

    const CachedSubGraph * subg = dynamic_cast<const CachedSubGraph *>(alg);
    if (!subg)
      return FAILED;

    nparsin.resize(1);
    nparsin[0] = subg->nparsin[0];
    nparsout = subg->nparsout;

    return SUCCESS;
  }

  Ret CachedFromSubGraph::evaluate(Context *,
                                   AlgInput xin[], Mod mod,
                                   AlgorithmData * datain,
                                   UInt xout[]) const
  {
    const Algorithm * alg = session_->algorithm(graphid_, nodeid_);
    if (!alg)
      return FAILED;

    unsigned nin = nparsin[0];
    const CachedSubGraph * subg = dynamic_cast<const CachedSubGraph *>(alg);
    if (!subg || subg->nparsin[0] != nin || subg->nparsout != nparsout)
      return FAILED;

    auto & data = *static_cast<CachedFromSubGraphData*>(datain);
    const UIntCache & cache = subg->main_cache();

    copy_input_params(xin[0], nin, mod, data.xin_);

    UInt * xo;
    bool lookup = cache.find(data.xin_.get(), &xo);

    if (!lookup)
      return FAILED;

    std::copy(xo, xo+nparsout, xout);

    return SUCCESS;
  }

  AlgorithmData::Ptr
  CachedFromSubGraph::clone_data(const AlgorithmData *) const
  {
    std::unique_ptr<CachedFromSubGraphData> ptr(new CachedFromSubGraphData());
    return std::move(ptr);
  }

} // namespace fflow
