#ifndef FFLOW_CACHED_SUBGRAPH_HH
#define FFLOW_CACHED_SUBGRAPH_HH

#include <fflow/subgraph.hh>
#include <fflow/function_cache.hh>

namespace fflow {


  class CachedSubGraphData : public SubGraphData {
  public:
    const UIntCache * cache() const
    {
      return cache_.get();
    }

    std::unique_ptr<UIntCache> && move_cache()
    {
      return std::move(cache_);
    }
  private:
    friend class CachedSubGraph;
    std::unique_ptr<UInt[]> xin_ = nullptr;
    std::unique_ptr<UIntCache> cache_ = nullptr;
  };


  class CachedSubGraph : public SubGraph {
  public:

    Ret init(Session & session, unsigned graphid,
             CachedSubGraphData & data,
             const unsigned npars[], unsigned npars_size);

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

    virtual Ret evaluateFromGraph(Context * ctxt,
                                  const MutAlgInput xin[], Mod mod,
                                  const Graph * graph, SubGraphData * data,
                                  UInt * xout) const override;

    void merge_caches(std::unique_ptr<UIntCache> caches[],
                      std::size_t n);

    void set_default_subcache_size(unsigned n);

    const UIntCache & main_cache() const
    {
      return cache_;
    }

  public:
    UIntCache cache_;
    unsigned default_subcache_size_ = 0;
  };



  class CachedFromSubGraph;
  class CachedFromSubGraphData;

  class CachedFromSubGraphData : public AlgorithmData {
  private:
    friend class CachedFromSubGraph;
  private:
    std::unique_ptr<UInt[]> xin_ = nullptr;
  };

  class CachedFromSubGraph : public Algorithm {
  public:

    Ret init(const Session * session, unsigned graphid, unsigned nodeid);

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

  private:
    const Session * session_;
    unsigned graphid_, nodeid_;
  };

} // namespace fflow


#endif // FFLOW_CACHED_SUBGRAPH_HH
