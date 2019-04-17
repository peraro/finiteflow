#ifndef FFLOW_SUBGRAPH_HH
#define FFLOW_SUBGRAPH_HH

#include <fflow/graph.hh>

namespace fflow {

  class SubGraphData : public AlgorithmData {
  public:

    GraphData * data()
    {
      return static_cast<GraphData*>(data_.get());
    }

    const GraphData * data() const
    {
      return static_cast<const GraphData*>(data_.get());
    }

  private:
    friend class SubGraph;
    std::unique_ptr<AlgorithmData> data_;
    std::unique_ptr<MutAlgInput[]> xin_;
    std::unique_ptr<UInt[]> xinval_;
  };


  class SubGraph : public Algorithm {
  public:

    virtual Ret learn(Context *, Input *, Mod, AlgorithmData *) override;

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

    void copy_data_into(const SubGraphData & datain,
                        SubGraphData & dataout) const;

    virtual Ret learnFromGraph(Context * ctxt,
                               const MutAlgInput xin[],  Mod mod,
                               const Graph * graph, SubGraphData * data);

    virtual Ret evaluateFromGraph(Context * ctxt,
                                  const MutAlgInput xin[], Mod mod,
                                  const Graph * graph, SubGraphData * data,
                                  UInt * xout) const = 0;

    const Graph * subgraph() const
    {
      return graph_.get();
    }

    virtual ~SubGraph();

  protected:

    Ret init_subgraph(Session & session, unsigned graphid,
                      SubGraphData & data,
                      const unsigned * param_nodes,
                      unsigned n_params_nodes,
                      unsigned n_aux_vars = 0);


  private:

    void chain_xin_(AlgInput xin[], MutAlgInput subxin[]) const;

    void set_xin_(SubGraphData & data) const;

    void dec_subgraph_count_();

  private:
    GraphPtr graph_;
    unsigned aux_vars_ = 0;
  };


  class SimpleSubGraph : public SubGraph {
  public:

    Ret init(Session & session, unsigned graphid,
             SubGraphData & data,
             const unsigned npars[], unsigned npars_size);

    virtual Ret evaluateFromGraph(Context * ctxt,
                                  const MutAlgInput xin[], Mod mod,
                                  const Graph * graph, SubGraphData * data,
                                  UInt * xout) const;

  };


  class MemoizedSubGraphData : public SubGraphData {
  private:
    friend class MemoizedSubGraph;
    std::unique_ptr<UInt[]> xin_;
    std::unique_ptr<UInt[]> xout_;
    UInt mod_ = 0;
  };

  class MemoizedSubGraph : public SubGraph {
  public:

    Ret init(Session & session, unsigned graphid,
             MemoizedSubGraphData & data,
             const unsigned npars[], unsigned npars_size);

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

    virtual Ret evaluateFromGraph(Context * ctxt,
                                  const MutAlgInput xin[], Mod mod,
                                  const Graph * graph, SubGraphData * data,
                                  UInt * xout) const override;

  };


  class SubGraphMap : public SubGraph {
  public:

    Ret init(Session & session, unsigned graphid,
             SubGraphData & data,
             const unsigned npars[], unsigned npars_size);

    virtual Ret evaluateFromGraph(Context * ctxt,
                                  const MutAlgInput xin[], Mod mod,
                                  const Graph * graph, SubGraphData * data,
                                  UInt * xout) const override;

  };

} // namespace fflow


#endif // FFLOW_SUBGRAPH_HH
