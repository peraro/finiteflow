#ifndef FFLOW_SUBGRAPH_RECONSTRUCT_HH
#define FFLOW_SUBGRAPH_RECONSTRUCT_HH

#include <fflow/subgraph.hh>

namespace fflow {


  class SubgraphRecData : public SubGraphData {
  private:
    friend class SubgraphRec;
  };

  class SubgraphRec : public SubGraph {
  public:

    Ret init(Session & session, unsigned graphid,
             SubgraphRecData & data,
             unsigned n_parsin, unsigned n_recvars,
             bool shiftvars);

    void set_default_shift();

    unsigned n_rec_vars() const
    {
      return subgraph()->nparsin[0] - nparsin[0];
    }

    const SparseRationalFunction * rec_function() const
    {
      return res_.data();
    }

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

    virtual Ret learnFromGraph(Context * ctxt,
                               const MutAlgInput xin[],  Mod mod,
                               const Graph * graph,
                               SubGraphData * data) override;

    virtual Ret evaluateFromGraph(Context * ctxt,
                                  const MutAlgInput xin[],  Mod mod,
                                  const Graph * graph, SubGraphData * data,
                                  UInt * xout) const override;

    virtual UInt min_learn_times() override { return 1; }

    void setOptions(ReconstructionOptions opt)
    {
      opt_ = opt;
    }

  private:

    void gen_sample_points_(SamplePointsVector & samples,
                            unsigned npars_in, unsigned npars_out,
                            const ReconstructionOptions & opt) const;
    void verify_sample_points_(SamplePointsVector & samples,
                               unsigned npars_in, unsigned npars_out,
                               const ReconstructionOptions & opt);
    Ret filter_points_(SamplePointsVector & samples,
                       unsigned npars_in, unsigned npars_out,
                       const ReconstructionOptions & opt) const;
    void remove_unneeded_points_(SamplePointsVector & samples,
                                 unsigned npars_in, unsigned npars_out,
                                 const ReconstructionOptions & opt) const;

  private:
    static unsigned get_start_mod_(Mod mod);
    static unsigned cache_size_(unsigned nparsin, unsigned nparsout,
                                const std::unique_ptr<UInt[]> * start,
                                const std::unique_ptr<UInt[]> * end);

  private:
    struct SubgraphToRecFun_;

  private:
    std::vector<SparseRationalFunction> res_;
    std::unique_ptr<UInt[]> shift_;
    Graph::AlgDegs degs_;
    ReconstructionOptions opt_;
    unsigned nsamples_ = 0;
    unsigned nfiltered_samples_ = 0;
    unsigned alloc_size_ = 0;
  };


} // namespace fflow


#endif // FFLOW_SUBGRAPH_RECONSTRUCT_HH
