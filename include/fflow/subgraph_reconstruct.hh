#ifndef FFLOW_SUBGRAPH_RECONSTRUCT_HH
#define FFLOW_SUBGRAPH_RECONSTRUCT_HH

#include <fflow/subgraph.hh>
#include <fflow/univariate_reconstruction.hh>
#include <fflow/function_cache.hh>

namespace fflow {


  // We have two different nodes for multivariate and univariate
  // SubgraphReconstruction, although this fact is hidden in the main
  // APIs (C,Python and Mathematica).


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



  // Univariate
  //
  // Since, internally, a univariate subgraph reconstruction and a
  // Laurent expansion algorithm are very similar, we implement both
  // in the same base class, controlled by a boolean flag.  We then
  // derive from it both the univariate SugraphRec algorithm and the
  // Laurent one.

  class ExpAlgorithmRatFun : public URatFun {
  public:

    void new_xpoint(UInt x[]);

    virtual UInt evaluate(const UInt x, Mod mod) override;

    void set_algorithm(const Algorithm & alg,
                       unsigned initial_sample_points = 40)
    {
      nparsout_ = std::max<unsigned>(alg.nparsout,1);
      cache_.init(1);
      cache_.reserve_more(initial_sample_points,
                          initial_sample_points*(1+nparsout_));
      alg_ = &alg;
    }

    void set_algorithm_data(AlgorithmData * alg_data)
    {
      alg_data_ = alg_data;
    }

    void set_context(Context * ctxt)
    {
      ctxt_ = ctxt;
    }

    const Algorithm * get_alg() const
    {
      return alg_;
    }

    std::size_t cache_size() const
    {
      return cache_.size();
    }

    void set_pref_exp(unsigned pref_exp)
    {
      pref_exp_ = pref_exp;
    }

  private:
    UIntCache cache_;
    UInt * x_ = nullptr;
    UInt * xin_ = nullptr;
    const Algorithm * alg_ = nullptr;
    AlgorithmData * alg_data_ = nullptr;
    Context * ctxt_ = nullptr;
    unsigned nparsout_ = 0;
    unsigned pref_exp_ = 0;

  public:
    unsigned idx = 0;
  };



  class BaseSubgraphUniRecData : public SubGraphData {
  private:
    friend class BaseSubgraphUniRec;

    void init_data_(const Algorithm & alg,
                    unsigned initial_sample_points = 40);

  private:
    URatFunNumDenRecoHighDegs rec_;
    ExpAlgorithmRatFun ufun_;
  };

  class BaseSubgraphUniRec : public SubGraph {
  public:

    Ret init(Session & session, unsigned graphid,
             BaseSubgraphUniRecData & data,
             const unsigned param_nodes[],
             unsigned n_params_nodes);

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

    virtual UInt min_learn_times() override { return 2; }

    unsigned output_size() const
    {
      return xout_size_;
    }

    const int * order() const
    {
      return order_.get();
    }

    int * order()
    {
      return order_.get();
    }

    const SparseRationalFunction * rec_function() const
    {
      return ratrec_.get();
    }

    unsigned n_rec_vars() const
    {
      return 1;
    }

    void prefactor_exponent(int pref[]) const;

  private:

    struct ExpRatFunData {
      int pref_ = 0;
      unsigned num_deg_ = 0, den_deg_ = 0;
    };

    void new_xpoint_(UInt * xin, AlgorithmData * data) const;

    unsigned n_outs_from_rec_(URationalFunction & recfun);
    Ret to_output(URationalFunction & recfun);

  private:
    std::unique_ptr<ExpRatFunData[]> rfdata_;
    std::unique_ptr<int[]> order_;
    std::unique_ptr<SparseRationalFunction[]> ratrec_;
    unsigned xout_size_ = ~unsigned(0);

  public:
    unsigned max_degree = RatFunReconstruction::DEFAULT_MAX_DEG;

  protected:
    bool laurent_expand_ = false;
  };


  typedef BaseSubgraphUniRecData SubgraphUniRecData;
  class SubgraphUniRec : public BaseSubgraphUniRec {};


} // namespace fflow


#endif // FFLOW_SUBGRAPH_RECONSTRUCT_HH
