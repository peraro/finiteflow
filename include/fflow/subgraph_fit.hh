#ifndef FFLOW_SUBGRAPH_FIT_HH
#define FFLOW_SUBGRAPH_FIT_HH

#include <fflow/subgraph.hh>
#include <fflow/alg_linear_fit.hh>

namespace fflow {

  class SubgraphFit;
  class SubgraphFitData;
  class SubgraphMultiFit;
  class SubgraphMultiFitData;


  namespace detail {

    class SubgraphSubFit;
    class SubgraphSubFitData;

    class SubgraphSubFitData : public LinearFitData  {
    private:
      friend class SubgraphSubFit;
      friend class ::fflow::SubgraphFit;
    private:
      Context * ctxt = nullptr;
      UInt * xin = nullptr;
      const Graph * graph = nullptr;
      GraphData * data = nullptr;
    };

    class SubgraphSubFit : public LinearFit {
    public:

      virtual Ret new_xpoint(Context * ctxt, AlgInput xi[], Mod mod,
                             AlgorithmData * data) const override;
      virtual Ret fill_equation(Context * ctxt,
                                AlgInput xi[], const UInt tau[], Mod mod,
                                AlgorithmData * data,
                                UInt res[]) const override;

    private:
      friend class ::fflow::SubgraphFit;
    };

  } // namespace detail



  class SubgraphFitData : public SubGraphData {
  private:
    friend class SubgraphFit;
  private:
    detail::SubgraphSubFitData impl_;
  };

  class SubgraphFit : public SubGraph {
  public:

    Ret init(Session & session, unsigned graphid,
             SubgraphFitData & data,
             const unsigned * param_nodes,
             unsigned n_params_nodes,
             std::size_t ncoeffs, std::size_t nsamplevars,
             unsigned * needed_coeffs, unsigned needed_size,
             unsigned extra_eqs);

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

    virtual UInt min_learn_times() override { return impl_.min_learn_times(); }

    const LinearFit & linear_fit() const
    {
      return impl_;
    }

  private:
    detail::SubgraphSubFit impl_;
  };


  // Multi-fit

  namespace detail {

    class SubgraphSubMultiFit;
    class SubgraphSubMultiFitData;

    class SubgraphSubMultiFitData : public LinearFitData  {
    private:
      friend class SubgraphSubMultiFit;
      friend class ::fflow::SubgraphMultiFit;
    private:
      Context * ctxt = nullptr;
      UInt * xin = nullptr;
      const Graph * graph = nullptr;
      GraphData * data = nullptr;
      UIntCache * cache = nullptr;
      UInt * gout = nullptr;
      const unsigned * take = nullptr;
      unsigned take_size = 0;
    };

    class SubgraphSubMultiFit : public LinearFit {
    public:

      virtual Ret new_xpoint(Context * ctxt, AlgInput xi[], Mod mod,
                             AlgorithmData * data) const override;
      virtual Ret fill_equation(Context * ctxt,
                                AlgInput xi[], const UInt tau[], Mod mod,
                                AlgorithmData * data,
                                UInt res[]) const override;

    private:
      friend class ::fflow::SubgraphMultiFit;
    };

  } // namespace detail



  class SubgraphMultiFitData : public SubGraphData {
  private:
    friend class SubgraphMultiFit;
  private:
    std::vector<detail::SubgraphSubMultiFitData> impl_;
    UIntCache cache_;
    std::unique_ptr<UInt[]> gout_;
  };

  class SubgraphMultiFit : public SubGraph {
  public:

    Ret init(Session & session, unsigned graphid,
             SubgraphMultiFitData & data,
             const unsigned * param_nodes,
             unsigned n_params_nodes,
             std::vector<std::vector<unsigned>> && take,
             std::size_t nsamplevars,
             std::vector<std::vector<unsigned>> & needed,
             unsigned extra_eqs);

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

    virtual UInt min_learn_times() override
    {
      return impl_.empty() ? 0 : impl_[0].min_learn_times();
    }

    const LinearFit & linear_fit(std::size_t j) const
    {
      return impl_[j];
    }

    std::size_t n_fits() const
    {
      return impl_.size();
    }

    unsigned n_sample_vars() const
    {
      return impl_.empty() ? 0 : impl_[0].n_sample_vars();
    }

  private:
    std::vector<detail::SubgraphSubMultiFit> impl_;
    std::vector<std::vector<unsigned>> take_;
  };

} // namespace fflow


#endif // FFLOW_SUBGRAPH_FIT_HH
