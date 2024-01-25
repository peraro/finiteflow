#ifndef FFLOW_ALG_LAURENT_HH
#define FFLOW_ALG_LAURENT_HH

#include <fflow/algorithm.hh>
#include <fflow/subgraph.hh>
#include <fflow/univariate_reconstruction.hh>
#include <fflow/function_cache.hh>

namespace fflow {

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



  class LaurentExpansionData : public SubGraphData {
  private:
    friend class LaurentExpansion;

    void init_data_(const Algorithm & alg,
                    unsigned initial_sample_points = 40);

  private:
    URatFunNumDenRecoHighDegs rec_;
    ExpAlgorithmRatFun ufun_;
  };

  class LaurentExpansion : public SubGraph {
  public:

    Ret init(Session & session, unsigned graphid,
             LaurentExpansionData & data,
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

    void prefactor_exponent(int pref[]) const;

  private:

    struct ExpRatFunData {
      int pref_ = 0;
      unsigned num_deg_ = 0, den_deg_ = 0;
    };

    void new_xpoint_(UInt * xin, AlgorithmData * data) const;

    Ret thiele_rec_(Context * ctxt, Input xin[], Mod mod,
                    AlgorithmData * data) const;

    Ret numden_rec_(Context * ctxt, Input xin[], Mod mod,
                    AlgorithmData * data) const;

  private:
    std::unique_ptr<ExpRatFunData[]> rfdata_;
    std::unique_ptr<int[]> order_;
    unsigned xout_size_ = ~unsigned(0);

  public:
    unsigned max_degree = 100;
  };

} // namespace fflow

#endif // FFLOW_ALG_LAURENT_HH
