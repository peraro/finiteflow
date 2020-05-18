#ifndef FFLOW_GRAPH_HH
#define FFLOW_GRAPH_HH

#include <fflow/algorithm.hh>
#include <fflow/alg_mp_reconstruction.hh>
#include <fflow/refcounted_ptr.hh>
#include <future>
#if FFLOW_THREAD_POOL
# include <fflow/thread_pool.hh>
#endif

namespace fflow {

  class Node;
  class Graph;
  class GraphData;
  class Context;
  class Session;
  class SubGraph;


  const unsigned ALG_NO_ID = ~unsigned(0);


  class Node {
  public:

    Node() = default;

    ~Node();

    unsigned id() const
    {
      return id_;
    }

    const Graph * graph() const
    {
      return graph_;
    }

    Graph * graph()
    {
      return graph_;
    }

    unsigned graph_id() const;

    const Session * session() const;

    Session * session();

    Algorithm * algorithm()
    {
      return alg_.get();
    }

    const Algorithm * algorithm() const
    {
      return alg_.get();
    }

  private:
    friend class Graph;
    friend class Session;

  private:
    AlgorithmPtr alg_;
    std::vector<unsigned> inputs_;
    Graph * graph_ = nullptr;
    unsigned id_ = ALG_NO_ID;
    unsigned depth_ = 0;

  public:
    LearningOptions learn_opt;
  };


  class GraphData : public AlgorithmData {
  private:
    void add_data_(unsigned id, std::nullptr_t);
    void add_data_(unsigned id, AlgorithmData::Ptr && data);
    void set_nparsout_(unsigned id, unsigned nvars);

    void delete_data_(unsigned id);

  private:
    friend class Graph;
    friend class Context;
    friend class Session;

  private:
    std::vector<std::unique_ptr<AlgorithmData>> algdata_;
    std::vector<MallocArray<UInt>> xout_;
  };


  class Graph : public Algorithm {
  public:

    Graph() = default;

    ~Graph();

    virtual UInt min_learn_times() override;

    virtual Ret learn(Context * ctxt, Input x[], Mod,
                      AlgorithmData * data) override;

    virtual Ret evaluate(Context * ctxt,
                         Input xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    Node * node(unsigned id)
    {
      if (!node_exists(id))
        return nullptr;
      return nodes_[id].get();
    }

    const Node * node(unsigned id) const
    {
      if (!node_exists(id))
        return nullptr;
      return nodes_[id].get();
    }

    Node * out_node()
    {
      return node(out_node_);
    }

    const Node * out_node() const
    {
      return node(out_node_);
    }

    unsigned out_node_id() const
    {
      return out_node_;
    }

    // These return a 2*#edges list of the form {in1,out1,in2,out2,...}
    // The return value is the id of the output node
    unsigned edges(std::vector<unsigned> & out) const;
    unsigned marked_edges(std::vector<unsigned> & out) const;

    unsigned nodes(std::vector<unsigned> & out) const;
    unsigned marked_nodes(std::vector<unsigned> & out) const;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * adata) const override;

    bool node_exists(unsigned id) const;

    unsigned new_node(Algorithm * alg,
                      std::unique_ptr<AlgorithmData> && algdata,
                      const unsigned * inputs);
    unsigned new_node(AlgorithmPtr && alg,
                      std::unique_ptr<AlgorithmData> && algdata,
                      const unsigned * inputs);

    // This will fail unless the node is mutable
    unsigned delete_node(unsigned id);

    unsigned set_input_vars(unsigned nvars);

    bool input_vars_not_set() const
    {
      return nodes_.empty();
    }

    unsigned set_output_node(unsigned id);
    bool check_output_node(unsigned id) const
    {
      return id == out_node_;
    }

    void invalidate_reconstruction_cache();

  private:

    unsigned get_free_node_id_();

    void mark_node_(unsigned id)
    {
      mark_nodes_[id] = 1;
    }
    bool node_is_marked_(unsigned id) const
    {
      return mark_nodes_[id];
    }

    void mark_node_and_deps_(unsigned id);

    void reset_node_marks_()
    {
      mark_nodes_.resize(0);
      mark_nodes_.resize(nodes_.size(),0);
    }

    unsigned delete_node_(unsigned id, bool force);

    Ret evaluate_comput_seq_(Context * ctxt,
                             Input xinin[], Mod mod, AlgorithmData * datain,
                             bool skip_out_node) const;

    GraphData & graph_data_();

    void sync_graph_data(const GraphData & adatain,
                         GraphData & adataout,
                         bool all_nodes = false) const;

    unsigned edges_(std::vector<unsigned> & out, bool marked_only) const;
    unsigned get_nodes_(std::vector<unsigned> & out, bool marked_only) const;

  public:

    struct AlgDegs {
      std::unique_ptr<unsigned[]> numdeg, dendeg;
      std::unique_ptr<RatFunVarDegrees[]> vdegs;
    };

    struct AlgRecData {
      CacheUPtr cache;
      std::unique_ptr<UInt[]> shift;
      AlgDegs degs;
      std::unique_ptr<bool[]> completed_;
      std::unique_ptr<MPReconstructedRatFun[]> recratfun_;
    };

    AlgRecData & rec_data();
    AlgDegs & degs_data();
    const AlgRecData & rec_data() const;
    const AlgDegs & degs_data() const;

    FFLOW_IMPLEM_REF_COUNT(refcount_)

  private:
    friend class Node;
    friend class Session;
    friend class SubGraph;

  private:
    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<unsigned> comput_seq_;
    std::vector<unsigned> free_node_slots_;
    unsigned out_node_ = ALG_NO_ID;
    std::vector<char> mark_nodes_;
    std::unique_ptr<AlgRecData> rec_data_;
    Session * session_ = nullptr;
    unsigned id_ = ALG_NO_ID;
    unsigned refcount_ = 0;
    unsigned subgraph_count_ = 0;
  };

  typedef RefCountedPtr<Graph> GraphPtr;


  inline bool Graph::node_exists(unsigned id) const
  {
    return nodes_.size() > id && nodes_[id].get() != nullptr;
  }


  class Context {
  private:
    friend class Graph;
    friend class Session;
  public:
    HornerWorkspacePtr ww;

    GraphData * graph_data(unsigned i)
    {
      return graph_data_[i].get();
    }

    bool has_graph_data(unsigned i) const;

  private:
    void delete_node_data_(unsigned graphid, unsigned nodeid);
    void delete_graph_data_(unsigned graphid);
    void allocate_graph_data_(unsigned id);
  private:
    std::vector<UInt*> xin_;
    std::vector<std::unique_ptr<GraphData>> graph_data_;
  };


  class Session {
  public:
    unsigned new_graph();

    Graph * graph(unsigned id);
    const Graph * graph(unsigned id) const;
    bool graph_exists(unsigned id) const;
    bool graph_can_be_evaluated(unsigned id) const;

    Node * node(unsigned graphid, unsigned nodeid);
    const Node * node(unsigned graphid, unsigned nodeid) const;
    bool node_exists(unsigned graphid, unsigned nodeid) const;

    const Algorithm * algorithm(unsigned graphid, unsigned nodeid) const;
    Algorithm * algorithm(unsigned graphid, unsigned nodeid);

    // This returns AlgorithmData from the main context.  If this is
    // modified in an incompatible way, one should first make sure the
    // algorithm is mutable (if needed) and then also call
    // invalidate_subctxt_alg_data(graphid, nodeid).  If any applied
    // change modifies the values of the output of a graph,
    // invalidate_reconstruction_cache should also be called.
    AlgorithmData * alg_data(unsigned graphid, unsigned nodeid);
    void invalidate_subctxt_alg_data(unsigned graphid, unsigned nodeid);
    void invalidate_reconstruction_cache(unsigned graphid);

    void delete_graph(unsigned graphid);
    unsigned delete_node(unsigned graphid, unsigned nodeid);

    unsigned set_output_node(unsigned graphid, unsigned nodeid);

    // Nodes and graphs which are not currently being used by other
    // algorithms, e.g. as inputs, can be made mutable again, even if
    // they had been declared immutable before
    unsigned set_node_mutable(unsigned graphid, unsigned nodeid);

    // Remove all nodes which are not needed to compute the output
    // node
    unsigned prune_graph(unsigned graphid);

    Ret learn(unsigned graphid);
    Ret degrees(unsigned graphid, const ReconstructionOptions & opt);
    Ret all_var_degrees(unsigned graphid, const ReconstructionOptions & opt);
    Ret all_degrees(unsigned graphid, const ReconstructionOptions & opt);
    Ret parallel_all_degrees(unsigned graphid, unsigned nthreads,
                             const ReconstructionOptions & opt);

    void sample(unsigned graphid, const ReconstructionOptions & opt,
                SamplePointsGenerator * samplegen = nullptr);
    void parallel_sample(unsigned graphid, unsigned nthreads,
                         const ReconstructionOptions & opt,
                         SamplePointsGenerator * samplegen = nullptr);

    Ret reconstruct(unsigned graphid, MPReconstructedRatFun res[],
                    const ReconstructionOptions & opt);
    Ret parallel_reconstruct(unsigned graphid, MPReconstructedRatFun res[],
                             unsigned nthreads,
                             const ReconstructionOptions & opt);
    Ret reconstruct_univariate(unsigned graphid, MPReconstructedRatFun res[],
                               const ReconstructionOptions & opt);
    Ret reconstruct_numeric(unsigned graphid, MPRational res[],
                            const ReconstructionOptions & opt);

    Ret full_reconstruction(unsigned graphid, MPReconstructedRatFun res[],
                            unsigned nthreads,
                            ReconstructionOptions opt,
                            unsigned min_primes = 1);

    Ret reconstruct_mod(unsigned graphid, MPReconstructedRatFun res[],
                        const ReconstructionOptions & opt);
    Ret parallel_reconstruct_mod(unsigned graphid, MPReconstructedRatFun res[],
                                 unsigned nthreads,
                                 const ReconstructionOptions & opt);
    Ret reconstruct_univariate_mod(unsigned graphid,
                                   MPReconstructedRatFun res[],
                                   const ReconstructionOptions & opt);
    Ret full_reconstruction_mod(unsigned graphid, MPReconstructedRatFun res[],
                                unsigned nthreads,
                                ReconstructionOptions opt);

    Ret evaluate_list(unsigned graphid, const SamplePointsVector & input,
                      SamplePointsVector & res,
                      unsigned nthreads = 0);

    // Checks that the output of an algorithm is independent of the
    // variable "var".  Returns 1=True, 0=false, or FAILED
    UInt independent_of_var(unsigned graphid, unsigned var);

    Ret dump_sample_points(unsigned graphid,
                           const ReconstructionOptions & opt,
                           const char * filename);

    Ret dump_evaluations(unsigned graphid,
                         const char * filename);
    Ret load_evaluations(unsigned graphid,
                         const char ** filename,
                         unsigned n_files);

    // Creates a dummy graph which can't be evaluated from a file
    // containing information about degrees.  It can be used in
    // combination with load_degrees and load_evaluations to perform a
    // reconstruction from previously cached and stored evaluations,
    // without having to evaluate the graph.
    unsigned new_dummy_graph(unsigned nparsin, unsigned nparsout);
    Ret dump_degrees(unsigned graphid, const char * file) const;
    Ret load_degrees(unsigned graphid, const char * file);

    // Returns the total number of evaluation points needed for the
    // reconstruction, using a number of primes, etc... specified in
    // the options.  If count is not null, it will write in
    // complexity[i] the number of points needed for reconstructing
    // the i-th element of the output only.  If file is not null, the
    // sample points are reads from the specified file, otherwise they
    // are generated by computing the degrees first.
    UInt count_sample_points(unsigned graphid,
                             const ReconstructionOptions & opt,
                             unsigned count[],
                             unsigned nthreads = 0,
                             const char * file = nullptr);

    Context * main_context()
    {
      return &ctxt_;
    }

    Context * subcontext(unsigned i)
    {
      return &(sub_ctxt_[i]);
    }

    unsigned subcontext_size()
    {
      return sub_ctxt_.size();
    }

    AlgorithmData * subctxt_alg_data(unsigned subctxt,
                                     unsigned graphid,unsigned nodeid);

    void make_reconstructible(unsigned graphid);

    ~Session();

    static unsigned default_nthreads();

  private:

    void delete_node_(unsigned graphid, unsigned nodeid);
    void delete_graph_(unsigned graphid);
    void delete_node_data_(unsigned graphid, unsigned nodeid);
    void delete_graph_data_(unsigned graphid);

    unsigned get_free_graph_id_();

    void make_reconstructible_(unsigned graphid);
    void set_default_shift_(unsigned graphid);

    static void init_cache_(CacheUPtr & cache,
                            unsigned nparsin, unsigned nparsout,
                            unsigned nsamples,
                            unsigned allocsize);
    static unsigned compute_cache_size_(unsigned nparsin, unsigned nparsout,
                                        const std::unique_ptr<UInt[]> * start,
                                        const std::unique_ptr<UInt[]> * end);
    void init_graph_cache_(unsigned graphid);

    void gen_sample_points_(unsigned graphid, SamplePointsVector & samples,
                            const ReconstructionOptions & opt);

    void set_up_to_rec_(unsigned id, std::vector<unsigned> & to_rec);
    void move_rec_(unsigned id, MPReconstructedRatFun res[]);

    Ret reconstruct_(unsigned graphid,
                     unsigned var_beg, unsigned step,
                     const std::vector<unsigned> & to_reconstruct,
                     const ReconstructionOptions & opt);

    Ret reconstruct_mod_(unsigned graphid,
                         unsigned var_beg, unsigned step,
                         const std::vector<unsigned> & to_reconstruct,
                         const ReconstructionOptions & opt,
                         Mod mod);

    void share_with_subctxt_(unsigned graphid, unsigned i);
    void share_with_subctxt_(unsigned graphid, Context & sctxt);

    AlgorithmData & algdata_(unsigned graphid, unsigned nodeid);
    static AlgorithmData & algdata_(Context & ctxt,
                                    unsigned graphid, unsigned nodeid);

    Ret var_degrees_(unsigned graphid, unsigned var,
                     const ReconstructionOptions & opt);

    void init_subcontexts_(unsigned n);

    void alloc_threads_(unsigned n);

    template<typename F, typename... Ts>
    auto enqueue_(F&& f, Ts&&... args)
      -> std::future<typename std::result_of<F(Ts...)>::type>;

  private:
    friend class Node;
    friend class Graph;
    friend class SubGraph;
    friend struct GraphParallelTotDegs;
    friend struct GraphParallelVarDegs;
    friend struct GraphParallelEvaluate;
    friend struct GraphParallelReconstruct;
    friend struct GraphParallelReconstructMod;

  private:
    std::vector<GraphPtr> graphs_;
    std::vector<unsigned> free_graph_slots_;
    Context ctxt_;
    std::vector<Context> sub_ctxt_;
#if FFLOW_THREAD_POOL
    ThreadPool pool_;
#endif
  };


} // namespace fflow


#endif // FFLOW_GRAPH_HH
