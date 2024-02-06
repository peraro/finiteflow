#include <algorithm>
#include <cstring>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include <stdlib.h>
#include <fflow/common.hh>
#include <fflow/primes.hh>
#include <fflow/capi.h>
#include <fflow/graph.hh>
#include <fflow/subgraph.hh>
#include <fflow/mp_common.hh>
#include <fflow/alg_functions.hh>
#include <fflow/analytic_solver.hh>
#include <fflow/numeric_solver.hh>
#include <fflow/node_solver.hh>
#include <fflow/alg_laurent.hh>
#include <fflow/alg_lists.hh>
#include <fflow/alg_mp_reconstruction.hh>
#include <fflow/mp_functions.hh>
#include <fflow/json.hh>
#include <fflow/ratfun_parser.hh>
#include <fflow/algorithm.hh>
using namespace fflow;

#define FF_MIN_ERROR (FF_ERROR - 10)

static fflow::Session session;

bool ffIsError(unsigned val)
{
  return val >= FF_MIN_ERROR;
}

template <typename IntT>
static unsigned * newU32Array(const IntT * vals, unsigned n_vals)
{
  unsigned * ret = (unsigned*)malloc(n_vals*sizeof(unsigned));
  std::copy(vals, vals+n_vals, ret);
  return ret;
}

#if 0
static char * toCStr(std::string str)
{
  char * ret = (char*)malloc(str.size() + 1);
  strcpy(ret, str.c_str());
  return ret;
}
#endif

static ReconstructionOptions toRecOpt(FFRecOptions options)
{
  ReconstructionOptions opt;
  opt.n_checks = RatFunReconstruction::DEFAULT_N_CHECKS;
  opt.n_uchecks = RatFunReconstruction::DEFAULT_N_UCHECKS;
  opt.n_singular = 0;
  opt.start_mod = options.start_mod;
  opt.max_primes = options.max_primes;
  opt.max_deg = options.max_deg;
  opt.dbginfo = options.dbginfo;
  opt.polymethod = options.polymethod;
  return opt;
}

// We check that columns are sorted and they don't exceed the matrix
// size
static FFStatus CopyNonZeroColumnEls(unsigned n_vars,
                                     const unsigned * elems,
                                     unsigned n_elems,
                                     unsigned * dest)
{
  int last_var = -1;
  for (unsigned j=0; j<n_elems; ++j, ++dest) {
    if (int(elems[j]) <= last_var || elems[j] >= n_vars + 1)
      return FF_ERROR;
    last_var = elems[j];
    *dest = elems[j];
  }
  return FF_SUCCESS;
}

static FFStatus CopyNonZeroColumnElsUnordered(unsigned n_cols,
                                              const unsigned * elems,
                                              unsigned n_elems,
                                              unsigned * dest)
{
  for (unsigned j=0; j<n_elems; ++j, ++dest) {
    if (elems[j] >= n_cols)
      return FF_ERROR;
    *dest = elems[j];
  }
  return FF_SUCCESS;
}


static FFStatus getPoly(unsigned n_vars,
                        unsigned n_terms,
                        FFCStr * coefficients,
                        const uint16_t * exponents,
                        MPReconstructedPoly & num)
{
  auto & mons = num.monomials_vec_();
  auto & coeffs = num.coeff_vec_();

  mons.resize(n_terms);
  coeffs.resize(n_terms);

  for (unsigned k=0; k<n_terms; ++k, ++coefficients) {

    if (!(coeffs[k].set(*coefficients) == 0))
      return FF_ERROR;

    auto & mon = mons[k];
    mon = Monomial(n_vars);
    mon.coeff() = 1;

    unsigned deg = 0;
    for (unsigned l=0; l<n_vars; ++l, ++exponents) {
      unsigned exp = mon.exponent(l) = *exponents;
      deg += exp;
    }
    mon.degree() = deg;
  }

  return FF_SUCCESS;
}


extern "C" {

  struct FFRatFunList {
    std::unique_ptr<MPReconstructedRatFun[]> rf;
    unsigned n_functions = 0;
    unsigned n_vars = 0;
  };

} // extern "C"


namespace {

  template <typename FunMap>
  void get_horner_ratfun(const FFRatFunList * rf,
                         unsigned idx,
                         HornerRatFunPtr & f,
                         FunMap & map,
                         HornerWorkspacePtr & workspace)
  {
    const MPReconstructedRatFun & rfun = rf->rf[idx];
    const MPReconstructedPoly & num = rfun.numerator();
    const MPReconstructedPoly & den = rfun.denominator();
    unsigned nterms_num = num.size();
    unsigned nterms_den = den.size();

    // numerator
    map.num_map.resize(nterms_num);
    auto num_c = map.num_map.coeff.get();
    for (unsigned j=0; j<nterms_num; ++j)
      num_c[j] = num.coeff(j);


    // denominator
    map.den_map.resize(nterms_den);
    auto den_c = map.den_map.coeff.get();
    for (unsigned j=0; j<nterms_den; ++j)
      den_c[j] = den.coeff(j);

    // build result
    f.from_sparse_poly(num.monomials(), num.size(),
                       den.monomials(), den.size(),
                       rf->n_vars, 0,
                       map.num_map.pos.get(), map.den_map.pos.get());

    std::size_t max_workspace = horner_required_workspace(f.num());
    max_workspace = std::max(horner_required_workspace(f.den()),
                             max_workspace);
    workspace.ensure_size(max_workspace);
  }

} // namespace


extern "C" {

  void ffInit(void)
  {
    // ...nothing for now...
  }

  void ffDeinit(void)
  {
    // ...nothing for now...
  }

  FFUInt ffMulInv(FFUInt z, unsigned prime_no)
  {
    if (prime_no >= BIG_UINT_PRIMES_SIZE)
      return FF_FAILED;
    Mod mod(BIG_UINT_PRIMES[prime_no]);
    return mul_inv(red_mod(z, mod), mod);
  }

  FFUInt ffPrimeNo(unsigned i)
  {
    if (i >= BIG_UINT_PRIMES_SIZE)
      return FF_FAILED;
    return BIG_UINT_PRIMES[i];
  }

  unsigned ffDefaultNThreads(void)
  {
    return Session::default_nthreads();
  }

  void ffFreeMemoryU32(unsigned * mem)
  {
    free(mem);
  }

  void ffFreeMemoryS32(int * mem)
  {
    free(mem);
  }

  void ffFreeMemoryU64(uint64_t * mem)
  {
    free(mem);
  }

  FFGraph ffNewGraph(void)
  {
    return session.new_graph();
  }

  FFGraph ffNewGraphWithInput(unsigned nvars, FFNode * node)
  {
    FFGraph graph = ffNewGraph();
    *node = ffSetGraphInput(graph, nvars);
    return graph;
  }

  FFGraph ffNewGraphDummy(unsigned n_in, unsigned n_out)
  {
    return session.new_dummy_graph(n_in, n_out);
  }

  FFStatus ffDeleteGraph(FFGraph graph)
  {
    if (!session.graph_exists(graph))
      return FF_ERROR;
    session.delete_graph(graph);
    return FF_SUCCESS;
  }

  FFStatus ffDeleteNode(FFGraph graph, FFNode node)
  {
    if (session.delete_node(graph, node) == ALG_NO_ID)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFStatus ffSetOutputNode(FFGraph graph, FFNode node)
  {
    if (session.set_output_node(graph, node) == ALG_NO_ID)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFNode ffSetGraphInput(FFGraph graph, unsigned n_vars)
  {
    Graph * graphptr = session.graph(graph);
    if (!graphptr || graphptr->set_input_vars(n_vars) == ALG_NO_ID)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  unsigned ffGraphNParsOut(FFGraph graph)
  {
    Graph * g = session.graph(graph);
    if (g == nullptr || !g->has_learned())
      return FF_ERROR;
    return g->nparsout;
  }

  unsigned ffNodeNParsOut(FFGraph graph, FFNode node)
  {
    Node * n = session.node(graph, node);
    if (n == nullptr || !n->algorithm()->has_learned())
      return FF_ERROR;
    return n->algorithm()->nparsout;
  }

  FFStatus ffMakeNodeMutable(FFGraph graph, FFNode node)
  {
    if (session.set_node_mutable(graph, node) == ALG_NO_ID)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFStatus ffPruneGraph(FFGraph graph)
  {
    if (session.prune_graph(graph) == ALG_NO_ID)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFStatus ffLearn(FFGraph graph)
  {
    if (session.learn(graph) != SUCCESS)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFUInt * ffEvaluateGraph(FFGraph graph,
                           const FFUInt * input, unsigned prime_no)
  {
    Graph * g = session.graph(graph);
    if (!session.graph_can_be_evaluated(graph) ||
        prime_no >= BIG_UINT_PRIMES_SIZE)
      return 0;

    Mod mod(BIG_UINT_PRIMES[prime_no]);
    unsigned nparsout = g->nparsout;

    FFUInt * output = (FFUInt*)malloc(sizeof(FFUInt)*nparsout);

    Context * ctxt = session.main_context();
    Ret ret = g->evaluate(ctxt, &input, mod, ctxt->graph_data(graph), output);

    if (ret != SUCCESS) {
      free(output);
      return 0;
    }

    return output;
  }

  unsigned ffSubgraphNParsout(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg)
      return FF_ERROR;

    SubGraph * subalg = dynamic_cast<SubGraph*>(alg);
    if (!subalg)
      return FF_NO_ALGORITHM;

    return subalg->subgraph()->nparsout;
  }

  int * ffLaurentMaxOrders(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);
    LaurentExpansion * lauralg = dynamic_cast<LaurentExpansion*>(alg);
    if (!lauralg)
      return 0;

    const unsigned nout = lauralg->subgraph()->nparsout;
    int * res = (int *)malloc(sizeof(int)*nout);
    std::copy(lauralg->order(), lauralg->order() + nout, res);

    return res;
  }

  int * ffLaurentMinOrders(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);
    LaurentExpansion * lauralg = dynamic_cast<LaurentExpansion*>(alg);
    if (!lauralg || !lauralg->has_learned())
      return 0;

    const unsigned nout = lauralg->subgraph()->nparsout;
    int * res = (int *)malloc(sizeof(int)*nout);
    lauralg->prefactor_exponent(res);

    return res;
  }



  // Algorithms

  FFNode ffAlgSimpleSubgraph(FFGraph graph,
                             const FFNode * in_nodes, unsigned n_in_nodes,
                             FFGraph subgraph)
  {
    typedef SubGraphData Data;
    std::unique_ptr<SimpleSubGraph> algptr(new SimpleSubGraph());
    std::unique_ptr<Data> data(new Data());
    std::vector<unsigned> nparsin(n_in_nodes);

    for (unsigned i=0; i<n_in_nodes; ++i) {
      Node * node = session.node(graph, in_nodes[i]);
      if (!node)
        return FF_ERROR;
      nparsin[i] = node->algorithm()->nparsout;
    }

    Ret ret = algptr->init(session, subgraph, *data,
                           nparsin.data(), nparsin.size());

    if (ret != SUCCESS)
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(data), in_nodes);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return FF_SUCCESS;
  }

  FFNode ffAlgMemoizedSubgraph(FFGraph graph,
                               const FFNode * in_nodes, unsigned n_in_nodes,
                               FFGraph subgraph)
  {
    typedef MemoizedSubGraphData Data;
    std::unique_ptr<MemoizedSubGraph> algptr(new MemoizedSubGraph());
    std::unique_ptr<Data> data(new Data());
    std::vector<unsigned> nparsin(n_in_nodes);

    for (unsigned i=0; i<n_in_nodes; ++i) {
      Node * node = session.node(graph, in_nodes[i]);
      if (!node)
        return FF_ERROR;
      nparsin[i] = node->algorithm()->nparsout;
    }

    Ret ret = algptr->init(session, subgraph, *data,
                           nparsin.data(), nparsin.size());

    if (ret != SUCCESS)
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(data), in_nodes);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return FF_SUCCESS;
  }

  FFNode ffAlgJSONSparseLSolve(FFGraph graph, FFNode in_node,
                               const char * json_file)
  {
    typedef AnalyticSparseSolverData Data;
    std::unique_ptr<AnalyticSparseSolver> algptr(new AnalyticSparseSolver());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;

    unsigned needed_workspace;
    Ret ret = json_sparse_system(json_file, alg, *data, needed_workspace);

    if (ret != SUCCESS || !session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(data), &in_node);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    session.main_context()->ww.ensure_size(needed_workspace);
    return id;
  }

  FFNode ffAlgAnalyticSparseLSolve(FFGraph graph, FFNode in_node,
                                   unsigned n_eqs, unsigned n_vars,
                                   const unsigned * n_non_zero,
                                   const unsigned * non_zero_els,
                                   const FFRatFunList * non_zero_coeffs,
                                   const unsigned * needed_vars,
                                   unsigned n_needed_vars)
  {
    typedef AnalyticSparseSolverData Data;
    std::unique_ptr<AnalyticSparseSolver> algptr(new AnalyticSparseSolver());
    std::unique_ptr<Data> dataptr(new Data());
    auto & sys = *algptr;
    auto & data = *dataptr;

    const unsigned n_rows = n_eqs;
    sys.rinfo.resize(n_rows);
    data.c.resize(n_rows);
    sys.cmap.resize(n_rows);

    const unsigned n_functions = non_zero_coeffs->n_functions;
    unsigned idx = 0;

    for (int i=0; i<n_rows; ++i) {
      const unsigned csize = n_non_zero[i];
      AnalyticSparseSolver::RowInfo & rinf = sys.rinfo[i];
      rinf.size = csize;
      rinf.cols.reset(new unsigned[csize]);
      data.c[i].reset(new HornerRatFunPtr[csize]);
      sys.cmap[i].reset(new MPHornerRatFunMap[csize]);
      FFStatus ret = CopyNonZeroColumnEls(n_vars, non_zero_els, csize,
                                          rinf.cols.get());
      if (ret != FF_SUCCESS)
        return ret;
      non_zero_els += csize;
      if (idx + csize > n_functions)
        return FF_ERROR;
      for (int j=0; j<csize; ++j, ++idx) {
        get_horner_ratfun(non_zero_coeffs, idx, data.c[i][j], sys.cmap[i][j],
                          session.main_context()->ww);
      }
    }

    sys.nparsin.resize(1);
    sys.nparsin[0] = non_zero_coeffs->n_vars;

    if (needed_vars) {
      sys.init(n_eqs, n_vars, needed_vars, n_needed_vars, data);
    } else {
      unsigned * needed = (unsigned*)malloc(n_vars*sizeof(n_vars));
      std::iota(needed, needed+n_vars, 0);
      sys.init(n_eqs, n_vars, needed, n_vars, data);
      free(needed);
    }

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(dataptr),
                              &in_node);
    return id;
  }

  FFNode ffAlgNumericSparseLSolve(FFGraph graph,
                                  unsigned n_eqs, unsigned n_vars,
                                  const unsigned * n_non_zero,
                                  const unsigned * non_zero_els,
                                  FFCStr * non_zero_coeffs,
                                  const unsigned * needed_vars,
                                  unsigned n_needed_vars)
  {
    typedef NumericSparseSolverData Data;
    std::unique_ptr<NumericSparseSolver> algptr(new NumericSparseSolver());
    std::unique_ptr<Data> dataptr(new Data());
    auto & sys = *algptr;
    auto & data = *dataptr;

    const unsigned n_rows = n_eqs;
    sys.rinfo.resize(n_rows);
    sys.c.resize(n_rows);

    unsigned idx = 0;

    for (int i=0; i<n_rows; ++i) {
      const unsigned csize = n_non_zero[i];
      NumericSparseSolver::RowInfo & rinf = sys.rinfo[i];
      rinf.size = csize;
      rinf.cols.reset(new unsigned[csize]);
      sys.c[i].reset(new MPRational[csize]);
      FFStatus ret = CopyNonZeroColumnEls(n_vars, non_zero_els, csize,
                                          rinf.cols.get());
      if (ret != FF_SUCCESS)
        return ret;
      non_zero_els += csize;
      for (int j=0; j<csize; ++j, ++idx)
        sys.c[i][j] = MPRational(non_zero_coeffs[idx]);
    }

    sys.nparsin.clear();

    if (needed_vars) {
      sys.init(n_eqs, n_vars, needed_vars, n_needed_vars, data);
    } else {
      unsigned * needed = (unsigned*)malloc(n_vars*sizeof(n_vars));
      std::iota(needed, needed+n_vars, 0);
      sys.init(n_eqs, n_vars, needed, n_vars, data);
      free(needed);
    }

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(dataptr),
                              nullptr);
    return id;
  }

  FFNode ffAlgNodeSparseLSolve(FFGraph graph, FFNode in_node,
                               unsigned n_eqs, unsigned n_vars,
                               const unsigned * n_non_zero,
                               const unsigned * non_zero_els,
                               const unsigned * needed_vars,
                               unsigned n_needed_vars)
  {
    typedef NodeSparseSolverData Data;
    std::unique_ptr<NodeSparseSolver> algptr(new NodeSparseSolver());
    std::unique_ptr<Data> dataptr(new Data());
    auto & sys = *algptr;
    auto & data = *dataptr;

    const unsigned n_rows = n_eqs;
    sys.rinfo.resize(n_rows);

    unsigned idx = 0;

    for (int i=0; i<n_rows; ++i) {
      const unsigned csize = n_non_zero[i];
      NodeSparseSolver::RowInfo & rinf = sys.rinfo[i];
      rinf.start = idx;
      rinf.size = csize;
      rinf.cols.reset(new unsigned[csize]);
      FFStatus ret = CopyNonZeroColumnEls(n_vars, non_zero_els, csize,
                                          rinf.cols.get());
      if (ret != FF_SUCCESS)
        return ret;
      non_zero_els += csize;
      idx += csize;
    }

    sys.nparsin.resize(1);
    sys.nparsin[0] = idx;

    if (needed_vars) {
      sys.init(n_eqs, n_vars, needed_vars, n_needed_vars, data);
    } else {
      unsigned * needed = (unsigned*)malloc(n_vars*sizeof(n_vars));
      std::iota(needed, needed+n_vars, 0);
      sys.init(n_eqs, n_vars, needed, n_vars, data);
      free(needed);
    }

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(dataptr),
                              &in_node);
    return id;
  }

  FFNode ffAlgJSONRatFunEval(FFGraph graph, FFNode in_node,
                             const char * json_file)
  {
    typedef AnalyticFunctionData Data;
    std::unique_ptr<AnalyticFunction> algptr(new AnalyticFunction());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;

    unsigned needed_workspace;
    Ret ret = json_sparse_ratfun(json_file, alg, *data, needed_workspace);

    if (ret != SUCCESS || !session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(data), &in_node);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    session.main_context()->ww.ensure_size(needed_workspace);
    return id;
  }

  FFNode ffAlgLaurent(FFGraph graph, FFNode in_node, FFGraph subgraph,
                      const int * order, int max_deg)
  {
    typedef LaurentExpansionData Data;
    std::unique_ptr<LaurentExpansion> algptr(new LaurentExpansion());
    std::unique_ptr<Data> data(new Data());

    unsigned nparsin = 0;
    {
      Node * node = session.node(graph, in_node);
      if (!node)
        return FF_ERROR;
      nparsin = node->algorithm()->nparsout;
    }

    Ret ret = algptr->init(session, subgraph, *data, &nparsin, 1);

    unsigned order_size = ffGraphNParsOut(subgraph);
    if (ffIsError(order_size) || ret != SUCCESS)
      return FF_ERROR;

    std::copy(order, order+order_size, algptr->order());

    if (max_deg >= 0)
      algptr->max_degree = max_deg;

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(data), &in_node);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }

  FFNode ffAlgLaurentConstOrder(FFGraph graph, FFNode in_node, FFGraph subgraph,
                                int order, int max_deg)
  {
    unsigned nout = ffGraphNParsOut(subgraph);
    if (ffIsError(nout))
      return FF_ERROR;
    std::vector<int> orders(nout);
    std::fill(orders.begin(), orders.end(), order);
    return ffAlgLaurent(graph, in_node, subgraph, orders.data(), max_deg);
  }

  FFNode ffAlgMatMul(FFGraph graph, FFNode in_node_a, FFNode in_node_b,
                     unsigned n_rows_a, unsigned n_cols_a, unsigned n_cols_b)
  {
    std::unique_ptr<MatrixMul> algptr(new MatrixMul());
    auto & alg = *algptr;

    alg.init(n_rows_a, n_cols_a, n_cols_b);

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    FFNode input[2] = {in_node_a, in_node_b};
    unsigned id = g->new_node(std::move(algptr), nullptr, input);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }

  FFNode ffAlgChain(FFGraph graph,
                    const FFNode * in_nodes, unsigned n_in_nodes)
  {
    std::unique_ptr<Chain> algptr(new Chain());
    Chain & alg = *algptr;
    std::vector<unsigned> nparsin(n_in_nodes);
    for (unsigned i=0; i<n_in_nodes; ++i) {
      Node * node = session.node(graph, in_nodes[i]);
      if (!node)
        return FF_ERROR;
      nparsin[i] = node->algorithm()->nparsout;
    }
    alg.init(nparsin.data(), nparsin.size());

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), nullptr, in_nodes);

    return id;
  }

  FFNode ffAlgTake(FFGraph graph,
                   const FFNode * in_nodes, unsigned n_in_nodes,
                   const unsigned * elems, unsigned n_elems)
  {
    std::unique_ptr<Take> algptr(new Take());
    Take & alg = *algptr;
    std::vector<unsigned> nparsin(n_in_nodes);
    for (unsigned i=0; i<n_in_nodes; ++i) {
      Node * node = session.node(graph, in_nodes[i]);
      if (!node)
        return FF_ERROR;
      nparsin[i] = node->algorithm()->nparsout;
    }

    std::vector<Take::InputEl> els;
    els.resize(n_elems);
    for (int j=0; j<n_elems; ++j, elems += 2) {

      els[j].list = elems[0];
      els[j].el = elems[1];

      // checks
      if (elems[0] >= n_in_nodes)
        return FF_ERROR;
      Node * node = session.node(graph, in_nodes[elems[0]]);
      if (!node || node->algorithm()->nparsout <= elems[1])
        return FF_ERROR;
    }

    alg.init(nparsin.data(), nparsin.size(), std::move(els));

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), nullptr, in_nodes);

    return id;
  }

  FFNode ffAlgSlice(FFGraph graph, FFNode in_node,
                    unsigned begin, int end)
  {
    std::unique_ptr<Slice> algptr(new Slice());
    Slice & alg = *algptr;

    Node * node = nullptr;
    if (!(node = session.node(graph, in_node)))
      return FF_ERROR;

    unsigned nparsin = node->algorithm()->nparsout;
    if (end < 0)
      end = nparsin;

    Ret ret = alg.init(nparsin, begin, end);
    if (ret != SUCCESS)
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), nullptr, &in_node);

    return id;
  }

  FFNode ffAlgAdd(FFGraph graph,
                  const FFNode * in_nodes, unsigned n_in_nodes)
  {
    std::unique_ptr<Add> algptr(new Add());
    Add & alg = *algptr;
    std::vector<unsigned> nparsin(n_in_nodes);
    unsigned len = 0;
    for (unsigned i=0; i<n_in_nodes; ++i) {
      Node * node = session.node(graph, in_nodes[i]);
      if (!node)
        return FF_ERROR;
      nparsin[i] = node->algorithm()->nparsout;
      if (i != 0) {
        if (nparsin[i] != len)
          return FF_ERROR;
      } else {
        len = nparsin[0];
      }
    }

    alg.init(nparsin.size(), len);

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), nullptr, in_nodes);

    return id;
  }

  FFNode ffAlgMul(FFGraph graph,
                  const FFNode * in_nodes, unsigned n_in_nodes)
  {
    std::unique_ptr<Mul> algptr(new Mul());
    Mul & alg = *algptr;
    std::vector<unsigned> nparsin(n_in_nodes);
    unsigned len = 0;
    for (unsigned i=0; i<n_in_nodes; ++i) {
      Node * node = session.node(graph, in_nodes[i]);
      if (!node)
        return FF_ERROR;
      nparsin[i] = node->algorithm()->nparsout;
      if (i != 0) {
        if (nparsin[i] != len)
          return FF_ERROR;
      } else {
        len = nparsin[0];
      }
    }

    alg.init(nparsin.size(), len);

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), nullptr, in_nodes);

    return id;
  }

  FFNode ffAlgTakeAndAdd(FFGraph graph,
                         const FFNode * in_nodes, unsigned n_in_nodes,
                         unsigned n_elems,
                         const unsigned * elems_len,
                         const unsigned * elems)
  {
    std::unique_ptr<TakeAndAdd> algptr(new TakeAndAdd());
    TakeAndAdd & alg = *algptr;
    std::vector<unsigned> nparsin(n_in_nodes);
    for (unsigned i=0; i<n_in_nodes; ++i) {
      Node * node = session.node(graph, in_nodes[i]);
      if (!node)
        return FF_ERROR;
      nparsin[i] = node->algorithm()->nparsout;
    }

    std::vector<std::vector<TakeAndAdd::InputEl>> els;
    els.resize(n_elems);
    for (int j=0; j<n_elems; ++j) {

      const unsigned this_len = elems_len[j];
      els[j].resize(this_len);

      for (int k=0; k<this_len; ++k, elems += 2) {

        els[j][k].list = elems[0];
        els[j][k].el = elems[1];

        // checks
        if (elems[0] >= n_in_nodes)
          return FF_ERROR;
        Node * node = session.node(graph, in_nodes[elems[0]]);
        if (!node || node->algorithm()->nparsout <= elems[1])
          return FF_ERROR;
      }
    }

    alg.init(nparsin.data(), nparsin.size(), std::move(els));

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), nullptr, in_nodes);

    return id;
  }

  FFNode ffAlgSparseMatMul(FFGraph graph, FFNode in_node_a, FFNode in_node_b,
                           unsigned n_rows_a, unsigned n_cols_a,
                           unsigned n_cols_b,
                           const unsigned * n_non_zero_a,
                           const unsigned * non_zero_els_a,
                           const unsigned * n_non_zero_b,
                           const unsigned * non_zero_els_b)
  {
    std::unique_ptr<SparseMatrixMul> algptr(new SparseMatrixMul());
    auto & alg = *algptr;

    alg.row1.resize(n_rows_a);
    for (unsigned j=0; j<n_rows_a; ++j) {
      unsigned len = alg.row1[j].size = n_non_zero_a[j];
      alg.row1[j].cols.reset(new unsigned[len]);
      unsigned * rcols = alg.row1[j].cols.get();
      Ret ret = CopyNonZeroColumnElsUnordered(n_cols_a,
                                              non_zero_els_a, len, rcols);
      if (ret != FF_SUCCESS)
        return FF_ERROR;
      non_zero_els_a += len;
    }

    const unsigned n_rows_b = n_cols_a;
    alg.row2.resize(n_rows_b);
    for (unsigned j=0; j<n_rows_b; ++j) {
      unsigned len = alg.row2[j].size = n_non_zero_b[j];
      alg.row2[j].cols.reset(new unsigned[len]);
      unsigned * rcols = alg.row2[j].cols.get();
      Ret ret = CopyNonZeroColumnElsUnordered(n_cols_b,
                                              non_zero_els_b, len, rcols);
      if (ret != FF_SUCCESS)
        return FF_ERROR;
      non_zero_els_b += len;
    }

    alg.init(n_rows_a, n_cols_a, n_cols_b);

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    FFNode input[2] = {in_node_a, in_node_b};
    unsigned id = g->new_node(std::move(algptr), nullptr, input);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }


  // Sparse systems

  FFStatus ffLSolveResetNeededVars(FFGraph graph, FFNode node,
                                   const unsigned * vars, unsigned n_vars)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg || !alg->is_mutable())
      return FF_ERROR;

    if (dynamic_cast<DenseLinearSolver *>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver *>(alg);
      Ret ret = ls.reset_needed(session.alg_data(graph, node), vars, n_vars);
      if (ret == SUCCESS)
        session.invalidate_subctxt_alg_data(graph, node);
      else
        return FF_ERROR;

    } else if(dynamic_cast<SparseLinearSolver *>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
      if (ls.marked_and_sweeped()) {
        return FF_ERROR;
      } else {
        Ret ret = ls.reset_needed(session.alg_data(graph, node), vars, n_vars);
        if (ret == SUCCESS)
          session.invalidate_subctxt_alg_data(graph, node);
        else
          return FF_ERROR;
      }

    } else {
      return FF_ERROR;
    }

    return FF_SUCCESS;
  }

  FFStatus ffLSolveOnlyHomogeneous(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg || !alg->is_mutable())
      return FF_ERROR;


    if (dynamic_cast<DenseLinearSolver *>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver *>(alg);
      ls.only_homogeneous();
      session.invalidate_subctxt_alg_data(graph, node);

    } else if (dynamic_cast<SparseLinearSolver *>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
      ls.only_homogeneous();
      session.invalidate_subctxt_alg_data(graph, node);

    } else {
      return FF_ERROR;
    }

    return FF_SUCCESS;
  }

  FFStatus ffLSolveSparseOutput(FFGraph graph, FFNode node, bool sparse)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg || !alg->is_mutable())
      return FF_ERROR;

    if (dynamic_cast<SparseLinearSolver *>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
      ls.sparse_output(sparse);
      session.invalidate_subctxt_alg_data(graph, node);
    } else {
      return FF_ERROR;
    }

    return FF_SUCCESS;
  }

  FFStatus ffLSolveMarkAndSweepEqs(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg || !alg->is_mutable())
      return FF_ERROR;

    if (dynamic_cast<SparseLinearSolver *>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
      if (ls.marked_and_sweeped())
        return FF_ERROR;
      ls.mark_and_sweep_eqs(session.alg_data(graph, node));
      session.invalidate_subctxt_alg_data(graph, node);
    } else {
      return FF_ERROR;
    }

    return FF_SUCCESS;
  }

  FFStatus ffLSolveDeleteUnneededEqs(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg) {
      return FF_ERROR;

    } else if (dynamic_cast<AnalyticSparseSolver *>(alg)) {
      AnalyticSparseSolver & ls = *static_cast<AnalyticSparseSolver *>(alg);
      ls.delete_unneeded_eqs(session.alg_data(graph, node));
      session.invalidate_subctxt_alg_data(graph, node);

    } else if (dynamic_cast<NumericSparseSolver *>(alg)) {
      NumericSparseSolver & ls = *static_cast<NumericSparseSolver *>(alg);
      ls.delete_unneeded_eqs(session.alg_data(graph, node));
      session.invalidate_subctxt_alg_data(graph, node);

    } else {
      return FF_ERROR;
    }

    return FF_SUCCESS;
  }

  unsigned ffLSolveNDepVars(FFGraph graph, FFNode node)
  {
    if (!session.node_exists(graph,node))
      return FF_ERROR;

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
      return ls.n_needed_depvars();
    } else if(dynamic_cast<SparseLinearSolver*>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
      return ls.n_needed_depvars();
    }

    return FF_ERROR;
  }

  unsigned * ffLSolveDepVars(FFGraph graph, FFNode node)
  {
    if (!session.node_exists(graph,node))
      return 0;

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
      return newU32Array(ls.needed_depvars(), ls.n_needed_depvars());
    } else if(dynamic_cast<SparseLinearSolver*>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
      return newU32Array(ls.needed_depvars(), ls.n_needed_depvars());
    }

    return 0;
  }

  unsigned ffLSolveNIndepVars(FFGraph graph, FFNode node, unsigned i)
  {
    if (!session.node_exists(graph,node))
      return FF_ERROR;

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
      return ls.n_needed_indepvars();
    } else if(dynamic_cast<SparseLinearSolver*>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
      if (!ls.output_is_sparse()) {
        return ls.n_needed_indepvars();
      } else {
        const auto & spdata = *(ls.sparse_output_data());
        if (i >= spdata.size())
          return FF_ERROR;
        return spdata[i].size();
      }
    }

    return FF_ERROR;
  }

  unsigned * ffLSolveIndepVars(FFGraph graph, FFNode node, unsigned i)
  {
    if (!session.node_exists(graph,node))
      return 0;

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
      return newU32Array(ls.needed_indepvars(), ls.n_needed_indepvars());
    } else if(dynamic_cast<SparseLinearSolver*>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
      if (!ls.output_is_sparse()) {
        return newU32Array(ls.needed_indepvars(), ls.n_needed_indepvars());
      } else {
        const auto & spdata = *(ls.sparse_output_data());
        if (i >= spdata.size())
          return 0;
        return newU32Array(spdata[i].data(), spdata[i].size());
      }
    }

    return 0;
  }


  // Rational functions

  void ffFreeRatFun(FFRatFunList * rf)
  {
    delete rf;
  }

  unsigned ffRatFunListSize(const FFRatFunList * rf)
  {
    return rf->n_functions;
  }

  unsigned ffRatFunListNVars(const FFRatFunList * rf)
  {
    return rf->n_vars;
  }

  FFStatus ffRatFunToJSON(const FFRatFunList * rf, const char * file)
  {
    Ret ret = json_write_ratfun(file, rf->rf.get(),
                                rf->n_functions, rf->n_vars);
    if (ret != SUCCESS)
      return FF_ERROR;

    return FF_SUCCESS;
  }


  FFRatFunList * ffParseRatFun(FFCStr * vars, unsigned n_vars,
                               FFCStr * inputs, unsigned n_functions)
  {
    unsigned * in_strlen = (unsigned*)malloc(n_functions * sizeof(unsigned));
    for (unsigned j=0; j<n_functions; ++j)
      in_strlen[j] = strlen(inputs[j]);
    FFRatFunList * ret = ffParseRatFunEx(vars, n_vars,
                                         inputs, in_strlen, n_functions);
    free(in_strlen);
    return ret;
  }


  FFRatFunList * ffParseRatFunEx(FFCStr * vars, unsigned n_vars,
                                 FFCStr * inputs, const unsigned * input_strlen,
                                 unsigned n_functions)
  {
    typedef MPReconstructedRatFun ResT;
    std::unique_ptr<ResT[]> res (new ResT[n_functions]);

    std::unique_ptr<std::string[]> svars(new std::string[n_vars]);
    for (unsigned j=0; j<n_vars; ++j) {
      svars[j] = vars[j];
    }

    for (unsigned j=0; j<n_functions; ++j) {
      res[j] = MPReconstructedRatFun(n_vars);
      Ret ret = parse_rat_fun(svars.get(), n_vars,
                              inputs[j], inputs[j] + input_strlen[j],
                              res[j].numerator().monomials_vec_(),
                              res[j].numerator().coeff_vec_(),
                              res[j].denominator().monomials_vec_(),
                              res[j].denominator().coeff_vec_());
      if (ret != SUCCESS)
        return 0;
    }

    FFRatFunList * ret = new FFRatFunList();
    ret->n_functions = n_functions;
    ret->n_vars = n_vars;
    ret->rf = std::move(res);

    return ret;
  }


  FFRatFunList * ffNewRatFunList(unsigned n_vars, unsigned n_functions,
                                 const unsigned * n_num_terms,
                                 const unsigned * n_den_terms,
                                 FFCStr * coefficients,
                                 const uint16_t * exponents)
  {
    typedef MPReconstructedRatFun ResT;
    std::unique_ptr<ResT[]> res (new ResT[n_functions]);

    for (unsigned j=0; j<n_functions; ++j) {
      res[j] = MPReconstructedRatFun(n_vars);

      // numerator
      FFStatus retn = getPoly(n_vars, n_num_terms[j], coefficients, exponents,
                              res[j].numerator());
      if (retn != FF_SUCCESS)
        return 0;
      coefficients += n_num_terms[j];
      exponents += n_num_terms[j]*n_vars;

      // denominator
      FFStatus retd = getPoly(n_vars, n_den_terms[j], coefficients, exponents,
                              res[j].denominator());
      if (retd != FF_SUCCESS)
        return 0;
      coefficients += n_den_terms[j];
      exponents += n_den_terms[j]*n_vars;
    }

    FFRatFunList * ret = new FFRatFunList();
    ret->n_functions = n_functions;
    ret->n_vars = n_vars;
    ret->rf = std::move(res);

    return ret;
  }


  FFNode ffAlgRatFunEval(FFGraph graph, FFNode in_node,
                         const FFRatFunList * rf)
  {
    typedef AnalyticFunctionData Data;
    std::unique_ptr<AnalyticFunction> algptr(new AnalyticFunction());
    std::unique_ptr<Data> data(new Data());
    AnalyticFunction & alg = *algptr;

    unsigned nparsin = rf->n_vars;
    unsigned nfunctions = rf->n_functions;
    data->f.reset(new HornerRatFunPtr[nfunctions]);
    alg.fmap.reset(new MPHornerRatFunMap[nfunctions]);
    for (int i=0; i<nfunctions; ++i)
      get_horner_ratfun(rf, i, data->f[i], alg.fmap[i],
                        session.main_context()->ww);

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);

    alg.init(nparsin, nfunctions, *data);
    unsigned id = g->new_node(std::move(algptr), std::move(data), &in_node);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }

  FFNode ffAlgRatNumEval(FFGraph graph, FFCStr * nums, unsigned n_nums)
  {
    typedef EvalRationalNumbersData Data;
    std::unique_ptr<EvalRationalNumbers> algptr(new EvalRationalNumbers());
    std::unique_ptr<Data> data(new Data());
    EvalRationalNumbers & alg = *algptr;

    std::vector<MPRational> numbers(n_nums);
    for (int i=0; i<n_nums; ++i)
      numbers[i] = MPRational(nums[i]);

    if (!session.graph_exists(graph))
      return FF_ERROR;

    Graph * g = session.graph(graph);
    alg.init(std::move(numbers), *data);
    unsigned id = g->new_node(std::move(algptr), std::move(data), nullptr);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }

  FFUInt * ffEvaluateRatFunList(const FFRatFunList * rf,
                                const FFUInt * input, unsigned prime_no)
  {
    if (prime_no >= BIG_UINT_PRIMES_SIZE)
      return 0;

    Mod mod(BIG_UINT_PRIMES[prime_no]);
    unsigned nparsout = rf->n_functions;

    FFUInt * output = (FFUInt*)malloc(sizeof(FFUInt)*nparsout);

    for (unsigned j=0; j<nparsout; ++j) {
      output[j] = rf->rf[j].eval(input, mod);
      if (output[j] == FAILED) {
        free(output);
        return 0;
      }
    }

    return output;
  }


  // Reconstruction

  FFStatus ffReconstructFunction(FFGraph graph, FFRecOptions options,
                                 FFRatFunList ** results)
  {
    if (!options.min_primes)
      options.min_primes = 1;

    if (!options.max_primes)
      options.max_primes = 1;

    if (!options.max_deg)
      options.max_deg = RatFunReconstruction::DEFAULT_MAX_DEG;

    ReconstructionOptions opt = toRecOpt(options);

    Graph * g = session.graph(graph);

    if (!g)
      return FF_ERROR;

    typedef MPReconstructedRatFun ResT;
    const unsigned nparsout = g->nparsout;
    std::unique_ptr<ResT[]> res (new ResT[nparsout]);

    Ret ret = session.full_reconstruction(graph, res.get(),
                                          options.n_threads, opt);

    // check success
    switch (ret) {

    case SUCCESS:
      {
        FFRatFunList * ret = new FFRatFunList();
        ret->n_functions = nparsout;
        ret->n_vars = g->nparsin[0];
        ret->rf = std::move(res);
        *results = ret;
        return FF_SUCCESS;
      }

    case MISSING_SAMPLES:
      return FF_MISSING_POINTS;

    case MISSING_PRIMES:
      return FF_MISSING_PRIMES;

    default:
      return FF_ERROR;

    }
  }

} // extern "C"
