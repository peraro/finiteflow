#include <algorithm>
#include <cstring>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <fflow/common.hh>
#include <fflow/primes.hh>
#include <fflow/capi.h>
#include <fflow/graph.hh>
#include <fflow/subgraph.hh>
#include <fflow/subgraph_reconstruct.hh>
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
#include <fflow/ratexpr_parser.hh>
#include <fflow/eval_count.hh>
#include <fflow/algorithm.hh>
#include <fflow/mp_gcd.hh>
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

static ReconstructionOptions toRecOpt(FFRecOptions & options)
{
  if (!options.min_primes)
    options.min_primes = 1;

  if (!options.max_primes)
    options.max_primes = 1;

  if (!options.max_deg)
    options.max_deg = RatFunReconstruction::DEFAULT_MAX_DEG;

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
    if (int(elems[j]) <= last_var || elems[j] >= n_vars + 1) {
      logerr("Invalid list of columns in node sparse solver");
      return FF_ERROR;
    }
    last_var = elems[j];
    *dest = elems[j];
  }
  return FF_SUCCESS;
}
static FFStatus CopyNonZeroIdxEls(std::size_t max_els,
                                  const std::size_t * elems,
                                  std::size_t n_elems,
                                  std::size_t * dest)
{
  for (unsigned j=0; j<n_elems; ++j, ++dest) {
    if (elems[j] >= max_els) {
      logerr("Indexes of non-zero coefficients out of bounds");
      return FF_ERROR;
    }
    *dest = elems[j];
  }
  return FF_SUCCESS;
}

static FFStatus CopyNonZeroIdxElsEx(const unsigned * nparsin,
                                    unsigned n_input_nodes,
                                    std::size_t max_els,
                                    const std::size_t * elems,
                                    const unsigned * weights,
                                    std::size_t n_elems,
                                    AnalyticSparseSolverEx::Weight * dest)
{
  for (unsigned j=0; j<n_elems; ++j, ++dest, weights += 2) {
    if (elems[j] >= max_els) {
      logerr("Indexes of non-zero coefficients out of bounds");
      return FF_ERROR;
    }
    if (weights[0]+1 >= n_input_nodes || weights[1] >= nparsin[weights[0]+1]) {
      logerr("Indexes of weights out of bounds");
      return FF_ERROR;
    }
    dest->idx = elems[j];
    dest->node = weights[0];
    dest->el = weights[1];
  }
  return FF_SUCCESS;
}

static FFStatus CopyNonZeroColumnElsUnordered(unsigned n_cols,
                                              const unsigned * elems,
                                              unsigned n_elems,
                                              unsigned * dest)
{
  for (unsigned j=0; j<n_elems; ++j, ++dest) {
    if (elems[j] >= n_cols) {
      logerr("Indexes of non-zero coefficients out of bounds");
      return FF_ERROR;
    }
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

    if (!(coeffs[k].set(*coefficients) == 0)) {
      logerr("Could not convert string to rational");
      return FF_ERROR;
    }

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

static char * writer_to_cstr(MemoryWriter & w)
{
  char * ret = (char*)malloc(w.size() + 1);
  strcpy(ret, w.c_str());
  return ret;
}

static char * cstr_to_cstr(const char * str, size_t n)
{
  char * ret = (char*)malloc(n + 1);
  strcpy(ret, str);
  return ret;
}

template <typename TPoly>
static uint16_t * ffRatNumRatPolyExps(const TPoly & poly)
{
  unsigned num_size = poly.size();
  unsigned nv = poly.nvars();

  uint16_t * ret = (uint16_t *)malloc(num_size * nv * sizeof(uint16_t));
  uint16_t * out = ret;
  for (unsigned j=0; j<num_size; ++j)
    for (unsigned k=0; k<nv; ++k, ++out)
      *out = poly.monomial(j).exponent(k);

  return ret;
}


extern "C" {

  struct FFRatFunList {
    std::unique_ptr<MPReconstructedRatFun[]> rf;
    size_t n_functions = 0;
    unsigned n_vars = 0;
  };

  struct FFIdxRatFunList {
    FFRatFunList rf;
    std::unique_ptr<std::size_t[]> idx;
    size_t n_indexes = 0;
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

  Ret get_horner_ratfun_coeffs(const FFRatFunList * rf,
                               unsigned ncoeffs,
                               unsigned idx,
                               HornerRatFunPtr & f,
                               CoeffHornerRatFunMap & map,
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
    for (unsigned j=0; j<nterms_num; ++j) {
      Int n, d;
      num.coeff(j).to_int(n, d);
      if (d != 1 || n < 0 || n >= ncoeffs) {
        logerr("RatFunEvalFromCoeffs coefficient index out of bounds");
        return FAILED;
      }
      num_c[j] = n;
    }


    // denominator
    map.den_map.resize(nterms_den);
    auto den_c = map.den_map.coeff.get();
    for (unsigned j=0; j<nterms_den; ++j) {
      Int n, d;
      den.coeff(j).to_int(n, d);
      if (d != 1 || n < 0 || n >= ncoeffs) {
        logerr("RatFunEvalFromCoeffs coefficient index out of bounds");
        return FAILED;
      }
      den_c[j] = n;
    }

    // build result
    f.from_sparse_poly(num.monomials(), num.size(),
                       den.monomials(), den.size(),
                       rf->n_vars, 0,
                       map.num_map.pos.get(), map.den_map.pos.get());

    std::size_t max_workspace = horner_required_workspace(f.num());
    max_workspace = std::max(horner_required_workspace(f.den()),
                             max_workspace);
    workspace.ensure_size(max_workspace);

    return SUCCESS;
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

  unsigned ffVersion(void)
  {
    return FFLOW_VERSION;
  }

  unsigned ffVersionMinor(void)
  {
    return FFLOW_VERSION_MINOR;
  }

  unsigned ffNAvailablePrimes(void)
  {
    return BIG_UINT_PRIMES_SIZE;
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

  void ffFreeMemoryU16(uint16_t * mem)
  {
    free(mem);
  }

  void ffFreeCStrArray(char ** mem)
  {
    if (!mem)
      return;
    char ** el = mem;
    for (; *el != 0; ++el)
      free(*el);
    free(mem);
  }

  void ffFreeCStr(char * mem)
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

  FFStatus ffGetOutputNode(FFGraph graph)
  {
    FFNode node = session.get_output_node(graph);
    if (node == ALG_NO_ID)
      return FF_ERROR;
    return node;
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

  FFStatus ffLearnEx(FFGraph graph, FFLearnOptions opt)
  {
    if (opt.prime_no) {
      FFNode node = ffGetOutputNode(graph);
      if (!session.node_exists(graph, node))
        return FF_ERROR;
      LearningOptions & lopt = session.node(graph, node)->learn_opt;
      lopt.prime_no = opt.prime_no;
    }
    return ffLearn(graph);
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

  FFUInt * ffEvaluatePoints(FFGraph graph,
                            const FFUInt * input, unsigned n_points,
                            unsigned n_threads)
  {
    Graph * g = session.graph(graph);
    if (!session.graph_can_be_evaluated(graph))
      return 0;

    unsigned nparsin = g->nparsin[0];
    unsigned nparsout = g->nparsout;

    SamplePointsVector xin(n_points);
    SamplePointsVector xout;

    for (unsigned j=0; j<n_points; ++j, input+=nparsin+1) {
      xin[j].reset(new UInt[nparsin+1]);
      std::copy(input, input+nparsin, xin[j].get());
      UInt prime_no = input[nparsin];
      if (prime_no >= BIG_UINT_PRIMES_SIZE) {
        logerr("Prime index out of bounds");
        return 0;
      }
      xin[j][nparsin] = BIG_UINT_PRIMES[prime_no];
    }

    Ret ret = session.evaluate_list(graph, xin, xout, n_threads);
    if (ret != SUCCESS)
      return 0;

    FFUInt * output = (FFUInt*)malloc(sizeof(FFUInt)*nparsout*n_points);
    FFUInt * dest = output;
    for (unsigned j=0; j<n_points; ++j, dest += nparsout)
      std::copy(xout[j].get(), xout[j].get()+nparsout, dest);

    return output;
  }

  unsigned ffSubgraphNParsout(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg)
      return FF_ERROR;

    SubGraph * subalg = dynamic_cast<SubGraph*>(alg);
    if (!subalg)
      return FF_ERROR;

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

    return id;
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

    return id;
  }

  FFNode ffAlgSubgraphMap(FFGraph graph,
                          const FFNode * in_nodes, unsigned n_in_nodes,
                          FFGraph subgraph)
  {
    typedef SubGraphData Data;
    std::unique_ptr<SubGraphMap> algptr(new SubGraphMap());
    std::unique_ptr<Data> data(new Data());

    std::vector<unsigned> nparsin(n_in_nodes);
    for (unsigned i=0; i<n_in_nodes; ++i) {
      Node * node = session.node(graph, in_nodes[i]);
      if (!node) {
        logerr("Input node does not exist");
        return FF_ERROR;
      }
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

    return id;
  }

  static FFNode ffAlgSubgraphUniRec(FFGraph graph, FFNode in_node,
                                    FFGraph subgraph)
  {
    typedef SubgraphUniRecData Data;
    std::unique_ptr<SubgraphUniRec> algptr(new SubgraphUniRec());
    std::unique_ptr<Data> data(new Data());

    Node * node = session.node(graph, in_node);
    if (!node) {
      logerr("Invalid input node");
      return FF_ERROR;
    }

    unsigned npars = node->algorithm()->nparsout;
    Ret ret = algptr->init(session, subgraph, *data, &npars, 1);
    if (ret != SUCCESS)
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(data), &in_node);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }

  FFNode ffAlgSubgraphRec(FFGraph graph, FFNode in_node,
                          FFGraph subgraph,
                          unsigned n_rec_vars, bool shift_vars)
  {
    if (n_rec_vars == 1)
      return ffAlgSubgraphUniRec(graph, in_node, subgraph);

    typedef SubgraphRecData Data;
    std::unique_ptr<SubgraphRec> algptr(new SubgraphRec());
    std::unique_ptr<Data> data(new Data());

    Node * node = session.node(graph, in_node);
    if (!node) {
      logerr("Invalid input node");
      return FF_ERROR;
    }

    unsigned npars = node->algorithm()->nparsout;
    Ret ret = algptr->init(session, subgraph, *data,
                           npars, n_rec_vars, shift_vars);
    if (ret != SUCCESS)
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), std::move(data), &in_node);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
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
                                   const size_t * non_zero_coeffs,
                                   const FFRatFunList * rat_functions,
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

    const std::size_t max_functions = rat_functions->n_functions;
    data.c.resize(max_functions);
    sys.cmap.resize(max_functions);

    for (int j=0; j<max_functions; ++j)
      get_horner_ratfun(rat_functions, j, data.c[j], sys.cmap[j],
                        session.main_context()->ww);

    for (int i=0; i<n_rows; ++i) {
      const unsigned csize = n_non_zero[i];
      AnalyticSparseSolver::RowInfo & rinf = sys.rinfo[i];
      rinf.size = csize;
      rinf.cols.reset(new unsigned[csize]);
      rinf.idx.reset(new std::size_t[csize]);
      FFStatus ret = CopyNonZeroColumnEls(n_vars, non_zero_els, csize,
                                          rinf.cols.get());
      if (ret != FF_SUCCESS)
        return ret;
      ret = CopyNonZeroIdxEls(max_functions, non_zero_coeffs, csize,
                              rinf.idx.get());
      non_zero_els += csize;
      non_zero_coeffs += csize;
    }

    sys.nparsin.resize(1);
    sys.nparsin[0] = rat_functions->n_vars;

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

  FFNode ffAlgAnalyticSparseLSolveIdx(FFGraph graph, FFNode in_node,
                                      unsigned n_eqs, unsigned n_vars,
                                      const unsigned * n_non_zero,
                                      const unsigned * non_zero_els,
                                      const FFIdxRatFunList * non_zero_functions,
                                      const unsigned * needed_vars,
                                      unsigned n_needed_vars)
  {
    return ffAlgAnalyticSparseLSolve(graph, in_node, n_eqs, n_vars,
                                     n_non_zero, non_zero_els,
                                     non_zero_functions->idx.get(),
                                     &non_zero_functions->rf,
                                     needed_vars, n_needed_vars);
  }

  FFNode ffAlgNumericSparseLSolve(FFGraph graph,
                                  unsigned n_eqs, unsigned n_vars,
                                  const unsigned * n_non_zero,
                                  const unsigned * non_zero_els,
                                  const size_t * non_zero_coeffs,
                                  FFCStr * rat_coeffs,
                                  size_t n_rat_coeffs,
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
    sys.c.resize(n_rat_coeffs);

    for (int i=0; i<n_rat_coeffs; ++i)
      sys.c[i] = MPRational(rat_coeffs[i]);

    for (int i=0; i<n_rows; ++i) {
      const unsigned csize = n_non_zero[i];
      NumericSparseSolver::RowInfo & rinf = sys.rinfo[i];
      rinf.size = csize;
      rinf.cols.reset(new unsigned[csize]);
      rinf.idx.reset(new std::size_t[csize]);
      FFStatus ret = CopyNonZeroColumnEls(n_vars, non_zero_els, csize,
                                          rinf.cols.get());
      if (ret != FF_SUCCESS)
        return ret;
      ret = CopyNonZeroIdxEls(n_rat_coeffs, non_zero_coeffs, csize,
                              rinf.idx.get());
      non_zero_els += csize;
      non_zero_coeffs += csize;
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

  FFNode ffAlgAnalyticSparseLSolveEx(FFGraph graph,
                                     FFNode * in_nodes, unsigned n_in_nodes,
                                     unsigned n_eqs, unsigned n_vars,
                                     const unsigned * n_non_zero,
                                     const unsigned * non_zero_els,
                                     const unsigned * n_weights,
                                     const unsigned * weights,
                                     const size_t * non_zero_coeffs,
                                     const FFRatFunList * rat_functions,
                                     const unsigned * needed_vars,
                                     unsigned n_needed_vars)
  {
    typedef AnalyticSparseSolverExData Data;
    std::unique_ptr<AnalyticSparseSolverEx> algptr(new AnalyticSparseSolverEx());
    std::unique_ptr<Data> dataptr(new Data());
    auto & sys = *algptr;
    auto & data = *dataptr;

    const unsigned n_rows = n_eqs;
    sys.rinfo.resize(n_rows);

    const std::size_t max_functions = rat_functions->n_functions;
    data.c.resize(max_functions);
    sys.cmap.resize(max_functions);

    for (int j=0; j<max_functions; ++j)
      get_horner_ratfun(rat_functions, j, data.c[j], sys.cmap[j],
                        session.main_context()->ww);

    if (n_in_nodes < 1) {
      logerr("At least one input node is required");
      return FF_ERROR;
    }
    sys.nparsin.resize(n_in_nodes);
    sys.nparsin[0] = rat_functions->n_vars;
    for (unsigned j=1; j<n_in_nodes; ++j) {
      unsigned node_nout = ffNodeNParsOut(graph, in_nodes[j]);
      if (ffIsError(node_nout)) {
        logerr("Invalid input node.");
        return FF_ERROR;
      }
      sys.nparsin[j] = node_nout;
    }

    for (int i=0; i<n_rows; ++i) {

      const unsigned csize = n_non_zero[i];
      AnalyticSparseSolverEx::RowInfo & rinf = sys.rinfo[i];
      rinf.size = csize;
      rinf.cols.reset(new unsigned[csize]);
      FFStatus ret = CopyNonZeroColumnEls(n_vars, non_zero_els, csize,
                                          rinf.cols.get());
      if (ret != FF_SUCCESS)
        return ret;

      auto & ws = rinf.w;
      ws.reset(new AnalyticSparseSolverEx::Weights[csize]);
      for (unsigned j=0; j<csize; ++j) {
        auto & w = ws[j];
        const unsigned n_ws = *n_weights;
        w.size = n_ws;
        w.w.reset(new AnalyticSparseSolverEx::Weight[n_ws]);
        ret = CopyNonZeroIdxElsEx(sys.nparsin.data(), n_in_nodes,
                                  max_functions,
                                  non_zero_coeffs, weights, n_ws,
                                  w.w.get());
        if (ret != FF_SUCCESS)
          return ret;
        non_zero_coeffs += n_ws;
        ++n_weights;
        weights += 2*n_ws;
      }

      non_zero_els += csize;
    }

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
                              in_nodes);
    return id;
  }

  FFNode ffAlgAnalyticSparseLSolveIdxEx(FFGraph graph,
                                        FFNode * in_nodes, unsigned n_in_nodes,
                                        unsigned n_eqs, unsigned n_vars,
                                        const unsigned * n_non_zero,
                                        const unsigned * non_zero_els,
                                        const unsigned * n_weights,
                                        const unsigned * weights,
                                        const FFIdxRatFunList * non_zero_funcs,
                                        const unsigned * needed_vars,
                                        unsigned n_needed_vars)
  {
    return ffAlgAnalyticSparseLSolveEx(graph, in_nodes, n_in_nodes,
                                       n_eqs, n_vars,
                                       n_non_zero, non_zero_els,
                                       n_weights, weights,
                                       non_zero_funcs->idx.get(),
                                       &non_zero_funcs->rf,
                                       needed_vars, n_needed_vars);
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
    }

    Ret ret = alg.init(nparsin.data(), nparsin.size(), std::move(els));

    if (!session.graph_exists(graph) || ret != SUCCESS)
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
      }

    }

    Ret ret = alg.init(nparsin.data(), nparsin.size(), std::move(els));

    if (!session.graph_exists(graph) || ret != SUCCESS)
      return FF_ERROR;

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), nullptr, in_nodes);

    return id;
  }

  FFNode ffAlgTakeAndAddBL(FFGraph graph,
                           const FFNode * in_nodes, unsigned n_in_nodes,
                           unsigned n_elems,
                           const unsigned * elems_len,
                           const unsigned * elems)
  {
    std::unique_ptr<TakeAndAddBL> algptr(new TakeAndAddBL());
    TakeAndAddBL & alg = *algptr;
    std::vector<unsigned> nparsin(n_in_nodes);
    for (unsigned i=0; i<n_in_nodes; ++i) {
      Node * node = session.node(graph, in_nodes[i]);
      if (!node)
        return FF_ERROR;
      nparsin[i] = node->algorithm()->nparsout;
    }

    std::vector<std::vector<TakeAndAddBL::InputEl>> els;
    els.resize(n_elems);
    for (int j=0; j<n_elems; ++j) {

      const unsigned this_len = elems_len[j];
      els[j].resize(this_len);

      for (int k=0; k<this_len; ++k, elems += 4) {
        els[j][k].list1 = elems[0];
        els[j][k].el1 = elems[1];
        els[j][k].list2 = elems[2];
        els[j][k].el2 = elems[3];
      }

    }

    Ret ret = alg.init(nparsin.data(), nparsin.size(), std::move(els));

    if (!session.graph_exists(graph) || ret != SUCCESS)
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

    Ret ret = alg.init(n_rows_a, n_cols_a, n_cols_b);

    if (!session.graph_exists(graph) || ret != SUCCESS)
      return FF_ERROR;

    Graph * g = session.graph(graph);
    FFNode input[2] = {in_node_a, in_node_b};
    unsigned id = g->new_node(std::move(algptr), nullptr, input);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }


  // Sparse systems

  FFStatus ffLSolveNEqsNVars(FFGraph graph, FFNode node, unsigned res[])
  {
    if (!session.node_exists(graph,node))
      return FF_ERROR;

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
      res[0] = ls.neqs();
      res[1] = ls.nvars();
      return FF_SUCCESS;
    } else if(dynamic_cast<SparseLinearSolver*>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
      res[0] = ls.neqs();
      res[1] = ls.nvars();
      return FF_SUCCESS;
    }

    logerr("Not a solver.");
    return FF_ERROR;
  }

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

  FFStatus ffLSolveOnlyNonHomogeneous(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg || !alg->is_mutable())
      return FF_ERROR;


    if (dynamic_cast<DenseLinearSolver *>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver *>(alg);
      if (ls.only_non_homogeneous() != SUCCESS)
        return FF_ERROR;
      session.invalidate_subctxt_alg_data(graph, node);

    } else if (dynamic_cast<SparseLinearSolver *>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
      if (ls.only_non_homogeneous() != SUCCESS)
        return FF_ERROR;
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

  FFStatus ffLSolveSparseOutputWithMaxCol(FFGraph graph, FFNode node,
                                          unsigned max_col,
                                          bool back_substitution,
                                          bool keep_full_output)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg || !alg->is_mutable())
      return FF_ERROR;

    if (dynamic_cast<SparseLinearSolver *>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
      ls.sparse_output_with_maxcol(max_col, back_substitution,
                                   keep_full_output);
      session.invalidate_subctxt_alg_data(graph, node);
    } else {
      return FF_ERROR;
    }

    return FF_SUCCESS;
  }

  FFStatus ffLSolveOutputIsSparse(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg || !alg->is_mutable())
      return FF_ERROR;

    if (dynamic_cast<SparseLinearSolver *>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
      return ls.output_is_sparse();
    }

    return FF_ERROR;
  }

  FFStatus ffLSolveEqWeight(FFGraph graph, FFNode node, const int * eq_weight)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg || !alg->is_mutable())
      return FF_ERROR;

    if (dynamic_cast<SparseLinearSolver *>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
      ls.set_eq_weight(eq_weight);
      session.invalidate_subctxt_alg_data(graph, node);
    } else {
      logerr("Not a sparse solver.");
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

    } else if (dynamic_cast<AnalyticSparseSolverEx *>(alg)) {
      AnalyticSparseSolverEx & ls = *static_cast<AnalyticSparseSolverEx *>(alg);
      ls.delete_unneeded_eqs(session.alg_data(graph, node));
      session.invalidate_subctxt_alg_data(graph, node);

    } else {
      return FF_ERROR;
    }

    return FF_SUCCESS;
  }

  unsigned ffLSolveIsImpossible(FFGraph graph, FFNode node)
  {
    if (!session.node_exists(graph,node))
      return FF_ERROR;

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
      return ls.is_impossible();
    } else if(dynamic_cast<SparseLinearSolver*>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
      return ls.is_impossible();
    }

    return FF_ERROR;
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

  unsigned ffLSolveNIndepEqs(FFGraph graph, FFNode node)
  {
    if (!session.node_exists(graph,node))
      return FF_ERROR;

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
      return ls.n_indep_eqs();
    } else if(dynamic_cast<SparseLinearSolver*>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
      return ls.n_indep_eqs();
    }

    return FF_ERROR;
  }

  unsigned * ffLSolveIndepEqs(FFGraph graph, FFNode node)
  {
    if (!session.node_exists(graph,node))
      return 0;

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
      DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
      return newU32Array(ls.indep_eqs(), ls.n_indep_eqs());
    } else if(dynamic_cast<SparseLinearSolver*>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
      return newU32Array(ls.indep_eqs(), ls.n_indep_eqs());
    }

    return 0;
  }

  FFStatus ffLSolveZeroVars(FFGraph graph, FFNode node,
                            unsigned ** zerovars, unsigned * n_zerovars)
  {
    if (!session.node_exists(graph,node))
      return FF_ERROR;

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<DenseLinearSolver*>(alg)) {

      DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
      std::vector<unsigned> zeroes;
      for (unsigned i=0; i<ls.nvars(); ++i)
        if (!(ls.xinfo()[i] & LSVar::IS_NON_ZERO))
          zeroes.push_back(i);
      if (zerovars)
        *zerovars = newU32Array(zeroes.data(), zeroes.size());
      if (n_zerovars)
        *n_zerovars = zeroes.size();
      return FF_SUCCESS;

    } else if (dynamic_cast<SparseLinearSolver*>(alg)) {

      SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
      std::vector<unsigned> zeroes;
      ls.zero_vars(zeroes);
      if (zerovars)
        *zerovars = newU32Array(zeroes.data(), zeroes.size());
      if (n_zerovars)
        *n_zerovars = zeroes.size();
      return FF_SUCCESS;

    }

    logerr("Not a solver");
    return FF_ERROR;
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

  void ffFreeIdxRatFun(FFIdxRatFunList * rf)
  {
    delete rf;
  }

  unsigned ffIdxRatFunListSize(const FFIdxRatFunList * rf)
  {
    return rf->n_indexes;
  }

  unsigned ffIdxRatFunListNFunctions(const FFIdxRatFunList * rf)
  {
    return rf->rf.n_functions;
  }

  unsigned ffIdxRatFunListNVars(const FFIdxRatFunList * rf)
  {
    return rf->rf.n_vars;
  }



  FFIdxRatFunList * ffMoveRatFunToIdx(FFRatFunList * rf,
                                      const size_t * idx, size_t n_indexes)

  {
    FFIdxRatFunList * ret = new FFIdxRatFunList();
    ret->rf.n_functions = rf->n_functions;
    ret->rf.n_vars = rf->n_vars;
    ret->n_indexes = n_indexes;
    ret->idx.reset(new std::size_t[n_indexes]);

    const size_t max_idx = rf->n_functions;
    for (size_t j=0; j<n_indexes; ++j) {
      if (idx[j] >= max_idx) {
        logerr("Index of IdxRatFun out of bounds");
        goto fail;
      }
      ret->idx[j] = idx[j];
    }

    ret->rf.rf = std::move(rf->rf);
    rf->n_functions = 0;
    return ret;

  fail:
    ffFreeIdxRatFun(ret);
    return 0;
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


  FFRatFunList * ffNewRatFunList(unsigned n_vars, size_t n_functions,
                                 const unsigned * n_num_terms,
                                 const unsigned * n_den_terms,
                                 FFCStr * coefficients,
                                 const uint16_t * exponents)
  {
    typedef MPReconstructedRatFun ResT;
    std::unique_ptr<ResT[]> res (new ResT[n_functions]);

    for (std::size_t j=0; j<n_functions; ++j) {
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

  unsigned ffRatFunNumNTerms(const FFRatFunList * rf, unsigned idx)
  {
    if (idx >= rf->n_functions)
      return FF_ERROR;
    return rf->rf[idx].numerator().size();
  }

  unsigned ffRatFunDenNTerms(const FFRatFunList * rf, unsigned idx)
  {
    if (idx >= rf->n_functions)
      return FF_ERROR;
    return rf->rf[idx].denominator().size();
  }

  static char ** ffRatNumRatPolyCoeff(const MPReconstructedPoly & poly)
  {
    unsigned poly_size = poly.size();

    char ** ccs = (char **)malloc((poly_size+1) * sizeof(char*));
    std::string strbuff;

    MemoryWriter w;
    for (unsigned j=0; j<poly_size; ++j) {
      w.clear();
      poly.coeff(j).print(w);
      ccs[j] = writer_to_cstr(w);
    }
    ccs[poly_size] = 0;

    return ccs;
  }

  char ** ffRatFunNumCoeffs(const FFRatFunList * rf, unsigned idx)
  {
    if (idx >= rf->n_functions)
      return 0;
    const MPReconstructedPoly & poly = rf->rf[idx].numerator();
    return ffRatNumRatPolyCoeff(poly);
  }

  char ** ffRatFunDenCoeffs(const FFRatFunList * rf, unsigned idx)
  {
    if (idx >= rf->n_functions)
      return 0;
    const MPReconstructedPoly & poly = rf->rf[idx].denominator();
    return ffRatNumRatPolyCoeff(poly);
  }

  uint16_t * ffRatFunNumExponents(const FFRatFunList * rf, unsigned idx)
  {
    if (idx >= rf->n_functions)
      return 0;
    const MPReconstructedPoly & poly = rf->rf[idx].numerator();
    return ffRatNumRatPolyExps(poly);
  }

  uint16_t * ffRatFunDenExponents(const FFRatFunList * rf, unsigned idx)
  {
    if (idx >= rf->n_functions)
      return 0;
    const MPReconstructedPoly & poly = rf->rf[idx].denominator();
    return ffRatNumRatPolyExps(poly);
  }

  static void ffMemoryWritePoly(MemoryWriter & w,
                                const MPReconstructedPoly & poly,
                                const FFCStr * vars)
  {
    const unsigned nv = poly.nvars();
    const unsigned poly_size = poly.size();

    if (!poly_size) {
      w << "0";
      return;
    }

    for (unsigned j=0; j<poly_size; ++j) {

      const auto & coeff = poly.coeff(j);
      if (j && coeff.sign() > 0)
        w << "+";
      coeff.print(w);

      const auto & mon = poly.monomial(j);
      for (unsigned k=0; k<nv; ++k) {
        unsigned e = mon.exponent(k);
        if (e) {
          w << "*" << vars[k];
          if (e > 1)
            w << "^" << e;
        }
      }
    }
  }

  char * ffRatFunToStr(const FFRatFunList * rf, unsigned idx,
                       const FFCStr * vars)
  {
    if (idx >= rf->n_functions)
      return 0;

    MemoryWriter w;

    const auto & ratfun = rf->rf[idx];
    const auto & num = ratfun.numerator();

    if (num.size() == 0)
      return cstr_to_cstr("0", 1);

    const auto & den = ratfun.denominator();
    if (den.size() == 1 && den.monomial(0).degree() == 0
        && den.coeff(0).cmp(1) == 0) {
      // trivial denominator
      ffMemoryWritePoly(w,num,vars);
      return writer_to_cstr(w);
    }

    w << "( ";
    ffMemoryWritePoly(w,num,vars);
    w << " )/( ";
    ffMemoryWritePoly(w,den,vars);
    w << " )";
    return writer_to_cstr(w);
  }


  FFIdxRatFunList * ffParseIdxRatFun(FFCStr * vars, unsigned n_vars,
                                     FFCStr * inputs, size_t n_functions,
                                     const size_t * idx,
                                     size_t n_indexes)
  {
    FFRatFunList * rf = ffParseRatFun(vars, n_vars, inputs, n_functions);
    return ffMoveRatFunToIdx(rf, idx, n_indexes);
  }

  FFIdxRatFunList * ffParseIdxRatFunEx(FFCStr * vars, unsigned n_vars,
                                       FFCStr * inputs,
                                       const unsigned * input_strlen,
                                       size_t n_functions,
                                       const size_t * idx,
                                       size_t n_indexes)
  {
    FFRatFunList * rf = ffParseRatFunEx(vars, n_vars, inputs, input_strlen,
                                        n_functions);
    return ffMoveRatFunToIdx(rf, idx, n_indexes);
  }


  FFIdxRatFunList * ffNewIdxRatFunList(unsigned n_vars, size_t n_functions,
                                       const unsigned * n_num_terms,
                                       const unsigned * n_den_terms,
                                       FFCStr * coefficients,
                                       const uint16_t * exponents,
                                       const size_t * idx,
                                       size_t n_indexes)
  {
    FFRatFunList * rf = ffNewRatFunList(n_vars, n_functions,
                                        n_num_terms, n_den_terms,
                                        coefficients, exponents);
    return ffMoveRatFunToIdx(rf, idx, n_indexes);
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

  FFNode ffAlgRatFunEvalFromCoeffs(FFGraph graph,
                                   FFNode coeffs_node,
                                   FFNode vars_node,
                                   const FFRatFunList * rf)
  {
    typedef FunctionFromCoeffsData Data;
    std::unique_ptr<FunctionFromCoeffs> algptr(new FunctionFromCoeffs());
    std::unique_ptr<Data> data(new Data());
    FunctionFromCoeffs & alg = *algptr;

    unsigned nparsin = rf->n_vars;
    unsigned nfunctions = rf->n_functions;

    unsigned ncoeffs = ffNodeNParsOut(graph, coeffs_node);
    if (ffIsError(ncoeffs))
      return FF_ERROR;

    data->f.reset(new HornerRatFunPtr[nfunctions]);
    alg.fmap.reset(new CoeffHornerRatFunMap[nfunctions]);
    Ret ret = 0;
    for (int i=0; i<nfunctions; ++i) {
      ret = get_horner_ratfun_coeffs(rf, ncoeffs, i, data->f[i], alg.fmap[i],
                                     session.main_context()->ww);
      if (ret != SUCCESS)
        return FF_ERROR;
    }

    Graph * g = session.graph(graph);
    if (!g) {
      logerr("Graph does not exist");
      return FF_ERROR;
    }

    FFNode inputnodes[2] = {coeffs_node, vars_node};
    alg.init(ncoeffs, nparsin, nfunctions, *data);
    unsigned id = g->new_node(std::move(algptr), std::move(data),
                              inputnodes);
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

  FFNode ffAlgEvalCount(FFGraph graph, FFNode input)
  {
    Node * node = nullptr;
    if (!(node = session.node(graph, input))) {
      logerr("Input node does not exist");
      return FF_ERROR;
    }

    unsigned nparsin = node->algorithm()->nparsout;
    std::unique_ptr<EvalCount> algptr(new EvalCount());
    EvalCount & alg = *algptr;
    alg.init(nparsin);

    Graph * g = session.graph(graph);
    unsigned id = g->new_node(std::move(algptr), nullptr, &input);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }

  FFUInt ffEvalCountGet(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);

    if (dynamic_cast<EvalCount*>(alg)) {
      EvalCount & c = *static_cast<EvalCount*>(alg);
      return c.getCount();
    }

    logerr("Algorithm is not of type EvalCount");
    return FF_FAILED;
  }

  FFUInt ffEvalCountReset(FFGraph graph, FFNode node, FFUInt count)
  {
    Algorithm * alg = session.algorithm(graph, node);

    if (dynamic_cast<EvalCount*>(alg)) {
      EvalCount & c = *static_cast<EvalCount*>(alg);
      return c.resetCount(count);
    }

    logerr("Algorithm is not of type EvalCount");
    return FF_FAILED;
  }


  // Reconstruction

  typedef enum {
    REC_FULL,
    REC_MOD,
    REC_FULL_FROM_CURRENT_EVALS,
    REC_MOD_FROM_CURRENT_EVALS
  } RecMode_;

  static Ret funrec_implem(unsigned graphid, MPReconstructedRatFun res[],
                           unsigned nthreads, ReconstructionOptions opt,
                           const unsigned * deg_data, bool mod,
                           unsigned min_primes=1)
  {
    if (!session.graph_can_be_evaluated(graphid))
      return FAILED;
    Graph & a = *session.graph(graphid);

    if (a.nparsin[0] < 1)
      return FAILED;

    if (a.nparsin[0] == 1)
      return session.reconstruct_univariate(graphid, res, opt);

    const unsigned max_primes = mod ? 1 : opt.max_primes;

    Ret ret = 0;
    if (!deg_data)
      ret = session.parallel_all_degrees(graphid, nthreads, opt);
    else
      ret = session.set_degrees(graphid, deg_data);

    if (ret != SUCCESS)
      return ret;

    ret = MISSING_PRIMES;
    if (mod)
      min_primes = 1;

    while (ret == MISSING_PRIMES && min_primes<=max_primes) {
      opt.max_primes = min_primes;

      session.parallel_sample(graphid, nthreads, opt);

      if (mod)
        ret = session.parallel_reconstruct_mod(graphid, res, nthreads, opt);
      else
        ret = session.parallel_reconstruct(graphid, res, nthreads, opt);

      if (ret != SUCCESS && ret != MISSING_PRIMES)
        return ret;

      ++min_primes;
    }

    return ret;
  }

  static FFStatus reconstruct_fun(FFGraph graph, FFRecOptions options,
                                  FFRatFunList ** results,
                                  RecMode_ rec_mode,
                                  const unsigned * deg_data = nullptr)
  {
    if (!options.max_primes)
      options.max_primes = DEFAULT_MAX_REC_PRIMES;
    ReconstructionOptions opt = toRecOpt(options);

    Graph * g = session.graph(graph);

    if (!g)
      return FF_ERROR;

    typedef MPReconstructedRatFun ResT;
    const unsigned nparsout = g->nparsout;
    std::unique_ptr<ResT[]> res (new ResT[nparsout]);

    Ret ret = 0;

    switch (rec_mode) {
    case REC_FULL:
      ret = funrec_implem(graph, res.get(), options.n_threads, opt,
                          deg_data, false);
      break;

    case REC_MOD:
      ret = funrec_implem(graph, res.get(), options.n_threads, opt,
                          deg_data, true);
      break;

    case REC_FULL_FROM_CURRENT_EVALS:
      ret = session.parallel_reconstruct(graph, res.get(),
                                         options.n_threads, opt);
      break;

    case REC_MOD_FROM_CURRENT_EVALS:
      ret = session.parallel_reconstruct_mod(graph, res.get(),
                                             options.n_threads, opt);
      break;
    }

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

  FFStatus ffReconstructFunction(FFGraph graph, FFRecOptions options,
                                 FFRatFunList ** results)
  {
    return reconstruct_fun(graph, options, results, REC_FULL);
  }

  FFStatus ffReconstructFunctionMod(FFGraph graph, FFRecOptions options,
                                    FFRatFunList ** results)
  {
    return reconstruct_fun(graph, options, results, REC_MOD);
  }

  FFStatus ffReconstructFromCurrentEvaluations(FFGraph graph,
                                               FFRecOptions options,
                                               FFRatFunList ** results)
  {
    return reconstruct_fun(graph, options, results, REC_FULL_FROM_CURRENT_EVALS);
  }

  FFStatus ffReconstructFromCurrentEvaluationsMod(FFGraph graph,
                                                  FFRecOptions options,
                                                  FFRatFunList ** results)
  {
    return reconstruct_fun(graph, options, results, REC_MOD_FROM_CURRENT_EVALS);
  }

  FFStatus ffReconstructFunctionWithDegrees(FFGraph graph,
                                            FFRecOptions options,
                                            const unsigned * deg_data,
                                            FFRatFunList ** results)
  {
    return reconstruct_fun(graph, options, results, REC_FULL, deg_data);
  }

  FFStatus ffReconstructFunctionWithDegreesMod(FFGraph graph,
                                               FFRecOptions options,
                                               const unsigned * deg_data,
                                               FFRatFunList ** results)
  {
    return reconstruct_fun(graph, options, results, REC_MOD, deg_data);
  }

  unsigned * ffAllDegrees(FFGraph graph, FFRecOptions options)
  {
    ReconstructionOptions opt = toRecOpt(options);

    Ret ret = session.parallel_all_degrees(graph, options.n_threads, opt);
    if (ret != SUCCESS)
      return 0;

    // If it hasn't failed the graph is valid and degrees are available
    const Graph * g = session.graph(graph);

    const unsigned * numdeg = g->degs_data().numdeg.get();
    const unsigned * dendeg = g->degs_data().dendeg.get();
    unsigned nparsout = g->nparsout;

    unsigned * degs = (unsigned*)malloc(2*nparsout*sizeof(unsigned));
    for (unsigned j=0; j<nparsout; ++j) {
      degs[2*j] = numdeg[j];
      degs[2*j+1] = dendeg[j];
    }

    return degs;
  }

  FFStatus ffDumpDegrees(FFGraph graph, FFCStr filename)
  {
    Ret ret = session.dump_degrees(graph, filename);

    if (ret != SUCCESS)
      return FF_ERROR;

    return FF_SUCCESS;
  }

  FFStatus ffNParsFromDegreeFile(FFCStr filename,
                                 unsigned * nparsin, unsigned * nparsout)
  {
    Ret ret = algorithm_npars_from_degree_info(filename, *nparsin, *nparsout);
    if (ret != SUCCESS)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFStatus ffLoadDegrees(FFGraph graph, FFCStr filename)
  {
    Ret ret = session.load_degrees(graph, filename);
    if (ret != SUCCESS)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFStatus ffSetDegrees(FFGraph graph, const unsigned * degdata)
  {
    Ret ret = session.set_degrees(graph, degdata);
    if (ret != SUCCESS)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFStatus ffLoadEvaluations(FFGraph graph, FFCStr * files, unsigned n_files)
  {
    Ret ret = session.load_evaluations(graph, files, n_files);
    if (ret != SUCCESS)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFStatus ffDumpSamplePoints(FFGraph graph, FFCStr filename,
                              FFRecOptions options)
  {
    ReconstructionOptions opt = toRecOpt(options);
    Ret ret = session.dump_sample_points(graph, opt, filename);
    if (ret != SUCCESS)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFUInt ffNSamplePointsInFile(FFCStr filename)
  {
    return samples_file_size(filename);
  }

  FFStatus ffEvaluatePointsInFile(FFGraph graph, FFCStr file,
                                  unsigned start, unsigned npoints,
                                  unsigned nthreads)
  {
    SamplePointsFromFile pts(file, start, npoints);
    Ret ret = session.parallel_sample(graph, nthreads,
                                      ReconstructionOptions(), &pts);

    if (ret == SUCCESS)
      return FF_SUCCESS;
    else
      return FF_ERROR;
  }

  FFStatus ffDumpEvaluations(FFGraph graph, FFCStr filename)
  {
    Ret ret = session.dump_evaluations(graph, filename);
    if (ret != SUCCESS)
      return FF_ERROR;
    return FF_SUCCESS;
  }

  FFStatus ffSample(FFGraph graph, FFRecOptions options)
  {
    unsigned n_threads = options.n_threads;
    ReconstructionOptions opt = toRecOpt(options);

    Ret ret = session.parallel_sample(graph, n_threads, opt);

    if (ret == SUCCESS)
      return FF_SUCCESS;
    else
      return FF_ERROR;
  }

  static FFStatus ffParallelReconstructUnivariate_(FFGraph graph,
                                                   FFRecOptions options,
                                                   FFRatFunList ** results,
                                                   bool mod)
  {
    if (!options.max_primes)
      options.max_primes = DEFAULT_MAX_REC_PRIMES;
    if (mod)
      options.max_primes = 1;
    if (!options.n_threads)
      options.n_threads = Session::default_nthreads();
    if (!options.min_deg)
      options.min_deg = std::max<int>(int(options.n_threads/2)-3,2);
    if (!options.max_deg)
      options.max_deg = std::max<int>(RatFunReconstruction::DEFAULT_MAX_DEG,
                                      options.min_deg);
    if (!options.deg_step)
      options.deg_step = std::max<int>(options.n_threads/2,1);

    unsigned np = options.min_primes;
    unsigned ret = FF_MISSING_POINTS;
    unsigned deg = options.min_deg;

    unsigned max_primes = options.max_primes;
    unsigned max_degree = options.max_deg;

    while (true) {

      if (np > max_primes || deg > max_degree)
        break;

      options.max_deg = deg;
      options.max_primes = np;

      if (ffSample(graph, options) != FF_SUCCESS) {
        ret = FF_ERROR;
        break;
      }

      if (!mod)
        ret = ffReconstructFromCurrentEvaluations(graph, options, results);
      else
        ret = ffReconstructFromCurrentEvaluationsMod(graph, options, results);

      if (ret == FF_SUCCESS) {
        break;
      } else if (!mod && ret == FF_MISSING_PRIMES) {
        ++np;
      } else if (deg < max_degree) {
        deg += options.deg_step;
        if (deg > max_degree)
          deg = max_degree;
      } else {
        deg += options.deg_step;
      }

    }

    return ret;
  }

  FFStatus ffParallelReconstructUnivariate(FFGraph graph,
                                           FFRecOptions options,
                                           FFRatFunList ** results)
  {
    return ffParallelReconstructUnivariate_(graph, options, results, false);
  }

  FFStatus ffParallelReconstructUnivariateMod(FFGraph graph,
                                              FFRecOptions options,
                                              FFRatFunList ** results)
  {
    return ffParallelReconstructUnivariate_(graph, options, results, true);
  }

  char ** ffReconstructNumeric(FFGraph graph, FFRecOptions options)
  {
    if (!options.max_primes)
      options.max_primes = DEFAULT_MAX_REC_PRIMES;
    ReconstructionOptions opt = toRecOpt(options);

    Graph * g = session.graph(graph);

    if (!g)
      return 0;

    typedef MPRational ResT;
    const unsigned nparsout = g->nparsout;
    std::unique_ptr<ResT[]> res(new ResT[nparsout]);

    Ret ret = session.reconstruct_numeric(graph, res.get(), opt);

    if (ret != SUCCESS)
      return 0;

    char ** ccs = (char **)malloc((nparsout+1) * sizeof(char*));
    std::string strbuff;

    MemoryWriter w;
    for (unsigned j=0; j<nparsout; ++j) {
      w.clear();
      res[j].print(w);
      ccs[j] = writer_to_cstr(w);
    }
    ccs[nparsout] = 0;

    return ccs;
  }


  // SubgraphRec

  static bool checkSubgraphRecIdx(const SubGraph * alg, unsigned idx)
  {
    if (idx >= alg->subgraph()->nparsout) {
      logerr("Function index of SubgraphRec out of bounds");
      return false;
    }
    return true;
  }

  unsigned ffSubgraphRecNVars(FFGraph graph, FFNode node)
  {
    Algorithm * alg = session.algorithm(graph, node);

    if (dynamic_cast<SubgraphRec*>(alg))
      return static_cast<SubgraphRec*>(alg)->n_rec_vars();

    if (dynamic_cast<SubgraphUniRec*>(alg))
      return static_cast<SubgraphUniRec*>(alg)->n_rec_vars();

    logerr("Not a SubgrahRec algorithm");
    return FF_ERROR;
  }

  static const SparseRationalFunction *
  getSubgrahRecFun(FFGraph graph, FFNode node, unsigned idx)
  {
    Algorithm * alg = session.algorithm(graph, node);
    if (!alg) {
      logerr("Node does not exist");
      return 0;
    }
    if (!alg->has_learned()) {
      logerr("Learning not complete");
      return 0;
    }

    if (dynamic_cast<SubgraphRec*>(alg)) {
      SubgraphRec * subrec = static_cast<SubgraphRec*>(alg);
      if (!checkSubgraphRecIdx(subrec,idx))
        return 0;
      return &subrec->rec_function()[idx];
    }

    if (dynamic_cast<SubgraphUniRec*>(alg)) {
      SubgraphUniRec * subrec = static_cast<SubgraphUniRec*>(alg);
      if (!checkSubgraphRecIdx(subrec,idx))
        return 0;
      return &subrec->rec_function()[idx];
    }

    logerr("Not a SubgrahRec algorithm");
    return 0;
  }

  unsigned ffSubgraphRecNumNTerms(FFGraph graph, FFNode node, unsigned idx)
  {
    const SparseRationalFunction * rf = getSubgrahRecFun(graph,node,idx);
    if (!rf)
      return FF_ERROR;
    return rf->numerator().size();
  }

  unsigned ffSubgraphRecDenNTerms(FFGraph graph, FFNode node, unsigned idx)
  {
    const SparseRationalFunction * rf = getSubgrahRecFun(graph,node,idx);
    if (!rf)
      return FF_ERROR;
    return rf->denominator().size();
  }

  uint16_t * ffSubgraphRecNumExponents(FFGraph graph, FFNode node,
                                       unsigned idx)
  {
    const SparseRationalFunction * rf = getSubgrahRecFun(graph,node,idx);
    if (!rf)
      return 0;
    return ffRatNumRatPolyExps(rf->numerator());
  }

  uint16_t * ffSubgraphRecDenExponents(FFGraph graph, FFNode node,
                                       unsigned idx)
  {
    const SparseRationalFunction * rf = getSubgrahRecFun(graph,node,idx);
    if (!rf)
      return 0;
    return ffRatNumRatPolyExps(rf->denominator());
  }


  // Chinese Remainder and RatRec

  char ** ffChineseRemainder(const FFCStr * z1, FFCStr mod1,
                             const FFUInt * z2, FFUInt mod2,
                             unsigned len)
  {
    MPInt p1(mod1);
    MPInt c1, c2, p12;
    chinese_remainder_coeffs(p1, mod2, c1, c2, p12);

    char ** out = (char**)malloc(sizeof(char*)*(len+2));

    MemoryWriter w;
    for (unsigned j=0; j<len; ++j) {
      MPInt zout(z1[j]);
      chinese_remainder_from_coeffs(zout, z2[j], c1, c2, p12, zout);
      w.clear();
      zout.print(w);
      out[j] = writer_to_cstr(w);
    }

    // write total mod in last element
    w.clear();
    p12.print(w);
    out[len] = writer_to_cstr(w);

    out[len+1] = 0;

    return out;
  }

  char ** ffChineseRemainderCoeffs(FFCStr mod1, FFUInt mod2)
  {
    MPInt p1(mod1);
    MPInt c1, c2, p12;
    chinese_remainder_coeffs(p1, mod2, c1, c2, p12);

    char ** out = (char**)malloc(sizeof(char*)*(3+1));

    MemoryWriter w;
    c1.print(w);
    out[0] = writer_to_cstr(w);
    w.clear();
    c2.print(w);
    out[1] = writer_to_cstr(w);
    w.clear();
    p12.print(w);
    out[2] = writer_to_cstr(w);

    out[3] = 0;

    return out;
  }

  char ** ffChineseRemainderFromCoeffs(const FFCStr * z1, const FFUInt * z2,
                                       unsigned len,
                                       FFCStr c1in, FFCStr c2in, FFCStr mod12)
  {
    MPInt c1(c1in), c2(c2in), p12(mod12);

    char ** out = (char**)malloc(sizeof(char*)*(len+1));

    MemoryWriter w;
    for (unsigned j=0; j<len; ++j) {
      MPInt zout(z1[j]);
      chinese_remainder_from_coeffs(zout, z2[j], c1, c2, p12, zout);
      w.clear();
      zout.print(w);
      out[j] = writer_to_cstr(w);
    }

    out[len] = 0;

    return out;
  }

  char ** ffRatRec(const FFCStr * z1, FFCStr mod, unsigned len)
  {
    MPInt p(mod);
    char ** out = (char**)malloc(sizeof(char*)*(len+1));

    MemoryWriter w;
    for (unsigned j=0; j<len; ++j) {
      MPInt z(z1[j]);
      MPRational q;
      rat_rec(z, p, q);
      w.clear();
      q.print(w);
      out[j] = writer_to_cstr(w);
    }
    out[len] = 0;

    return out;
  }

  char ** ffParallelRatRec(const FFCStr * z1, FFCStr mod, unsigned len,
                           unsigned n_threads)
  {
    MPInt p(mod);
    std::vector<MPRational> q;

    {
      std::vector<MPInt> z(len);
      for (unsigned j=0; j<len; ++j)
        z[j] = MPInt(z1[j]);

      q.resize(len);
      session.parallel_rat_rec(z.data(), len, p, q.data(), n_threads);
    }

    char ** out = (char**)malloc(sizeof(char*)*(len+1));

    MemoryWriter w;
    for (unsigned j=0; j<len; ++j) {
      w.clear();
      q[j].print(w);
      out[j] = writer_to_cstr(w);
    }
    out[len] = 0;

    return out;
  }

  FFNode ffAlgRatExprEvalEx(FFGraph graph, FFNode in_node,
                            FFCStr * vars, unsigned n_vars,
                            FFCStr var_prefix,
                            FFCStr * functions,
                            const unsigned * functions_len,
                            unsigned n_functions)
  {
    if (!vars && !var_prefix) {
      logerr("Either a list of variables of a variable prefix "
             "needs to be specified.");
      return FF_ERROR;
    }

    std::unique_ptr<std::string[]> var_names;
    if (vars) {
      var_names.reset(new std::string[n_vars]);
      for (unsigned j=0; j<n_vars; ++j)
        var_names[j] = vars[j];
    }
    std::string varpref;
    if (var_prefix)
      varpref = var_prefix;

    std::vector<std::vector<Instruction>> bytecode(n_functions);
    std::vector<MPRational> numbers;
    std::vector<AnalyticExpression::VarPow> varpows;

    Ret ret = parse_ratexpr_list(n_vars, varpref, var_names.get(),
                                 functions, functions_len, n_functions,
                                 bytecode, numbers, varpows);
    if (ret != SUCCESS)
      return FF_ERROR;

    Graph * g = session.graph(graph);
    if (!g)
      return FF_ERROR;

    typedef AnalyticExpressionData Data;
    std::unique_ptr<AnalyticExpression> algptr(new AnalyticExpression());
    std::unique_ptr<Data> data(new Data());
    AnalyticExpression & alg = *algptr;

    alg.init(n_vars, std::move(bytecode), std::move(numbers),
             std::move(varpows), *data);
    FFNode id = g->new_node(std::move(algptr), std::move(data), &in_node);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }

  FFNode ffAlgRatExprEval(FFGraph graph, FFNode in_node,
                          FFCStr * vars, unsigned n_vars,
                          FFCStr var_prefix,
                          FFCStr * functions, unsigned n_functions)
  {
    std::unique_ptr<unsigned []> functions_len(new unsigned[n_functions]);
    for (unsigned j=0; j<n_functions; ++j)
      functions_len[j] = strlen(functions[j]);
    return ffAlgRatExprEvalEx(graph, in_node, vars, n_vars, var_prefix,
                              functions, functions_len.get(), n_functions);
  }

} // extern "C"
