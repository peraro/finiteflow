/*
 * C-API for interfacing FiniteFlow and other programs/languages.
 *
 * Eventually it will become a simple and stable API for using
 * FiniteFlow from other programs and programming languages.
 *
 * WARNING: It is currently a work in progress and the API is not
 * stable yet.
 */


#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

  // Used for FFStatus
#define FF_SUCCESS (0)
#define FF_ERROR (~(unsigned)(0))
#define FF_MISSING_POINTS (FF_ERROR - 1)
#define FF_MISSING_PRIMES (FF_ERROR - 2)

  // Used for missing FFGraph, FFNode or things returning unsigned
#define FF_NO_ALGORITHM (~(unsigned)(0))

  // Used for errors in functions returning 64-bit unsigned ints
#define FF_FAILED (~(FFUInt)(0))


  /* API begin */

  typedef uint64_t FFUInt;
  typedef unsigned FFGraph;
  typedef unsigned FFNode;
  typedef unsigned FFStatus;
  typedef const char * FFCStr;

  typedef struct FFRatFunList FFRatFunList;

  // Zero-initialize this to pick default options
  typedef struct {
    unsigned start_mod;
    unsigned min_primes;
    unsigned max_primes;
    unsigned max_deg;
    unsigned dbginfo;
    unsigned polymethod;
    unsigned n_threads;
  } FFRecOptions;

  void ffInit(void);
  void ffDeinit(void);

  // checks if a FFStatus, FFGraph or FFNode represents an error
  bool ffIsError(unsigned val);

  FFUInt ffMulInv(FFUInt z, unsigned prime_no);
  FFUInt ffPrimeNo(unsigned i);

  unsigned ffDefaultNThreads(void);

  // Frees arrays of 32-bit integers returned by several routines
  void ffFreeMemoryU32(unsigned * mem);
  void ffFreeMemoryS32(int * mem);
  // Frees arrays of 64-bit integers returned by several routines
  void ffFreeMemoryU64(uint64_t * mem);

  FFGraph ffNewGraph(void);
  FFGraph ffNewGraphWithInput(unsigned nvars, FFNode * node);
  FFGraph ffNewGraphDummy(unsigned n_in, unsigned n_out);

  FFStatus ffDeleteGraph(FFGraph graph);
  FFStatus ffDeleteNode(FFGraph graph, FFNode node);

  FFStatus ffSetOutputNode(FFGraph graph, FFNode node);

  FFNode ffSetGraphInput(FFGraph graph, unsigned n_vars);

  unsigned ffGraphNParsOut(FFGraph graph);
  unsigned ffNodeNParsOut(FFGraph graph, FFNode node);
  FFStatus ffMakeNodeMutable(FFGraph graph, FFNode node);

  FFStatus ffPruneGraph(FFGraph graph);

  FFStatus ffLearn(FFGraph graph);

  // output can be freed with ffFreeMemoryU64
  FFUInt * ffEvaluateGraph(FFGraph graph,
                           const FFUInt * input, unsigned prime_no);

  // Get the output length of a graph from a subgraph node that uses
  // it or FF_ERROR if this is not a subgraph node
  unsigned ffSubgraphNParsout(FFGraph graph, FFNode node);


  ////////////////
  // Algorithms //
  ////////////////

  FFNode ffAlgSimpleSubgraph(FFGraph graph,
                             const FFNode * in_nodes, unsigned n_in_nodes,
                             FFGraph subgraph);
  FFNode ffAlgMemoizedSubgraph(FFGraph graph,
                               const FFNode * in_nodes, unsigned n_in_nodes,
                               FFGraph subgraph);

  /*
   * Solves the system A.x = b with n_eqs equations and n_vars unknown
   * varibles x.  The entries of the matrix a A and the vector b are
   * rational functions of the free parameters returned by in_node.
   *
   * The input are the non vanishing entries of the n_eqs*vars+1
   * dimensional matrix (A|b).  More precisely:
   * - n_non_zero[i] is the number of non-zero entries of row i
   * - n_zero_cols is a list of the columns of the non-zero entries
   *   for each row, stores in a contigous array (one row after the
   *   other)
   * - non_zero_coeffs is a list of rational functions representing
   *   the non-vanishing coefficients in the matrix.
   *
   * Finally needed_vars is a list of variables whose solution will be
   * in the output, if found.  If NULL, then all unknown variables are
   * needed (in which case n_needed_vars is ignored).
   */
  FFNode ffAlgAnalyticSparseLSolve(FFGraph graph, FFNode in_node,
                                   unsigned n_eqs, unsigned n_vars,
                                   const unsigned * n_non_zero,
                                   const unsigned * non_zero_els,
                                   const FFRatFunList * non_zero_coeffs,
                                   const unsigned * needed_vars,
                                   unsigned n_needed_vars);

  /*
   * Same as ffAlgAnalyticSparseLSolve but the entries (A|b) are
   * rational numbers passed as strings in non_zero_coeffs.
   */
  FFNode ffAlgNumericSparseLSolve(FFGraph graph,
                                  unsigned n_eqs, unsigned n_vars,
                                  const unsigned * n_non_zero,
                                  const unsigned * non_zero_els,
                                  FFCStr * non_zero_coeffs,
                                  const unsigned * needed_vars,
                                  unsigned n_needed_vars);

  /*
   * Same as ffAlgAnalyticSparseLSolve but the non-vanishing entries
   * of (A|b) are taken from the input node in_node.
   */
  FFNode ffAlgNodeSparseLSolve(FFGraph graph, FFNode in_node,
                               unsigned n_eqs, unsigned n_vars,
                               const unsigned * n_non_zero,
                               const unsigned * non_zero_els,
                               const unsigned * needed_vars,
                               unsigned n_needed_vars);

  FFNode ffAlgJSONSparseLSolve(FFGraph graph, FFNode in_node,
                               FFCStr json_file);

  FFNode ffAlgJSONRatFunEval(FFGraph graph, FFNode in_node,
                             FFCStr json_file);
  FFNode ffAlgRatFunEval(FFGraph graph, FFNode in_node,
                         const FFRatFunList * rf);
  FFNode ffAlgRatNumEval(FFGraph graph, FFCStr * nums, unsigned n_nums);

  // The array pointed by order must have length ==
  // ffNodeNParsOut(subgraph).  max_deg is the maximum degree in the
  // expansion variable of the functions to be expanded.  If max_deg <
  // 0 a default value will be used.
  FFNode ffAlgLaurent(FFGraph graph, FFNode in_node, FFGraph subgraph,
                      const int * order, int max_deg);

  // Equivalent to
  // ffAlgLaurent(graph,in_node,subgraph,{order,...,order},max_deg)
  FFNode ffAlgLaurentConstOrder(FFGraph graph, FFNode in_node, FFGraph subgraph,
                                int order, int max_deg);

  FFNode ffAlgMatMul(FFGraph graph, FFNode in_node_a, FFNode in_node_b,
                     unsigned n_rows_a, unsigned n_cols_a, unsigned n_cols_b);

  FFNode ffAlgChain(FFGraph graph,
                    const FFNode * in_nodes, unsigned n_in_nodes);

  // This algorithm simply copies some elements of the input nodes to
  // its output.
  //
  // elems is an array of (n_elems * 2) integers.  The j-th output
  // element of this node will be the elems[2*j+1]-th output element
  // of the elems[2*j]-th input node (as listed in in_nodes), for
  // j=0,...,n_elems-1
  FFNode ffAlgTake(FFGraph graph,
                   const FFNode * in_nodes, unsigned n_in_nodes,
                   const unsigned * elems, unsigned n_elems);

  // Slices the input to its [begin, end) subset
  // If end is negative, the slice ends at the end of the input
  FFNode ffAlgSlice(FFGraph graph, FFNode in_node,
                    unsigned begin, int end);

  FFNode ffAlgAdd(FFGraph graph,
                  const FFNode * in_nodes, unsigned n_in_nodes);
  FFNode ffAlgMul(FFGraph graph,
                  const FFNode * in_nodes, unsigned n_in_nodes);

  FFNode ffAlgTakeAndAdd(FFGraph graph,
                         const FFNode * in_nodes, unsigned n_in_nodes,
                         unsigned n_elems, const unsigned * elems_len,
                         const unsigned * elems);

  FFNode ffAlgSparseMatMul(FFGraph graph, FFNode in_node_a, FFNode in_node_b,
                           unsigned n_rows_a, unsigned n_cols_a,
                           unsigned n_cols_b,
                           const unsigned * n_non_zero_a,
                           const unsigned * non_zero_els_a,
                           const unsigned * n_non_zero_b,
                           const unsigned * non_zero_els_b);


  /////////////////////////////
  // Info about learned algs //
  /////////////////////////////

  // Orders beyond which expansions are truncated (including the last
  // one which is kept in).
  // The length of the returned list is ffSubgraphNParsout(graph,node).
  // The returned list can be freed with ffFreeMemoryS32.
  int * ffLaurentMaxOrders(FFGraph graph, FFNode node);

  // Orders at which expansions start.
  // The length of the returned list is ffSubgraphNParsout(graph,node).
  // The returned list can be freed with ffFreeMemoryS32.
  int * ffLaurentMinOrders(FFGraph graph, FFNode node);


  ////////////////////
  // Sparse systems //
  ////////////////////

  FFStatus ffLSolveResetNeededVars(FFGraph graph, FFNode node,
                                   const unsigned * vars, unsigned n_vars);
  FFStatus ffLSolveOnlyHomogeneous(FFGraph graph, FFNode node);
  FFStatus ffLSolveSparseOutput(FFGraph graph, FFNode node, bool sparse);
  FFStatus ffLSolveMarkAndSweepEqs(FFGraph graph, FFNode node);

  // Only for Analytic and Numeric solvers (fails for Node solvers)
  FFStatus ffLSolveDeleteUnneededEqs(FFGraph graph, FFNode node);

  // no. of dependent unknowns
  unsigned ffLSolveNDepVars(FFGraph graph, FFNode node);

  // Lists the dependent unknowns.  The returned array has legth =
  // ffSparseLSolveNDepVars(graph,node) and its memory can be freed
  // using the ffFreeMemoryU32 function.
  unsigned * ffLSolveDepVars(FFGraph graph, FFNode node);

  // no. of dependent unknowns on the r.h.s. of i-th dep. variable (as
  // returned by ffSparseLSolveDepVars).  If output is not sparse, i
  // is ignored.  If output is sparse and non-homogeneous part is
  // requested, it is added to the number of independent variables.
  unsigned ffLSolveNIndepVars(FFGraph graph, FFNode node, unsigned i);

  // List of dependent unknowns on the r.h.s. of i-th dep. variable
  // (as returned by ffSparseLSolveDepVars).  If output is not sparse,
  // i is ignored.  The array can be freed by ffFreeMemoryU32.
  unsigned * ffLSolveIndepVars(FFGraph graph, FFNode node, unsigned i);


  ////////////////////////////////
  // Rational functions, etc... //
  ////////////////////////////////

  unsigned ffRatFunListSize(const FFRatFunList * rf);
  unsigned ffRatFunListNVars(const FFRatFunList * rf);
  void ffFreeRatFun(FFRatFunList * rf);
  FFStatus ffRatFunToJSON(const FFRatFunList * rf, FFCStr file);


  // This uses a simple and limited parser of rational functions:
  // - functions must be collected under common denominator
  // - numerators and denominators must be in expanded form
  // - rational coefficients must be in front of their monomials
  //   (e.g. "1/2 z" is ok but "z/2" is not)
  //
  // If these lmitations are too restrictive, consider using the
  // parser of a proper CAS and then pass the functions to fflow using
  // ffNewRatFunList instead.
  FFRatFunList * ffParseRatFun(FFCStr * vars, unsigned n_vars,
                               FFCStr * inputs, unsigned n_functions);
  FFRatFunList * ffParseRatFunEx(FFCStr * vars, unsigned n_vars,
                                 FFCStr * inputs, const unsigned * input_strlen,
                                 unsigned n_functions);


  // API for creating a list of n_functions rational functions in
  // n_vars variables.
  //
  // n_num_terms[i] and n_den_terms[i] must contain the no. of terms
  // in the numerator and denominator of the i-th function.
  //
  // The total number of terms is therefore
  //
  //    tot_terms = \sum_i (n_num_terms[i] + n_den_terms[i])
  //
  // with i=0,...,n_functions-1.
  //
  // coefficients is an array of tot_terms strings to be parsed as
  // rational numbers.  They are sorted by function index first.
  // Coefficients belonging to the same function must be sorted with
  // numerator coefficients before denominator coefficients.  The
  // relative order of coefficients belonging to the same polynomial
  // (numerator or denominator) is not important.
  //
  // exponents is an array of tot_terms*n_nvars integers representing
  // exponents, namely it is a pointer to tot_terms arrays of n_nvars
  // indexes, stored contigously in memory.  Each n_vars-dimensional
  // array represents the exponents of a single term.  These must be
  // sorted the same way as their respective coefficients.
  FFRatFunList * ffNewRatFunList(unsigned n_vars, unsigned n_functions,
                                 const unsigned * n_num_terms,
                                 const unsigned * n_den_terms,
                                 FFCStr * coefficients,
                                 const uint16_t * exponents);

  // output must be freed using ffFreeMemoryU64
  FFUInt * ffEvaluateRatFunList(const FFRatFunList * rf,
                                const FFUInt * input, unsigned prime_no);


  ////////////////////
  // Reconstruction //
  ////////////////////

  // on success, results will point to a list of rational functions,
  // which can be cleared with ffFreeRatFun
  FFStatus ffReconstructFunction(FFGraph graph, FFRecOptions options,
                                 FFRatFunList ** results);

  /* API end */

#ifdef __cplusplus
}
#endif
