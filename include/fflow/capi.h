/*
 * C-API for interfacing FiniteFlow and other programs/languages.
 *
 * Eventually it will become a simple and stable API for using
 * FiniteFlow from other programs and programming languages.
 *
 * WARNING: It is currently a work in progress and the API is not
 * stable yet.
 */


#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <fflow/config.hh>

#ifdef __cplusplus
extern "C" {
#endif

  // Possible return values of FFStatus.  FF_ERROR is also used on
  // error when FFGraph or FFNode are expected.
#define FF_SUCCESS (0)
#define FF_ERROR (~(unsigned)(0))
#define FF_MISSING_POINTS (FF_ERROR - 1)
#define FF_MISSING_PRIMES (FF_ERROR - 2)

  // Used for errors in functions returning 64-bit unsigned ints
#define FF_FAILED (~(FFUInt)(0))


  /* API begin */

  typedef uint64_t FFUInt;
  typedef unsigned FFGraph;
  typedef unsigned FFNode;
  typedef unsigned FFStatus;
  typedef const char * FFCStr;

  // An opaque object representing a list of rational functions
  typedef struct FFRatFunList FFRatFunList;

  // An opaque object representing an indexed list of rational
  // functions, namely a list of functions and a list of indexes into
  // the list of functions.  The indexes can be used to avoid
  // repeating equal entries in the list.
  //
  // It is used to define matrices for sparse linear systems, which
  // often have many repeated entries.
  //
  // As an example using
  //     functions = [f0,f1,f2]
  //     indexes = [0,0,1,2,2,0]
  // would be effectively equivalent to the list
  //     functions = [f0,f0,f1,f2,f2,f0]
  // but it stores the functions [f0,f1,f2] only once.
  typedef struct FFIdxRatFunList FFIdxRatFunList;

  // Zero-initialize this to pick default options
  typedef struct {
    unsigned start_mod;
    unsigned min_primes;
    unsigned max_primes;
    unsigned max_deg;
    unsigned min_deg;
    unsigned deg_step;
    unsigned dbginfo;
    unsigned polymethod;
    unsigned n_threads;
  } FFRecOptions;

  typedef struct {
    unsigned prime_no;
  } FFLearnOptions;

  void ffInit(void);
  void ffDeinit(void);

  // These return FFLOW_VERSION and FFLOW_VERSION_MINOR respectively
  // (note: in C/C++ one should just use the macros, these functions
  // are meant for interfaces with other programs or languages).
  unsigned ffVersion(void);
  unsigned ffVersionMinor(void);

  // checks if a FFStatus, FFGraph or FFNode represents an error
  bool ffIsError(unsigned val);

  FFUInt ffMulInv(FFUInt z, unsigned prime_no);
  FFUInt ffPrimeNo(unsigned i);
  unsigned ffNAvailablePrimes(void);

  unsigned ffDefaultNThreads(void);

  // Frees arrays of 32-bit integers returned by several routines
  void ffFreeMemoryU32(unsigned * mem);
  void ffFreeMemoryS32(int * mem);
  // Frees arrays of 64-bit integers returned by several routines
  void ffFreeMemoryU64(uint64_t * mem);
  // Frees arrays of 16-bit integers
  void ffFreeMemoryU16(uint16_t * mem);
  // Frees string
  void ffFreeCStr(char * mem);
  // Frees null-terminated arrays of null-terminated strings
  void ffFreeCStrArray(char ** mem);

  FFGraph ffNewGraph(void);
  FFGraph ffNewGraphWithInput(unsigned nvars, FFNode * node);
  FFGraph ffNewGraphDummy(unsigned n_in, unsigned n_out);

  FFStatus ffDeleteGraph(FFGraph graph);
  FFStatus ffDeleteNode(FFGraph graph, FFNode node);

  FFNode ffGetOutputNode(FFGraph graph);
  FFStatus ffSetOutputNode(FFGraph graph, FFNode node);

  FFNode ffSetGraphInput(FFGraph graph, unsigned n_vars);

  unsigned ffGraphNParsOut(FFGraph graph);
  unsigned ffNodeNParsOut(FFGraph graph, FFNode node);
  FFStatus ffMakeNodeMutable(FFGraph graph, FFNode node);

  FFStatus ffPruneGraph(FFGraph graph);

  FFStatus ffLearn(FFGraph graph);
  FFStatus ffLearnEx(FFGraph graph, FFLearnOptions opt);

  // output can be freed with ffFreeMemoryU64
  FFUInt * ffEvaluateGraph(FFGraph graph,
                           const FFUInt * input, unsigned prime_no);

  // Get the output length of a graph from a subgraph node that uses
  // it or FF_ERROR if this is not a subgraph node
  unsigned ffSubgraphNParsout(FFGraph graph, FFNode node);

  // The array `input` must be a flattened `(nparsin+1)*n_points`
  // matrix of integers, with `nparsin` being the number of input
  // parameters of the graph (which can be obtained as
  // ffNodeNParsOut(graph,0)).  Each row of the matrix corresponding
  // to a point. The first `nparsin` entries of each row are the
  // values of the input parameters while the last is the index of the
  // prime.
  // The output is a flattened `nparsout*n_points` array, whose memory
  // must be freed using ffFreeMemoryU64.
  FFUInt * ffEvaluatePoints(FFGraph graph,
                            const FFUInt * input, unsigned n_points,
                            unsigned n_threads);

  // Returns the FFNode id that would be assigned to a new node of the
  // graph, if defined immediately after this function call.  It can
  // be useful in APIs or applications requiring to assign unique
  // identifiers to nodes before their creation.
  FFNode ffPeekNewNodeId(FFGraph graph);


  ////////////////
  // Algorithms //
  ////////////////

  FFNode ffAlgSimpleSubgraph(FFGraph graph,
                             const FFNode * in_nodes, unsigned n_in_nodes,
                             FFGraph subgraph);
  FFNode ffAlgMemoizedSubgraph(FFGraph graph,
                               const FFNode * in_nodes, unsigned n_in_nodes,
                               FFGraph subgraph);
  FFNode ffAlgSubgraphMap(FFGraph graph,
                          const FFNode * in_nodes, unsigned n_in_nodes,
                          FFGraph subgraph);
  FFNode ffAlgSubgraphRec(FFGraph graph, FFNode in_node,
                          FFGraph subgraph,
                          unsigned n_rec_vars, bool shift_vars);

  /*
   * Solves the system A.x = b with n_eqs equations and n_vars unknown
   * variables x.  The entries of the matrix a A and the vector b are
   * rational functions of the free parameters returned by in_node.
   *
   * The inputs are the non vanishing entries of the n_eqs*(n_vars+1)
   * dimensional matrix (A|b).  More precisely:
   * - n_non_zero[i] is the number of non-zero entries of row i
   * - non_zero_els lists the columns with non-zero entries for each
   *   row, stored in a contiguous array (one row after the other)
   * - non_zero_coeffs is a list of indexes in the array of rational
   *   functions rat_functions representing the non-vanishing
   *   coefficients in the matrix.
   *
   * Finally, needed_vars is a list of variables whose solution will
   * be in the output, if found.  If NULL, then all unknown variables
   * are needed (in which case n_needed_vars is ignored).
   */
  FFNode ffAlgAnalyticSparseLSolve(FFGraph graph, FFNode in_node,
                                   unsigned n_eqs, unsigned n_vars,
                                   const unsigned * n_non_zero,
                                   const unsigned * non_zero_els,
                                   const size_t * non_zero_coeffs,
                                   const FFRatFunList * rat_functions,
                                   const unsigned * needed_vars,
                                   unsigned n_needed_vars);

  /*
   * Same as ffAlgAnalyticSparseLSolve, but the indexes in the array
   * of rational functions and the functions themselves are passed as
   * the argument as a FFIdxRatFunList in non_zero_functions.
   */
  FFNode ffAlgAnalyticSparseLSolveIdx(FFGraph graph, FFNode in_node,
                                      unsigned n_eqs, unsigned n_vars,
                                      const unsigned * n_non_zero,
                                      const unsigned * non_zero_els,
                                      const FFIdxRatFunList * non_zero_functions,
                                      const unsigned * needed_vars,
                                      unsigned n_needed_vars);

  /*
   * Same as ffAlgAnalyticSparseLSolve but the entries (A|b) are
   * indexes non_zero_coeffs[] into the array of rational numbers
   * rat_coeffs[] of length n_rat_coeffs.
   */
  FFNode ffAlgNumericSparseLSolve(FFGraph graph,
                                  unsigned n_eqs, unsigned n_vars,
                                  const unsigned * n_non_zero,
                                  const unsigned * non_zero_els,
                                  const size_t * non_zero_coeffs,
                                  FFCStr * rat_coeffs,
                                  size_t n_rat_coeffs,
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

  // Ex version of analytic solver
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
                                     unsigned n_needed_vars);
  FFNode ffAlgAnalyticSparseLSolveIdxEx(FFGraph graph,
                                        FFNode * in_nodes, unsigned n_in_nodes,
                                        unsigned n_eqs, unsigned n_vars,
                                        const unsigned * n_non_zero,
                                        const unsigned * non_zero_els,
                                        const unsigned * n_weights,
                                        const unsigned * weights,
                                        const FFIdxRatFunList * non_zero_funcs,
                                        const unsigned * needed_vars,
                                        unsigned n_needed_vars);

  FFNode ffAlgJSONSparseLSolve(FFGraph graph, FFNode in_node,
                               FFCStr json_file);

  FFNode ffAlgJSONRatFunEval(FFGraph graph, FFNode in_node,
                             FFCStr json_file);
  FFNode ffAlgRatFunEval(FFGraph graph, FFNode in_node,
                         const FFRatFunList * rf);
  FFNode ffAlgRatFunEvalFromCoeffs(FFGraph graph,
                                   FFNode coeffs_node, FFNode vars_node,
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
  FFNode ffAlgTakeAndAddBL(FFGraph graph,
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

  FFNode ffAlgEvalCount(FFGraph graph, FFNode input);
  FFUInt ffEvalCountGet(FFGraph graph, FFNode node);
  FFUInt ffEvalCountReset(FFGraph graph, FFNode node, FFUInt count);


  // At least one between vars and var_prefix must be specified.  The
  // latter allows to refer to the variables in the expression as
  // "`var_prefix``idx`" where `idx` is their index (zero-based).
  FFNode ffAlgRatExprEval(FFGraph graph, FFNode in_node,
                          FFCStr * vars, unsigned n_vars,
                          FFCStr var_prefix,
                          FFCStr * functions, unsigned n_functions);
  FFNode ffAlgRatExprEvalEx(FFGraph graph, FFNode in_node,
                            FFCStr * vars, unsigned n_vars,
                            FFCStr var_prefix,
                            FFCStr * functions,
                            const unsigned * functions_len,
                            unsigned n_functions);


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

  // Information about reconstructed functions in a SubgraphRec node
  // (same meaning and conventions of similarly named ffRatFun*
  // functions below).
  unsigned ffSubgraphRecNVars(FFGraph graph, FFNode node);
  unsigned ffSubgraphRecNumNTerms(FFGraph graph, FFNode node, unsigned idx);
  unsigned ffSubgraphRecDenNTerms(FFGraph graph, FFNode node, unsigned idx);
  uint16_t * ffSubgraphRecNumExponents(FFGraph graph, FFNode node,
                                       unsigned idx);
  uint16_t * ffSubgraphRecDenExponents(FFGraph graph, FFNode node,
                                       unsigned idx);


  ////////////////////
  // Sparse systems //
  ////////////////////

  FFStatus ffLSolveNEqsNVars(FFGraph graph, FFNode node, unsigned res[]);
  FFStatus ffLSolveResetNeededVars(FFGraph graph, FFNode node,
                                   const unsigned * vars, unsigned n_vars);
  FFStatus ffLSolveOnlyHomogeneous(FFGraph graph, FFNode node);
  FFStatus ffLSolveOnlyNonHomogeneous(FFGraph graph, FFNode node);
  FFStatus ffLSolveSparseOutput(FFGraph graph, FFNode node, bool sparse);
  FFStatus ffLSolveSparseOutputWithMaxCol(FFGraph graph, FFNode node,
                                          unsigned max_col,
                                          bool back_substitution,
                                          bool keep_full_output);
  FFStatus ffLSolveEqWeight(FFGraph graph, FFNode node, const int * eq_weight);
  FFStatus ffLSolveMarkAndSweepEqs(FFGraph graph, FFNode node);

  // Only for Analytic and Numeric solvers (fails for Node solvers)
  FFStatus ffLSolveDeleteUnneededEqs(FFGraph graph, FFNode node);

  // Check if system is impossible (the system needs to have completed
  // learning before this call).  Returns 1 is system is impossible, 0
  // if it is not and FF_ERROR if an error occurred.
  unsigned ffLSolveIsImpossible(FFGraph graph, FFNode node);

  // Output: 0 if false, 1 if true, FF_ERROR if an error occurred.
  unsigned ffLSolveOutputIsSparse(FFGraph graph, FFNode node);

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

  // no. of dependent equations
  unsigned ffLSolveNIndepEqs(FFGraph graph, FFNode node);

  // Lists the independent equations.  The returned array has length =
  // ffLSolveNIndepEqs(graph,node) and its memory can be freed using
  // the ffFreeMemoryU32 function.
  unsigned * ffLSolveIndepEqs(FFGraph graph, FFNode node);

  // On success, *zerovars will point to the list of zero variables
  // (including needed and non-needed) while *n_zerovars will be its
  // length.  The allocated array can be freed using ffFreeMemoryU32.
  FFStatus ffLSolveZeroVars(FFGraph graph, FFNode node,
                            unsigned ** zerovars, unsigned * n_zerovars);


  ////////////////////////////////
  // Rational functions, etc... //
  ////////////////////////////////

  unsigned ffRatFunListSize(const FFRatFunList * rf);
  unsigned ffRatFunListNVars(const FFRatFunList * rf);
  void ffFreeRatFun(FFRatFunList * rf);
  FFStatus ffRatFunToJSON(const FFRatFunList * rf, FFCStr file);

  // These return FF_ERROR if idx is out of bounds
  unsigned ffRatFunNumNTerms(const FFRatFunList * rf, unsigned idx);
  unsigned ffRatFunDenNTerms(const FFRatFunList * rf, unsigned idx);

  // These return a null-terminated array of null-terminated strings
  // representing the rational coefficients
  char ** ffRatFunNumCoeffs(const FFRatFunList * rf, unsigned idx);
  char ** ffRatFunDenCoeffs(const FFRatFunList * rf, unsigned idx);

  // These return an nterms*nvars array of exponents (see
  // ffRatFunListNVars, ffRatFunNumNTerms and ffRatFunDenNTerms)
  uint16_t * ffRatFunNumExponents(const FFRatFunList * rf, unsigned idx);
  uint16_t * ffRatFunDenExponents(const FFRatFunList * rf, unsigned idx);

  // The returned string must be freed using ffFreeCStr
  char * ffRatFunToStr(const FFRatFunList * rf, unsigned idx,
                       const FFCStr * vars);

  unsigned ffIdxRatFunListSize(const FFIdxRatFunList * rf);
  unsigned ffIdxRatFunListNFunctions(const FFIdxRatFunList * rf);
  unsigned ffIdxRatFunListNVars(const FFIdxRatFunList * rf);
  void ffFreeIdxRatFun(FFIdxRatFunList * rf);

  // When successful, this invalidates the input rf and creates a
  // `FFIdxRatFunList` from it, using the provided list of indexes.
  FFIdxRatFunList * ffMoveRatFunToIdx(FFRatFunList * rf,
                                      const size_t * idx, size_t n_indexes);


  // This uses a simple and limited parser of rational functions:
  // - functions must be collected under common denominator
  // - numerators and denominators must be in expanded form
  // - rational coefficients must be in front of their monomials,
  //   except for their denominator (e.g. "3/2 z" and "3z/2" are both
  //   ok but "z 3/2" is not)
  // This may also return valid rational functions for some invalid
  // inputs.
  //
  // If these limitations are too restrictive, consider using the
  // parser of a proper CAS and then pass the functions to fflow using
  // ffNewRatFunList instead.
  FFRatFunList * ffParseRatFun(FFCStr * vars, unsigned n_vars,
                               FFCStr * inputs, unsigned n_functions);
  FFRatFunList * ffParseRatFunEx(FFCStr * vars, unsigned n_vars,
                                 FFCStr * inputs, const unsigned * input_strlen,
                                 unsigned n_functions);

  /*
   * Same as ffParseRatFun and ffParseRatFunEx but returns a
   * FFIdxRatFunList based on the provided list of indexes.
   */
  FFIdxRatFunList * ffParseIdxRatFun(FFCStr * vars, unsigned n_vars,
                                     FFCStr * inputs, size_t n_functions,
                                     const size_t * idx,
                                     size_t n_indexes);
  FFIdxRatFunList * ffParseIdxRatFunEx(FFCStr * vars, unsigned n_vars,
                                       FFCStr * inputs,
                                       const unsigned * input_strlen,
                                       size_t n_functions,
                                       const size_t * idx,
                                       size_t n_indexes);


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
  // indexes, stored contiguously in memory.  Each n_vars-dimensional
  // array represents the exponents of a single term.  These must be
  // sorted the same way as their respective coefficients.
  FFRatFunList * ffNewRatFunList(unsigned n_vars, size_t n_functions,
                                 const unsigned * n_num_terms,
                                 const unsigned * n_den_terms,
                                 FFCStr * coefficients,
                                 const uint16_t * exponents);
  FFIdxRatFunList * ffNewIdxRatFunList(unsigned n_vars, size_t n_functions,
                                       const unsigned * n_num_terms,
                                       const unsigned * n_den_terms,
                                       FFCStr * coefficients,
                                       const uint16_t * exponents,
                                       const size_t * idx,
                                       size_t n_indexes);

  // The returned list must be freed using ffFreeMemoryU64
  FFUInt * ffEvaluateRatFunList(const FFRatFunList * rf,
                                const FFUInt * input, unsigned prime_no);


  ////////////////////
  // Reconstruction //
  ////////////////////

  // On success, results will point to a list of rational functions,
  // which can be cleared with ffFreeRatFun
  FFStatus ffReconstructFunction(FFGraph graph, FFRecOptions options,
                                 FFRatFunList ** results);
  FFStatus ffReconstructFunctionMod(FFGraph graph, FFRecOptions options,
                                    FFRatFunList ** results);
  FFStatus ffParallelReconstructUnivariate(FFGraph graph,
                                           FFRecOptions options,
                                           FFRatFunList ** results);
  FFStatus ffParallelReconstructUnivariateMod(FFGraph graph,
                                              FFRecOptions options,
                                              FFRatFunList ** results);

  // On success, it returns a null-terminated list of null-terminated
  // strings representing the reconstructed rational numbers.  The
  // returned list can be freed using ffFreeCStrArray.
  char ** ffReconstructNumeric(FFGraph graph, FFRecOptions options);

  // Tools for custom reconstruction algorithms
  char ** ffChineseRemainder(const FFCStr * z1, FFCStr mod1,
                             const FFUInt * z2, FFUInt mod2,
                             unsigned len);
  char ** ffChineseRemainderCoeffs(FFCStr mod1, FFUInt mod2);
  char ** ffChineseRemainderFromCoeffs(const FFCStr * z1, const FFUInt * z2,
                                       unsigned len,
                                       FFCStr c1in, FFCStr c2in, FFCStr mod12);
  char ** ffRatRec(const FFCStr * z1, FFCStr mod, unsigned len);
  char ** ffParallelRatRec(const FFCStr * z1, FFCStr mod, unsigned len,
                           unsigned n_threads);


  // Routines to split evaluations and reconstruction into several
  // jobs. Note that:
  // - "sample points" are numerical inputs for a graph
  // - "evaluations" are points that have been evaluated
  // - "nthreads=0" picks a default value
  // - ffAllDegrees, on success, returns the array
  //     [numdeg[0],dendeg[0],...,numdeg[n-1],dendeg[n-1]]
  //   with n == ffGraphNParsOut(graph).  The array must be freed with
  //   ffFreeMemoryU32. On failure it returns null.
  unsigned * ffAllDegrees(FFGraph graph, FFRecOptions options);
  FFStatus ffDumpDegrees(FFGraph graph, FFCStr filename);
  FFStatus ffNParsFromDegreeFile(FFCStr filename,
                                 unsigned * nparsin, unsigned * nparsout);
  FFStatus ffLoadDegrees(FFGraph graph, FFCStr filename);
  FFStatus ffSetDegrees(FFGraph graph, const unsigned * degdata); // see below
  FFStatus ffLoadEvaluations(FFGraph graph, FFCStr * files, unsigned n_files);
  FFStatus ffDumpSamplePoints(FFGraph graph, FFCStr filename,
                              FFRecOptions options);
  FFUInt ffNSamplePointsInFile(FFCStr filename);
  FFStatus ffEvaluatePointsInFile(FFGraph graph, FFCStr file,
                                  unsigned start, unsigned npoints,
                                  unsigned nthreads);
  FFStatus ffDumpEvaluations(FFGraph graph, FFCStr filename);
  FFStatus ffReconstructFromCurrentEvaluations(FFGraph graph,
                                               FFRecOptions options,
                                               FFRatFunList ** results);
  FFStatus ffReconstructFromCurrentEvaluationsMod(FFGraph graph,
                                                  FFRecOptions options,
                                                  FFRatFunList ** results);

  /*
   * Reconstruct by setting the degrees manually (as if using
   * ffSetDegrees) from `deg_data` which must contain, in row-major
   * order, the entries of a
   *
   *   nparsout x (2+4*nparsin)
   *
   * matrix with the `j`-th row specifying the degrees of the j-th
   * output element in the following form:
   *
   *   [numdeg, dendeg, nM1, nm1, dM1, dm1, nM2, nm2, dM2, dm2, ...]
   *
   * where `numdeg` and `dendeg` are the total degrees of numerator
   * and denominator, while the other entries must be read according
   * to:
   * - n = numerator's, d = denominator's
   * - M = maximum degree, m = minimum degree
   * - `i` (integer) = with respect to the `i`-th input parameter
   */
  FFStatus ffReconstructFunctionWithDegrees(FFGraph graph,
                                            FFRecOptions options,
                                            const unsigned * deg_data,
                                            FFRatFunList ** results);
  FFStatus ffReconstructFunctionWithDegreesMod(FFGraph graph,
                                               FFRecOptions options,
                                               const unsigned * deg_data,
                                               FFRatFunList ** results);

  /* API end */

#ifdef __cplusplus
}
#endif
