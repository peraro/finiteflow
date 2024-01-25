#include <algorithm>
#include <fflow/gcd.hh>
#include <fflow/mp_functions.hh>
#include <fflow/polynomial.hh>
#include <fflow/mp_multivariate_reconstruction.hh>
#include <fflow/mp_gcd.hh>
#include <fflow/primes.hh>
#include <fflow/ostream.h>
#include <fflow/json.hh>
#include <fflow/graph.hh>
#include <fflow/subgraph.hh>
#include <fflow/alg_lists.hh>
#include <fflow/alg_functions.hh>
#include <fflow/analytic_solver.hh>
#include <fflow/analytic_fit.hh>
#include <fflow/numeric_solver.hh>
#include <fflow/numeric_fit.hh>
#include <fflow/alg_laurent.hh>
#include <fflow/subgraph_fit.hh>
#include <fflow/subgraph_reconstruct.hh>
#include <fflow/node_solver.hh>
#include <fflow/cached_subgraph.hh>
#include <mathlink.h>
#include <WolframLibrary.h>
using namespace fflow;

namespace fflow {

  void fflowml_print(WolframLibraryData libData, const std::string & str)
  {
    MLINK link = libData->getWSLINK(libData);
    MLPutFunction(link, "EvaluatePacket", 1);
    MLPutFunction(link, "Print", 1);
    MLPutString(link, str.c_str());
    libData->processWSLINK(link);

    if (MLNextPacket(link) == RETURNPKT)
      MLNewPacket(link);
  }

  void fflowml_evaluate(WolframLibraryData libData, const std::string & str)
  {
    MLINK link = libData->getWSLINK(libData);
    MLPutFunction(link, "EvaluatePacket", 1);
    MLPutFunction(link, "ToExpression", 1);
    MLPutString(link, str.c_str());
    libData->processWSLINK(link);

    if (MLNextPacket(link) == RETURNPKT)
      MLNewPacket(link);
  }

} // namespace fflow


#ifndef FFLOW_NO_DBG

namespace  {

  class MathDbgPrint : public fflow::DBGPrint {
  public:

    explicit MathDbgPrint(WolframLibraryData libData)
      : libData_(libData) {}

    virtual void print(const std::string & msg)
    {
      fflowml_print(libData_, msg);
    }

  private:
    WolframLibraryData libData_;
  };

#define FFLOWML_SET_DBGPRINT() \
  MathDbgPrint math_dbg_print_(libData); \
  fflow::DBGPrintSet math_dbg_print_setter_(math_dbg_print_)

} // namespace

#else

#define FFLOWML_SET_DBGPRINT()

#endif


namespace  {

  // global data

  //const UInt DEFAULT_MOD = BIG_UINT_PRIMES[BIG_UINT_PRIMES_SIZE-1];

  // mathematica global context
  Session session;


  // internal classes, functions and methods

  struct MathGetHornerFun {

    template <typename FunMap>
    void get_horner_ratfun(MLINK mlp,
                           HornerRatFunPtr & f,
                           FunMap & map,
                           HornerWorkspacePtr & workspace);

    void get_horner_poly(MLINK mlp,
                         HornerPtr & f,
                         MPHornerMap & map,
                         HornerWorkspacePtr & workspace);

    void get_exponents_and_coeffs_(MLINK mlp, int nterms,
                                   Monomial * data, MPRational * coeffs,
                                   unsigned & nvars);

    void get_exponents_and_coeffs_(MLINK mlp, int nterms,
                                   Monomial * data, unsigned * coeffs,
                                   unsigned & nvars);

    unsigned nvars;
  };

  struct MathGetHornerCutFun {

    void get_horner_ratfun(MLINK mlp,
                           HornerRatFunPtr & f,
                           HornerHornerRatFunMap & map,
                           std::vector<MPHornerMap> & mpmap,
                           HornerWorkspacePtr & workspace);

    void get_exponents_and_coeffs_(MLINK mlp, int nterms,
                                   Monomial * data,
                                   HornerPtr * coeffs,
                                   MPHornerMap * mpmap,
                                   unsigned & nvars,
                                   HornerWorkspacePtr & ww);

    void get_horner_subpoly(MLINK mlp,
                            HornerPtr & poly, MPHornerMap & map,
                            HornerWorkspacePtr & ww);

    unsigned sub_nvars;
    unsigned nvars;
  };

  void get_rational(MLINK mlp, MPRational & c)
  {
    const char * coeff;
    MLGetString(mlp, &coeff);
    c = MPRational(coeff);
    MLReleaseString(mlp, coeff);
  }

  void get_integer(MLINK mlp, MPInt & c)
  {
    const char * coeff;
    MLGetString(mlp, &coeff);
    c = MPInt(coeff);
    MLReleaseString(mlp, coeff);
  }


  void MathGetHornerCutFun::get_horner_subpoly(MLINK mlp,
                                               HornerPtr & poly,
                                               MPHornerMap & map,
                                               HornerWorkspacePtr & ww)
  {
    MathGetHornerFun getpoly;
    getpoly.get_horner_poly(mlp, poly, map, ww);
    sub_nvars = getpoly.nvars;
  }


  void MathGetHornerFun::get_exponents_and_coeffs_(MLINK mlp, int nterms,
                                                   Monomial * monomials,
                                                   MPRational* coeffs,
                                                   unsigned & nvars)
  {
    long two;
    int * exponents;
    int n;
    unsigned this_deg;

    for (long i=0; i<nterms; ++i) {
      MLCheckFunction(mlp, "List", &two);

      MLGetInteger32List(mlp, &exponents, &n);
      monomials[i] = Monomial(n);
      std::copy(exponents, exponents+n, monomials[i].exponents());
      this_deg = std::accumulate(exponents, exponents+n, 0);
      monomials[i].degree() = this_deg;

      monomials[i].coeff() = 1;
      get_rational(mlp, coeffs[i]);

      MLReleaseInteger32List(mlp, exponents, n);
    }

    nvars = n;
  }


  template <typename FunMap>
  void MathGetHornerFun::get_horner_ratfun(MLINK mlp,
                                           HornerRatFunPtr & f,
                                           FunMap & map,
                                           HornerWorkspacePtr & ww)
  {
    long two;
    long nterms_num;
    long nterms_den;
    std::vector<Monomial> monomials_num;
    std::vector<Monomial> monomials_den;

    MLCheckFunction(mlp, "List", &two);

    // special case: a length-one list represents a vanishing function
    // containing the number of variables
    if (two == 1) {
      MLGetInteger32(mlp, (int*)&nvars);
      return;
    }

    // numerator

    MLCheckFunction(mlp, "List", &nterms_num);
    monomials_num.resize(nterms_num);
    map.num_map.resize(nterms_num);

    get_exponents_and_coeffs_(mlp, nterms_num,
                              monomials_num.data(),
                              map.num_map.coeff.get(), nvars);

    // denominator

    MLCheckFunction(mlp, "List", &nterms_den);
    monomials_den.resize(nterms_den);
    map.den_map.resize(nterms_den);

    get_exponents_and_coeffs_(mlp, nterms_den,
                              monomials_den.data(),
                              map.den_map.coeff.get(), nvars);

    // build result

    f.from_sparse_poly(monomials_num.data(), monomials_num.size(),
                       monomials_den.data(), monomials_den.size(),
                       nvars, 0,
                       map.num_map.pos.get(), map.den_map.pos.get());

    std::size_t max_workspace = horner_required_workspace(f.num());
    max_workspace = std::max(horner_required_workspace(f.den()),
                             max_workspace);
    ww.ensure_size(max_workspace);
  }


  void MathGetHornerFun::get_exponents_and_coeffs_(MLINK mlp, int nterms,
                                                   Monomial * monomials,
                                                   unsigned * coeffs,
                                                   unsigned & nvars)
  {
    long two;
    int * exponents;
    int n;
    unsigned this_deg;

    for (long i=0; i<nterms; ++i) {
      MLCheckFunction(mlp, "List", &two);

      MLGetInteger32List(mlp, &exponents, &n);
      monomials[i] = Monomial(n);
      std::copy(exponents, exponents+n, monomials[i].exponents());
      this_deg = std::accumulate(exponents, exponents+n, 0);
      monomials[i].degree() = this_deg;

      monomials[i].coeff() = 1;
      MLGetInteger32(mlp, (int*)(&coeffs[i]));

      MLReleaseInteger32List(mlp, exponents, n);
    }

    nvars = n;
  }


  void MathGetHornerFun::get_horner_poly(MLINK mlp,
                                         HornerPtr & f,
                                         MPHornerMap & map,
                                         HornerWorkspacePtr & ww)
  {
    long nterms;
    std::vector<Monomial> monomials;

    MLCheckFunction(mlp, "List", &nterms);
    monomials.resize(nterms);
    map.resize(nterms);

    get_exponents_and_coeffs_(mlp, nterms,
                              monomials.data(),
                              map.coeff.get(), nvars);

    f = hornerptr_from_sparse_poly(monomials.data(), nterms, nvars, 0,
                                   map.pos.get());
    ww.ensure_size(horner_required_workspace(f.get()));
  }


  void MathGetHornerCutFun::get_exponents_and_coeffs_(MLINK mlp, int nterms,
                                                      Monomial * monomials,
                                                      HornerPtr * coeffs,
                                                      MPHornerMap * mpmap,
                                                      unsigned & nvars,
                                                      HornerWorkspacePtr & ww)
  {
    long two;
    int * exponents;
    int n;
    unsigned this_deg;

    for (long i=0; i<nterms; ++i) {
      MLCheckFunction(mlp, "List", &two);

      MLGetInteger32List(mlp, &exponents, &n);
      monomials[i] = Monomial(n);
      std::copy(exponents, exponents+n, monomials[i].exponents());
      this_deg = std::accumulate(exponents, exponents+n, 0);
      monomials[i].degree() = this_deg;

      monomials[i].coeff() = 1;
      get_horner_subpoly(mlp, coeffs[i], mpmap[i], ww);

      MLReleaseInteger32List(mlp, exponents, n);
    }

    nvars = n;
  }


  void MathGetHornerCutFun::get_horner_ratfun(MLINK mlp,
                                              HornerRatFunPtr & f,
                                              HornerHornerRatFunMap & map,
                                              std::vector<MPHornerMap> & mpmap,
                                              HornerWorkspacePtr & ww)
  {
    long two;
    long nterms_num;
    long nterms_den;
    std::size_t mpmap_pos;
    std::vector<Monomial> monomials_num;
    std::vector<Monomial> monomials_den;

    MLCheckFunction(mlp, "List", &two);

    // numerator

    MLCheckFunction(mlp, "List", &nterms_num);
    monomials_num.resize(nterms_num);
    map.num_map.resize(nterms_num);
    mpmap_pos = mpmap.size(); // <-- This is zero
    mpmap.resize(mpmap_pos + nterms_num);

    get_exponents_and_coeffs_(mlp, nterms_num,
                              monomials_num.data(),
                              map.num_map.coeff.get(),
                              mpmap.data() + mpmap_pos,
                              nvars, ww);

    // denominator

    MLCheckFunction(mlp, "List", &nterms_den);
    monomials_den.resize(nterms_den);
    map.den_map.resize(nterms_den);
    mpmap_pos = mpmap.size();
    mpmap.resize(mpmap_pos+ nterms_den);

    get_exponents_and_coeffs_(mlp, nterms_den,
                              monomials_den.data(),
                              map.den_map.coeff.get(),
                              mpmap.data() + mpmap_pos,
                              nvars, ww);

    // build result

    f.from_sparse_poly(monomials_num.data(), monomials_num.size(),
                       monomials_den.data(), monomials_den.size(),
                       nvars, 0,
                       map.num_map.pos.get(), map.den_map.pos.get());

    std::size_t max_workspace = horner_required_workspace(f.num());
    max_workspace = std::max(horner_required_workspace(f.den()),
                             max_workspace);
    ww.ensure_size(max_workspace);
  }


  // to math

  void put_mprat(MLINK mlp, const MPRational & q,
                 void (*gmpfreefunc) (void *, size_t))
  {
    if (q.sign() == 0) {
      MLPutInteger(mlp, 0);
    } else {
      char * coeff = mpq_get_str(0, 10, q.get());
      MLPutString(mlp, coeff);
      (*gmpfreefunc)(coeff, std::strlen(coeff)+1);
    }
  }


  void put_sparse_poly(MLINK mlp, const MPReconstructedPoly & p)
  {
    const std::size_t psize = p.size();
    const std::size_t nvars = p.nvars();

    void (*gmpfreefunc) (void *, size_t);
    mp_get_memory_functions (0, 0, &gmpfreefunc);

    if (psize == 0) {
      Monomial mono(nvars);
      MLPutFunction(mlp, "List", 1);
      MLPutFunction(mlp, "List", 2);
      MLPutInteger16List(mlp, mono.exponents(), nvars);
      MLPutString(mlp, "0");
      return;
    }

    MLPutFunction(mlp, "List", psize);

    for (unsigned i=0; i<psize; ++i) {

      MLPutFunction(mlp, "List", 2);

      MLPutInteger16List(mlp, p.monomial(i).exponents(), nvars);

      char * coeff = mpq_get_str(0, 10, p.coeff(i).get());
      MLPutString(mlp, coeff);
      (*gmpfreefunc)(coeff, std::strlen(coeff)+1);
    }
  }

  void put_sparse_ratfun(MLINK mlp, const MPReconstructedRatFun & f)
  {
    MLPutFunction(mlp, "List", 2);
    put_sparse_poly(mlp, f.numerator());
    put_sparse_poly(mlp, f.denominator());
  }


  void put_alg_degree(MLINK mlp, unsigned numdeg, unsigned dendeg)
  {
    MLPutFunction(mlp, "FiniteFlow`FFFunDeg", 2);
    MLPutInteger32(mlp, numdeg);
    MLPutInteger32(mlp, dendeg);
  }


  template <typename T>
  inline void set_bolean_flag(T & opt, unsigned flag, int val)
  {
    if (val) {
      if (val == 1)
        opt |= flag;
      else
        opt &= ~flag;
    }
  }

  template <typename T, typename U>
  inline void set_uint_opt(T & opt, std::int64_t val, U defaultval)
  {
    if (val < 0)
      opt = defaultval;
    else
      opt = val;
  }


  void get_dense_system_data(MLINK mlp,
                             DynamicMatrixT<HornerRatFunPtr> & c,
                             DynamicMatrixT<MPHornerRatFunMap> & cmap,
                             HornerWorkspacePtr & ww,
                             unsigned & nparsin)
  {
    long n_rows, n_cols;
    int n_vars;
    MLGetInteger32(mlp, (int*)(&nparsin));
    MLGetInteger32(mlp, &n_vars);
    MLCheckFunction(mlp, "List", &n_rows);
    c.resize(n_rows, n_vars+1);
    cmap.resize(n_rows, n_vars+1);
    for (int i=0; i<n_rows; ++i) {
      MLCheckFunction(mlp, "List", &n_cols);
      for (int j=0; j<n_cols; ++j) {
        MathGetHornerFun getfun;
        getfun.get_horner_ratfun(mlp, c(i,j), cmap(i,j), ww);
        nparsin = getfun.nvars;
      }
    }
  }


  void get_sparse_system_data(MLINK mlp,
                              AnalyticSparseSolver & sys,
                              AnalyticSparseSolverData & data,
                              HornerWorkspacePtr & ww,
                              unsigned & nparsin)
  {
    long n_rows, two, lcsize;
    int csize;

    MLCheckFunction(mlp, "List", &n_rows);
    sys.rinfo.resize(n_rows);
    data.c.resize(n_rows);
    sys.cmap.resize(n_rows);
    for (int i=0; i<n_rows; ++i) {
      MLCheckFunction(mlp, "List", &two);
      {
        int * crinfo;
        MLGetInteger32List(mlp, &crinfo, &csize);
        AnalyticSparseSolver::RowInfo & rinf = sys.rinfo[i];
        rinf.size = csize;
        rinf.cols.reset(new unsigned[csize]);
        data.c[i].reset(new HornerRatFunPtr[csize]);
        sys.cmap[i].reset(new MPHornerRatFunMap[csize]);
        std::copy(crinfo, crinfo+csize, rinf.cols.get());
        MLReleaseInteger32List(mlp, crinfo, csize);
      }
      MLCheckFunction(mlp, "List", &lcsize);
      for (int j=0; j<csize; ++j) {
        MathGetHornerFun getfun;
        getfun.get_horner_ratfun(mlp, data.c[i][j], sys.cmap[i][j], ww);
        nparsin = getfun.nvars;
      }
    }
  }


  void get_num_dense_system_data(MLINK mlp,
                                 DynamicMatrixT<MPRational> & c)
  {
    long n_rows, n_cols;
    int n_vars;

    MLGetInteger32(mlp, &n_vars);
    MLCheckFunction(mlp, "List", &n_rows);
    c.resize(n_rows, n_vars+1);
    for (int i=0; i<n_rows; ++i) {
      MLCheckFunction(mlp, "List", &n_cols);
      for (int j=0; j<n_cols; ++j)
        get_rational(mlp, c(i,j));
    }
  }


  void get_num_sparse_system_data(MLINK mlp,
                                  NumericSparseSolver & sys)
  {
    long n_rows, two, lcsize;
    int csize;

    MLCheckFunction(mlp, "List", &n_rows);
    sys.rinfo.resize(n_rows);
    sys.c.resize(n_rows);
    for (int i=0; i<n_rows; ++i) {
      MLCheckFunction(mlp, "List", &two);
      {
        int * crinfo;
        MLGetInteger32List(mlp, &crinfo, &csize);
        NumericSparseSolver::RowInfo & rinf = sys.rinfo[i];
        rinf.size = csize;
        rinf.cols.reset(new unsigned[csize]);
        sys.c[i].reset(new MPRational[csize]);
        std::copy(crinfo, crinfo+csize, rinf.cols.get());
        MLReleaseInteger32List(mlp, crinfo, csize);
      }
      MLCheckFunction(mlp, "List", &lcsize);
      for (int j=0; j<csize; ++j)
        get_rational(mlp, sys.c[i][j]);
    }
  }


  void get_cut_solve_horner_data(MLINK mlp,
                                 std::vector<HornerRatFunPtr> & sfun,
                                 std::vector<HornerHornerRatFunMap> & map,
                                 std::vector<std::vector<MPHornerMap>> & mpmap,
                                 unsigned & nvars, unsigned & sub_nvars,
                                 HornerWorkspacePtr & ww)
  {
    sfun.clear();

    long n_cols;

    MLCheckFunction(mlp, "List", &n_cols);
    sfun.resize(n_cols);
    map.resize(n_cols);
    mpmap.resize(n_cols);

    for (int j=0; j<n_cols; ++j) {
        MathGetHornerCutFun getfun;
        getfun.get_horner_ratfun(mlp, sfun[j], map[j], mpmap[j], ww);
        nvars = getfun.nvars;
        sub_nvars = getfun.sub_nvars;
      }
  }

  void get_cut_solve_numeric_data(MLINK mlp,
                                  std::vector<HornerRatFunPtr> & sfun,
                                  std::vector<MPHornerRatFunMap> & map,
                                  unsigned & nvars,
                                  HornerWorkspacePtr & ww)
  {
    sfun.clear();

    long n_cols;

    MLCheckFunction(mlp, "List", &n_cols);
    sfun.resize(n_cols);
    map.resize(n_cols);

    for (int j=0; j<n_cols; ++j) {
      MathGetHornerFun getfun;
      getfun.get_horner_ratfun(mlp, sfun[j], map[j], ww);
      nvars = getfun.nvars;
    }
  }


  void put_alg_info(MLINK mlp, unsigned graphid, unsigned id)
  {
    if (!session.node_exists(graphid,id)) {
      MLPutSymbol(mlp, "$Failed");
    }

    Algorithm * alg = session.node(graphid,id)->algorithm();
    if (!alg) {
      MLPutSymbol(mlp, "Null");
    }

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
        DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);

        MLPutFunction(mlp, "List", 3);

        if (sizeof(std::size_t) == 8)
          MLPutInteger64List(mlp,
                             (mlint64*)ls.needed_depvars(),
                             ls.n_needed_depvars());
        else
          MLPutInteger32List(mlp,
                             (int*)ls.needed_depvars(),
                             ls.n_needed_depvars());

        if (sizeof(std::size_t) == 8)
          MLPutInteger64List(mlp,
                             (mlint64*)ls.needed_indepvars(),
                             ls.n_needed_indepvars());
        else
          MLPutInteger32List(mlp,
                             (int*)ls.needed_indepvars(),
                             ls.n_needed_indepvars());

        std::vector<int> zeroes;
        zeroes.reserve(20);
        for (unsigned i=0; i<ls.nvars(); ++i)
          if (!(ls.xinfo()[i] & LSVar::IS_NON_ZERO))
            zeroes.push_back(i);
        MLPutInteger32List(mlp, zeroes.data(), zeroes.size());

    } else if(dynamic_cast<SparseLinearSolver*>(alg)) {
        SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);

        MLPutFunction(mlp, "List", 3);

        if (sizeof(std::size_t) == 8)
          MLPutInteger64List(mlp,
                             (mlint64*)ls.needed_depvars(),
                             ls.n_needed_depvars());
        else
          MLPutInteger32List(mlp,
                             (int*)ls.needed_depvars(),
                             ls.n_needed_depvars());

        if (!ls.output_is_sparse()) {
          if (sizeof(std::size_t) == 8)
            MLPutInteger64List(mlp,
                               (mlint64*)ls.needed_indepvars(),
                               ls.n_needed_indepvars());
          else
            MLPutInteger32List(mlp,
                               (int*)ls.needed_indepvars(),
                               ls.n_needed_indepvars());
        } else {
          const auto & spdata = *(ls.sparse_output_data());
          MLPutFunction(mlp, "List", spdata.size());
          for (const auto & elist : spdata) {
            MLPutInteger32List(mlp,
                               (int*)elist.data(),
                               elist.size());
          }
        }

        MLPutInteger32(mlp, int(ls.output_is_sparse()));

    } else if (dynamic_cast<LinearFit*>(alg) ||
               dynamic_cast<SubgraphFit*>(alg)) {
        const LinearFit & ls = (dynamic_cast<LinearFit*>(alg) ?
                                *static_cast<LinearFit*>(alg)
                                : static_cast<SubgraphFit*>(alg)->linear_fit());

        MLPutFunction(mlp, "List", 3);

        if (sizeof(std::size_t) == 8)
          MLPutInteger64List(mlp,
                             (mlint64*)ls.needed_depvars(),
                             ls.n_needed_depvars());
        else
          MLPutInteger32List(mlp,
                             (int*)ls.needed_depvars(),
                             ls.n_needed_depvars());

        if (sizeof(std::size_t) == 8)
          MLPutInteger64List(mlp,
                             (mlint64*)ls.needed_indepvars(),
                             ls.n_needed_indepvars());
        else
          MLPutInteger32List(mlp,
                             (int*)ls.needed_indepvars(),
                             ls.n_needed_indepvars());

        std::vector<int> zeroes;
        zeroes.reserve(20);
        for (unsigned i=0; i<ls.ncoeffs(); ++i)
          if (!(ls.xinfo()[i] & LSVar::IS_NON_ZERO))
            zeroes.push_back(i);
        MLPutInteger32List(mlp, zeroes.data(), zeroes.size());

    } else if (dynamic_cast<SubgraphMultiFit*>(alg)) {
        const auto & multifit = *static_cast<SubgraphMultiFit*>(alg);
        std::size_t nfits = multifit.n_fits();
        MLPutFunction(mlp, "List", nfits);
        for (unsigned j=0; j<nfits; ++j) {
          const LinearFit & ls = multifit.linear_fit(j);

          MLPutFunction(mlp, "List", 3);

          if (sizeof(std::size_t) == 8)
            MLPutInteger64List(mlp,
                               (mlint64*)ls.needed_depvars(),
                               ls.n_needed_depvars());
          else
            MLPutInteger32List(mlp,
                               (int*)ls.needed_depvars(),
                               ls.n_needed_depvars());

          if (sizeof(std::size_t) == 8)
            MLPutInteger64List(mlp,
                               (mlint64*)ls.needed_indepvars(),
                               ls.n_needed_indepvars());
          else
            MLPutInteger32List(mlp,
                               (int*)ls.needed_indepvars(),
                               ls.n_needed_indepvars());

          std::vector<int> zeroes;
          zeroes.reserve(20);
          for (unsigned i=0; i<ls.ncoeffs(); ++i)
            if (!(ls.xinfo()[i] & LSVar::IS_NON_ZERO))
              zeroes.push_back(i);
          MLPutInteger32List(mlp, zeroes.data(), zeroes.size());
        }

    } else if (dynamic_cast<NonZeroes*>(alg)) {
        NonZeroes & salg = *static_cast<NonZeroes*>(alg);
        MLPutFunction(mlp, "List", 2);
        MLPutInteger32(mlp, salg.nparsin[0]);
        MLPutInteger32List(mlp, (int*)(salg.non_zeroes()), salg.nparsout);

    } else if (dynamic_cast<LaurentExpansion*>(alg)) {
        auto & salg = *static_cast<LaurentExpansion*>(alg);
        MLPutFunction(mlp, "List", 2);
        std::vector<int> pref_exp(salg.subgraph()->nparsout);
        salg.prefactor_exponent(pref_exp.data());
        MLPutInteger32List(mlp, pref_exp.data(), pref_exp.size());
        MLPutInteger32List(mlp, salg.order(), pref_exp.size());

    } else if(dynamic_cast<SubgraphRec*>(alg)) {
      SubgraphRec & ls = *static_cast<SubgraphRec*>(alg);
      const SparseRationalFunction * fun = ls.rec_function();
      const unsigned nout = ls.subgraph()->nparsout;
      const unsigned nrec = ls.n_rec_vars();
      MLPutFunction(mlp, "List", nout);
      for (unsigned j=0; j<nout; ++j) {
        MLPutFunction(mlp, "List", 2);
        MLPutFunction(mlp, "List", fun[j].numerator().size());
        for (const auto & mon : fun[j].numerator())
          MLPutInteger16List(mlp, mon.exponents(), nrec);
        MLPutFunction(mlp, "List", fun[j].denominator().size());
        for (const auto & mon : fun[j].denominator())
          MLPutInteger16List(mlp, mon.exponents(), nrec);
      }

    } else {
      MLPutSymbol(mlp, "Null");
    }
  }


  void get_input_nodes(MLINK mlp, std::vector<unsigned> & nodes)
  {
    int * data;
    int count=0;
    MLGetInteger32List(mlp, &data, &count);
    nodes.resize(count);
    std::copy(data, data+count, nodes.data());
    MLReleaseInteger32List(mlp, data, count);
  }


  void set_first_nparsin(Algorithm & alg, unsigned nparsin)
  {
    if (!alg.nparsin.size())
      alg.nparsin.resize(1);
    alg.nparsin[0] = nparsin;
  }


  ReconstructionOptions get_rec_opt(MLINK mlp, ReconstructionOptions opt)
  {
    int * opt_data;
    int opt_size;

    MLGetInteger32List(mlp, &opt_data, &opt_size);

    if (opt_data[0] >= 0)
      opt.n_checks = opt_data[0];
    if (opt_data[1] >= 0)
      opt.n_uchecks = opt_data[1];
    if (opt_data[2] >= 0)
      opt.n_singular = opt_data[2];
    if (opt_data[3] >= 0)
      opt.start_mod = opt_data[3];
    if (opt_data[4] >= 0)
      opt.max_primes = opt_data[4];
    if (opt_data[5] >= 0)
      opt.max_deg = opt_data[5];
    if (opt_data[6] >= 0)
      opt.dbginfo = opt_data[6];
    if (opt_data[7] >= 0)
      opt.polymethod = opt_data[7];

    MLReleaseInteger32List(mlp, opt_data, opt_size);

    return opt;
  }

  ReconstructionOptions get_rec_opt(MLINK mlp)
  {
    return get_rec_opt(mlp, ReconstructionOptions());
  }


} // namespace


extern "C" {

  mint WolframLibrary_getVersion()
  {
	return WolframLibraryVersion;
  }


  int WolframLibrary_initialize(WolframLibraryData libData)
  {
    (void)(libData);
	return 0;
  }


  int fflowml_mul_inv(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    mlint64 a, mod;
    int two;

    MLTestHead( mlp, "List", &two);
    MLGetInteger64(mlp, &a);
    MLGetInteger64(mlp, &mod);
    MLNewPacket(mlp);

    Mod modn(mod);
    MLPutInteger64(mlp, mul_inv(red_mod(a, modn), modn));

    return LIBRARY_NO_ERROR;
  }


  int fflowml_prime_no(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int prime_no;
    int two;

    MLTestHead( mlp, "List", &two);
    MLGetInteger32(mlp, &prime_no);
    MLNewPacket(mlp);

    if (prime_no < 0)
      prime_no = 0;

    if (unsigned(prime_no) < BIG_UINT_PRIMES_SIZE) {
      MLPutInteger64(mlp, BIG_UINT_PRIMES[prime_no]);
    } else {
      MLPutSymbol(mlp, "$Failed");
    }

    return LIBRARY_NO_ERROR;
  }


#if 0
  int fflowml_tree_level(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();
    int n_args, n_legs, n_flavs, n_vars;
    int * hels;
    mlint64 * flavs;

    // input
    MLTestHead(mlp, "List", &n_args);

    // number of legs
    MLGetInteger32(mlp, &n_legs);

    // external variables
    MLGetInteger32(mlp, &n_vars);

    // external particle types
    MLGetInteger64List(mlp, &flavs, &n_flavs);

    // define amplitude object
    MathTreeLevelAmp amp(n_legs, (UInt*)(flavs), n_vars);

    // clean up flavours
    MLReleaseInteger64List(mlp, flavs, n_flavs);

    // helicities
    MLGetInteger32List(mlp, &hels, &n_flavs);
    amp.setHelicities(hels);
    MLReleaseInteger32List(mlp, hels, n_flavs);

    // expressions for spinor variables
    get_external_spinors(mlp,
                         amp.sp, amp.spmap,
                         amp.pref, amp.prefmap, amp.ww);

    // reconstruction options
    std::vector<UInt> shift;
    UInt shift0 = 299159;
    MPRatFunReconstruction & rec = amp.getRatFunReconstructon();

    // get shift
    get_reconstruction_shift_flags(mlp, shift, shift0);
    adjust_shift(shift, n_vars);
    rec.setShift(shift.data());
    rec.setT0(shift0);

    // get uint32flags
    get_reconstriction_uint32_flags(mlp, rec);

    // finished reading input
    MLNewPacket(mlp);

    // reconstruct
    MPReconstructedRatFun fun;
    Ret ret = amp.get_tree(fun);

    if (ret != SUCCESS) {

      MLPutSymbol(mlp, "$Failed");

    } else {

      put_sparse_ratfun(mlp, fun);

    }

    return LIBRARY_NO_ERROR;
  }
#endif // 0


  // algorithms interface

  int fflowml_default_nthreads(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int zero;
    MLTestHead( mlp, "List", &zero);
    MLNewPacket(mlp);
    MLPutInteger32(mlp, Session::default_nthreads());

    return LIBRARY_NO_ERROR;
  }

  int fflowml_graph_new(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLTestHead( mlp, "List", &one);
    MLNewPacket(mlp);

    unsigned id = session.new_graph();
    MLPutInteger32(mlp, id);

    return LIBRARY_NO_ERROR;
  }

  int fflowml_graph_dummy(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLTestHead( mlp, "List", &one);
    int nin, nout;
    MLGetInteger32(mlp, &nin);
    MLGetInteger32(mlp, &nout);
    MLNewPacket(mlp);

    unsigned id = session.new_dummy_graph(nin, nout);
    MLPutInteger32(mlp, id);

    return LIBRARY_NO_ERROR;
  }

  int fflowml_graph_delete(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLTestHead( mlp, "List", &one);
    int id;
    MLGetInteger32(mlp, &id);
    MLNewPacket(mlp);

    session.delete_graph(id);

    MLPutSymbol(mlp, "Null");

    return LIBRARY_NO_ERROR;
  }

  int fflowml_node_delete(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLTestHead( mlp, "List", &one);
    int gid, nid;
    MLGetInteger32(mlp, &gid);
    MLGetInteger32(mlp, &nid);
    MLNewPacket(mlp);

    if (session.delete_node(gid,nid) == ALG_NO_ID)
      MLPutSymbol(mlp, "$Failed");
    else
      MLPutSymbol(mlp, "Null");

    return LIBRARY_NO_ERROR;
  }

  int fflowml_graph_set_out_node(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid, nodeid;
    MLGetInteger32(mlp, &graphid);
    MLGetInteger32(mlp, &nodeid);

    MLNewPacket(mlp);

    if (session.set_output_node(graphid, nodeid)==ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutSymbol(mlp, "Null");
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_graph_input_vars(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid, nvars;
    MLGetInteger32(mlp, &graphid);
    MLGetInteger32(mlp, &nvars);

    Graph * graph = session.graph(graphid);
    MLNewPacket(mlp);

    if (!graph || graph->set_input_vars(nvars) == ALG_NO_ID)
      MLPutSymbol(mlp, "$Failed");
    else
      MLPutSymbol(mlp, "Null");

    return LIBRARY_NO_ERROR;
  }

  int fflowml_graph_nparsout(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid;
    MLGetInteger32(mlp, &graphid);

    MLNewPacket(mlp);

    Graph * g = session.graph(graphid);

    if (g == nullptr || !g->has_learned()) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, g->nparsout);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_node_nparsout(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid, nodeid;
    MLGetInteger32(mlp, &graphid);
    MLGetInteger32(mlp, &nodeid);

    MLNewPacket(mlp);

    Node * node = session.node(graphid, nodeid);

    if (node == nullptr || !node->algorithm()->has_learned()) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, node->algorithm()->nparsout);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_node_set_mutable(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid, nodeid;
    MLGetInteger32(mlp, &graphid);
    MLGetInteger32(mlp, &nodeid);

    MLNewPacket(mlp);

    unsigned ret = session.set_node_mutable(graphid, nodeid);

    if (ret == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutSymbol(mlp, "Null");
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_graph_prune(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid;
    MLGetInteger32(mlp, &graphid);

    MLNewPacket(mlp);

    unsigned ret = session.prune_graph(graphid);

    if (ret == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutSymbol(mlp, "Null");
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_graph_edges(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid, markedonly;
    MLGetInteger32(mlp, &graphid);
    MLGetInteger32(mlp, &markedonly);
    markedonly = (markedonly == 1);

    MLNewPacket(mlp);

    Graph * graph = session.graph(graphid);
    if (!graph) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    std::vector<unsigned> edges;
    if (markedonly)
      graph->marked_edges(edges);
    else
      graph->edges(edges);

    MLPutInteger32List(mlp, (int*)(edges.data()), edges.size());

    return LIBRARY_NO_ERROR;
  }

  int fflowml_graph_nodes(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid, markedonly;
    MLGetInteger32(mlp, &graphid);
    MLGetInteger32(mlp, &markedonly);
    markedonly = (markedonly == 1);

    MLNewPacket(mlp);

    Graph * graph = session.graph(graphid);
    if (!graph) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    std::vector<unsigned> nodes;
    if (markedonly)
      graph->marked_nodes(nodes);
    else
      graph->nodes(nodes);

    MLPutInteger32List(mlp, (int*)(nodes.data()), nodes.size());

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_simple_subgraph(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);
    int graphid, subgraphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &subgraphid);

    typedef SubGraphData Data;
    std::unique_ptr<SimpleSubGraph> algptr(new SimpleSubGraph());
    std::unique_ptr<Data> data(new Data());

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    Ret ret = algptr->init(session, subgraphid, *data,
                           nparsin.data(), nparsin.size());

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_memoized_subgraph(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);
    int graphid, subgraphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &subgraphid);

    typedef MemoizedSubGraphData Data;
    std::unique_ptr<MemoizedSubGraph> algptr(new MemoizedSubGraph());
    std::unique_ptr<Data> data(new Data());

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    Ret ret = algptr->init(session, subgraphid, *data,
                           nparsin.data(), nparsin.size());

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_subgraph_map(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);
    int graphid, subgraphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &subgraphid);

    typedef SubGraphData Data;
    std::unique_ptr<SubGraphMap> algptr(new SubGraphMap());
    std::unique_ptr<Data> data(new Data());

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    Ret ret = algptr->init(session, subgraphid, *data,
                           nparsin.data(), nparsin.size());

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_dense_system(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);

    get_input_nodes(mlp, inputnodes);

    unsigned nparsin = 0;

    typedef AnalyticDenseSolverData Data;
    std::unique_ptr<AnalyticDenseSolver> algptr(new AnalyticDenseSolver());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;
    get_dense_system_data(mlp, data->c, alg.cmap,
                          session.main_context()->ww, nparsin);
    set_first_nparsin(alg,nparsin);

    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);
    alg.init(data->c.rows(), data->c.columns()-1,
             (unsigned*)neededv, needed_size, *data);
    MLReleaseInteger32List(mlp, neededv, needed_size);
    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_sparse_system(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs, nvars;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);

    get_input_nodes(mlp, inputnodes);

    unsigned nparsin = 0;

    MLGetInteger32(mlp, (int*)(&nparsin));
    MLGetInteger32(mlp, &nvars);

    typedef AnalyticSparseSolverData Data;
    std::unique_ptr<AnalyticSparseSolver> algptr(new AnalyticSparseSolver());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;
    get_sparse_system_data(mlp, alg,
                           *data, session.main_context()->ww, nparsin);
    set_first_nparsin(alg,nparsin);

    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);
    alg.init(data->c.size(), nvars, (unsigned*)neededv, needed_size, *data);
    MLReleaseInteger32List(mlp, neededv, needed_size);
    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_num_dense_system(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);

    get_input_nodes(mlp, inputnodes);

    typedef NumericDenseSolverData Data;
    std::unique_ptr<NumericDenseSolver> algptr(new NumericDenseSolver());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;
    get_num_dense_system_data(mlp, alg.c);

    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);
    alg.init(alg.c.rows(), alg.c.columns()-1, (unsigned*)neededv, needed_size,
             *data);
    MLReleaseInteger32List(mlp, neededv, needed_size);
    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_node_dense_system(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);

    get_input_nodes(mlp, inputnodes);

    typedef NodeDenseSolverData Data;
    std::unique_ptr<NodeDenseSolver> algptr(new NodeDenseSolver());
    std::unique_ptr<Data> data(new Data());
    int neqs = 0, nvars = 0;
    MLGetInteger32(mlp, &neqs);
    MLGetInteger32(mlp, &nvars);
    auto & alg = *algptr;

    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);
    alg.init_node_solver(neqs, nvars, (unsigned*)neededv, needed_size, *data);
    MLReleaseInteger32List(mlp, neededv, needed_size);
    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_node_sparse_system(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs, nvars;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);

    get_input_nodes(mlp, inputnodes);

    MLGetInteger32(mlp, &nvars);

    typedef NodeSparseSolverData Data;
    std::unique_ptr<NodeSparseSolver> algptr(new NodeSparseSolver());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;
    auto & sys = alg;

    long n_rows;
    int csize;
    MLCheckFunction(mlp, "List", &n_rows);
    sys.rinfo.resize(n_rows);
    unsigned ncoeffs = 0;
    for (int i=0; i<n_rows; ++i) {
      {
        int * crinfo;
        MLGetInteger32List(mlp, &crinfo, &csize);
        NodeSparseSolver::RowInfo & rinf = sys.rinfo[i];
        rinf.start = ncoeffs;
        rinf.size = csize;
        rinf.cols.reset(new unsigned[csize]);
        std::copy(crinfo, crinfo+csize, rinf.cols.get());
        ncoeffs += csize;
        MLReleaseInteger32List(mlp, crinfo, csize);
      }
    }

    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);
    alg.init_node_solver(sys.rinfo.size(), nvars,
                         (unsigned*)neededv, needed_size, *data);
    MLReleaseInteger32List(mlp, neededv, needed_size);
    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_num_sparse_system(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs, nvars;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);
    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);

    get_input_nodes(mlp, inputnodes);

    MLGetInteger32(mlp, &nvars);

    typedef NumericSparseSolverData Data;
    std::unique_ptr<NumericSparseSolver> algptr(new NumericSparseSolver());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;
    get_num_sparse_system_data(mlp, alg);

    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);
    alg.init(alg.c.size(), nvars, (unsigned*)neededv, needed_size, *data);
    MLReleaseInteger32List(mlp, neededv, needed_size);
    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_json_sparse_system(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);

    get_input_nodes(mlp, inputnodes);

    typedef AnalyticSparseSolverData Data;
    std::unique_ptr<AnalyticSparseSolver> algptr(new AnalyticSparseSolver());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;

    const char * filename;
    MLGetString(mlp, &filename);
    unsigned needed_workspace;
    Ret ret = json_sparse_system(filename, alg, *data, needed_workspace);
    MLReleaseString(mlp, filename);

    MLNewPacket(mlp);

    if (ret != SUCCESS || !session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
      session.main_context()->ww.ensure_size(needed_workspace);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_json_ratfun(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);

    get_input_nodes(mlp, inputnodes);

    typedef AnalyticFunctionData Data;
    std::unique_ptr<AnalyticFunction> algptr(new AnalyticFunction());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;

    const char * filename;
    MLGetString(mlp, &filename);
    unsigned needed_workspace;
    Ret ret = json_sparse_ratfun(filename, alg, *data, needed_workspace);
    MLReleaseString(mlp, filename);

    MLNewPacket(mlp);

    if (ret != SUCCESS || !session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
      session.main_context()->ww.ensure_size(needed_workspace);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_system_reset_neeed(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int id, nodeid, nargs;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);
    bool okay = true;

    Algorithm * alg = session.algorithm(id, nodeid);
    if (!alg || !alg->is_mutable())
      okay = false;

    if (okay) {
      if (dynamic_cast<DenseLinearSolver *>(alg)) {
          DenseLinearSolver & ls = *static_cast<DenseLinearSolver *>(alg);
          Ret ret = ls.reset_needed(session.alg_data(id, nodeid),
                                    (unsigned*)neededv, needed_size);
          if (ret == SUCCESS)
            session.invalidate_subctxt_alg_data(id, nodeid);
          else
            okay = false;

      } else if(dynamic_cast<SparseLinearSolver *>(alg)) {
          SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
          if (ls.marked_and_sweeped()) {
            okay = false;
          } else {
            Ret ret = ls.reset_needed(session.alg_data(id, nodeid),
                                      (unsigned*)neededv, needed_size);
            if (ret == SUCCESS)
              session.invalidate_subctxt_alg_data(id, nodeid);
            else
              okay = false;
          }

      } else {
        okay = false;
      }
    }

    MLReleaseInteger32List(mlp, neededv, needed_size);
    MLNewPacket(mlp);

    if (okay)
      MLPutSymbol(mlp, "Null");
    else
      MLPutSymbol(mlp, "$Failed");

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_system_only_homogeneous(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int id, nodeid, nargs;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    bool okay = true;

    Algorithm * alg = session.algorithm(id, nodeid);
    if (!alg || !alg->is_mutable())
      okay = false;

    if (okay) {
      if (dynamic_cast<DenseLinearSolver *>(alg)) {
          DenseLinearSolver & ls = *static_cast<DenseLinearSolver *>(alg);
          ls.only_homogeneous();
          session.invalidate_subctxt_alg_data(id, nodeid);

      } else if (dynamic_cast<SparseLinearSolver *>(alg)) {
          SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
          ls.only_homogeneous();
          session.invalidate_subctxt_alg_data(id, nodeid);

      } else {
        okay = false;
      }
    }
    MLNewPacket(mlp);

    if (okay)
      MLPutSymbol(mlp, "Null");
    else
      MLPutSymbol(mlp, "$Failed");

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_system_sparse_output(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int id, nodeid, nargs;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    bool okay = true;

    Algorithm * alg = session.algorithm(id, nodeid);
    if (!alg || !alg->is_mutable())
      okay = false;

    if (okay) {
      if (dynamic_cast<SparseLinearSolver *>(alg)) {
        SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
        ls.sparse_output();
        session.invalidate_subctxt_alg_data(id, nodeid);
      } else {
        okay = false;
      }
    }
    MLNewPacket(mlp);

    if (okay)
      MLPutSymbol(mlp, "Null");
    else
      MLPutSymbol(mlp, "$Failed");

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_system_sparse_output_with_maxrow(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int id, nodeid, nargs, maxrow, back_subst_flag;
    bool back_subst = true;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    MLGetInteger32(mlp, &maxrow);
    MLGetInteger32(mlp, &back_subst_flag);
    bool okay = true;

    if(back_subst_flag == -1)
      back_subst = false;

    Algorithm * alg = session.algorithm(id, nodeid);
    if (!alg || !alg->is_mutable())
      okay = false;

    if (okay) {
      if (dynamic_cast<SparseLinearSolver *>(alg)) {
        SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
        Ret ret = ls.sparse_output_with_maxrow(maxrow, back_subst);
        if (ret == SUCCESS)
          session.invalidate_subctxt_alg_data(id, nodeid);
        else
          okay = false;
      } else {
        okay = false;
      }
    }
    MLNewPacket(mlp);

    if (okay)
      MLPutSymbol(mlp, "Null");
    else
      MLPutSymbol(mlp, "$Failed");

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_mark_and_sweep_eqs(WolframLibraryData libData, MLINK mlp)
  {
    // note: The algorithm doesn't need to be mutable, because
    // input/output are unchanged.  AlgorithmData is however
    // invalidated in subcontexts

    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int id, nodeid, nargs;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    Algorithm * alg = session.algorithm(id, nodeid);
    MLNewPacket(mlp);

    if (!alg) {
      MLPutSymbol(mlp, "$Failed");
    } else if (dynamic_cast<SparseLinearSolver *>(alg)) {
      SparseLinearSolver & ls = *static_cast<SparseLinearSolver *>(alg);
      if (ls.marked_and_sweeped()) {
        MLPutSymbol(mlp, "$Failed");
      } else {
        ls.mark_and_sweep_eqs(session.alg_data(id, nodeid));
        session.invalidate_subctxt_alg_data(id, nodeid);
        MLPutInteger32(mlp, ls.n_indep_eqs());
      }
    } else {
      MLPutSymbol(mlp, "$Failed");
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_delete_unneeded_eqs(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int id, nodeid, nargs;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    MLNewPacket(mlp);

    Algorithm * alg = session.algorithm(id, nodeid);
    if (!alg) {
      MLPutSymbol(mlp, "$Failed");

    } else if (dynamic_cast<AnalyticSparseSolver *>(alg)) {
      AnalyticSparseSolver & ls = *static_cast<AnalyticSparseSolver *>(alg);
      ls.delete_unneeded_eqs(session.alg_data(id, nodeid));
      session.invalidate_subctxt_alg_data(id, nodeid);
      MLPutSymbol(mlp, "Null");

    } else if (dynamic_cast<NumericSparseSolver *>(alg)) {
      NumericSparseSolver & ls = *static_cast<NumericSparseSolver *>(alg);
      ls.delete_unneeded_eqs(session.alg_data(id, nodeid));
      session.invalidate_subctxt_alg_data(id, nodeid);
      MLPutSymbol(mlp, "Null");

    } else {
      MLPutSymbol(mlp, "$Failed");
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_learn(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &nargs);

    int id;
    MLGetInteger32(mlp, &id);
    MLNewPacket(mlp);

    Ret ret = session.learn(id);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    put_alg_info(mlp, id, session.graph(id)->out_node()->id());

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_get_info(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &nargs);

    int id, nodeid;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    MLNewPacket(mlp);

    if (!session.graph_exists(id)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    if (nodeid < 0)
      nodeid = session.graph(id)->out_node_id();

    put_alg_info(mlp, id, nodeid);

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_set_learning_options(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &nargs);

    int id, nodeid, primeno, maxsingular;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    MLGetInteger32(mlp, &primeno);
    MLGetInteger32(mlp, &maxsingular);
    MLNewPacket(mlp);

    if (!session.node_exists(id, nodeid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    LearningOptions & opt = session.node(id, nodeid)->learn_opt;
    if (primeno >= 0)
      opt.prime_no = primeno;
    if (maxsingular >= 0)
      opt.n_singular = maxsingular;

    MLPutSymbol(mlp, "Null");

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_n_indep_eqs(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &nargs);

    int id, nodeid;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    MLNewPacket(mlp);

    Algorithm * alg = session.algorithm(id, nodeid);

    if (!session.graph_exists(id)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
        DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
        MLPutInteger32(mlp, ls.n_indep_eqs());

    } else if (dynamic_cast<SparseLinearSolver*>(alg)) {
        SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
        MLPutInteger32(mlp, ls.n_indep_eqs());

    } else {
        MLPutSymbol(mlp, "$Failed");
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_indep_eqs(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &nargs);

    int id, nodeid;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    MLNewPacket(mlp);

    Algorithm * alg = session.algorithm(id, nodeid);

    if (!session.graph_exists(id)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    if (dynamic_cast<DenseLinearSolver*>(alg)) {
        DenseLinearSolver & ls = *static_cast<DenseLinearSolver*>(alg);
        std::vector<int> indep_eqs(ls.n_indep_eqs());
        std::copy(ls.indep_eqs(), ls.indep_eqs()+ls.n_indep_eqs(),
                  indep_eqs.data());
        MLPutInteger32List(mlp, indep_eqs.data(), indep_eqs.size());

    } else if (dynamic_cast<SparseLinearSolver*>(alg)) {
        SparseLinearSolver & ls = *static_cast<SparseLinearSolver*>(alg);
        std::vector<int> indep_eqs(ls.n_indep_eqs());
        std::copy(ls.indep_eqs(), ls.indep_eqs()+ls.n_indep_eqs(),
                  indep_eqs.data());
        MLPutInteger32List(mlp, indep_eqs.data(), indep_eqs.size());

    } else {
      MLPutSymbol(mlp, "$Failed");

    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_degrees(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &nargs);

    int id;
    MLGetInteger32(mlp, &id);
    ReconstructionOptions opt = get_rec_opt(mlp);
    MLNewPacket(mlp);

    Ret ret = session.degrees(id, opt);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    const Graph & g = *session.graph(id);

    const unsigned * numdeg = g.degs_data().numdeg.get();
    const unsigned * dendeg = g.degs_data().dendeg.get();
    unsigned nparsout = g.nparsout;

    MLPutFunction(mlp, "List", nparsout);
    for (unsigned i=0; i<nparsout; ++i)
      put_alg_degree(mlp, numdeg[i], dendeg[i]);

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_vars_degrees(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &one);

    int id;
    MLGetInteger32(mlp, &id);
    ReconstructionOptions opt = get_rec_opt(mlp);
    MLNewPacket(mlp);

    Ret ret = session.all_var_degrees(id, opt);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    MLPutSymbol(mlp, "Null");

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_all_degrees(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &nargs);

    int id, nthreads;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nthreads);
    if (nthreads < 0)
      nthreads = 0;
    ReconstructionOptions opt = get_rec_opt(mlp);
    MLNewPacket(mlp);

    Ret ret = session.parallel_all_degrees(id, nthreads, opt);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    const Graph & g = *session.graph(id);

    const unsigned * numdeg = g.degs_data().numdeg.get();
    const unsigned * dendeg = g.degs_data().dendeg.get();
    unsigned nparsout = g.nparsout;

    MLPutFunction(mlp, "List", nparsout);
    for (unsigned i=0; i<nparsout; ++i)
      put_alg_degree(mlp, numdeg[i], dendeg[i]);

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_sample(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int two;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &two);

    int id, nthreads;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nthreads);
    if (nthreads < 0)
      nthreads = 0;
    ReconstructionOptions opt = get_rec_opt(mlp);
    MLNewPacket(mlp);

    session.parallel_sample(id, nthreads, opt);

    MLPutSymbol(mlp, "Null");

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_sample_from_points(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int two;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &two);

    int id, nthreads, samples_start, samples_size;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nthreads);
    if (nthreads < 0)
      nthreads = 0;
    ReconstructionOptions opt = get_rec_opt(mlp);

    std::string file;
    {
      const char * filename;
      MLGetString(mlp, &filename);
      file = std::string(filename);
      MLReleaseString(mlp, filename);
    }
    MLGetInteger32(mlp, &samples_start);
    MLGetInteger32(mlp, &samples_size);

    MLNewPacket(mlp);

    SamplePointsFromFile pts(file.c_str(), samples_start, samples_size);
    session.parallel_sample(id, nthreads, opt, &pts);

    MLPutSymbol(mlp, "Null");

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_reconstruct(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &one);

    int id, nthreads;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nthreads);
    if (nthreads < 0)
      nthreads = 0;
    ReconstructionOptions opt = get_rec_opt(mlp);
    MLNewPacket(mlp);

    Graph * g = session.graph(id);

    if (!g) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    typedef MPReconstructedRatFun ResT;
    const unsigned nparsout = g->nparsout;
    std::unique_ptr<ResT[]> res (new ResT[nparsout]);

    Ret ret = session.parallel_reconstruct(id, res.get(), nthreads, opt);

    if (ret != SUCCESS) {

      if (ret == MISSING_SAMPLES) {
        MLPutSymbol(mlp, "FiniteFlow`FFMissingPoints");
        return LIBRARY_NO_ERROR;
      }

      if (ret == MISSING_PRIMES) {
        MLPutSymbol(mlp, "FiniteFlow`FFMissingPrimes");
        return LIBRARY_NO_ERROR;
      }

      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    MLPutFunction(mlp, "List", nparsout);
    for (unsigned i=0; i<nparsout; ++i)
      put_sparse_ratfun(mlp, res[i]);

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_reconstruct_mod(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &one);

    int id, nthreads;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nthreads);
    if (nthreads < 0)
      nthreads = 0;
    ReconstructionOptions opt = get_rec_opt(mlp);
    MLNewPacket(mlp);

    Graph * g = session.graph(id);

    if (!g) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    typedef MPReconstructedRatFun ResT;
    const unsigned nparsout = g->nparsout;
    std::unique_ptr<ResT[]> res (new ResT[nparsout]);

    Ret ret = session.parallel_reconstruct_mod(id, res.get(), nthreads, opt);

    if (ret != SUCCESS) {

      if (ret == MISSING_SAMPLES) {
        MLPutSymbol(mlp, "FiniteFlow`FFMissingPoints");
        return LIBRARY_NO_ERROR;
      }

      if (ret == MISSING_PRIMES) {
        MLPutSymbol(mlp, "FiniteFlow`FFMissingPrimes");
        return LIBRARY_NO_ERROR;
      }

      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    MLPutFunction(mlp, "List", nparsout);
    for (unsigned i=0; i<nparsout; ++i)
      put_sparse_ratfun(mlp, res[i]);

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_univariate_reconstruct(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &one);

    int id;
    MLGetInteger32(mlp, &id);
    ReconstructionOptions opt = get_rec_opt(mlp);
    MLNewPacket(mlp);

    Graph * g = session.graph(id);

    if (!g) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    typedef MPReconstructedRatFun ResT;
    const unsigned nparsout = g->nparsout;
    std::unique_ptr<ResT[]> res (new ResT[nparsout]);

    Ret ret = session.reconstruct_univariate(id, res.get(), opt);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    MLPutFunction(mlp, "List", nparsout);
    for (unsigned i=0; i<nparsout; ++i)
      put_sparse_ratfun(mlp, res[i]);

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_univariate_reconstruct_mod(WolframLibraryData libData,
                                             MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &one);

    int id;
    MLGetInteger32(mlp, &id);
    ReconstructionOptions opt = get_rec_opt(mlp);
    MLNewPacket(mlp);

    Graph * g = session.graph(id);

    if (!g) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    typedef MPReconstructedRatFun ResT;
    const unsigned nparsout = g->nparsout;
    std::unique_ptr<ResT[]> res (new ResT[nparsout]);

    Ret ret = session.reconstruct_univariate_mod(id, res.get(), opt);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    MLPutFunction(mlp, "List", nparsout);
    for (unsigned i=0; i<nparsout; ++i)
      put_sparse_ratfun(mlp, res[i]);

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_numeric_reconstruct(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int one;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &one);

    int id;
    MLGetInteger32(mlp, &id);
    ReconstructionOptions opt;
    opt.max_primes = 20;
    opt = get_rec_opt(mlp, opt);
    MLNewPacket(mlp);

    Graph * g = session.graph(id);

    if (!g) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    typedef MPRational ResT;
    const unsigned nparsout = g->nparsout;
    std::unique_ptr<ResT[]> res (new ResT[nparsout]);

    Ret ret = session.reconstruct_numeric(id, res.get(), opt);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    void (*gmpfreefunc) (void *, size_t);
    mp_get_memory_functions(0, 0, &gmpfreefunc);

    MLPutFunction(mlp, "List", nparsout);
    for (unsigned i=0; i<nparsout; ++i)
      put_mprat(mlp, res[i], gmpfreefunc);

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_independent_of_var(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int id, var, n_args;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &var);
    MLNewPacket(mlp);

    UInt ret = session.independent_of_var(id, var);

    if (ret == 0)
      MLPutSymbol(mlp, "False");
    else if (ret == 1)
      MLPutSymbol(mlp, "True");
    else if (ret == FAILED)
      MLPutSymbol(mlp, "$Failed");

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_nonzero(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid) || inputnodes.size() != 1) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);

    std::unique_ptr<NonZeroes> algptr(new NonZeroes());
    auto & alg = *algptr;

    Node * innode = session.node(graphid, inputnodes[0]);
    if (!innode) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    alg.init(innode->algorithm()->nparsout);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_ratfun_eval(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args, nfunctions;
    int nparsin=0;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    typedef AnalyticFunctionData Data;
    std::unique_ptr<AnalyticFunction> algptr(new AnalyticFunction());
    std::unique_ptr<Data> data(new Data());
    AnalyticFunction & alg = *algptr;
    MLGetInteger32(mlp, &nparsin);
    MLTestHead( mlp, "List", &nfunctions);
    data->f.reset(new HornerRatFunPtr[nfunctions]);
    alg.fmap.reset(new MPHornerRatFunMap[nfunctions]);
    for (int i=0; i<nfunctions; ++i) {
      MathGetHornerFun getfun;
      getfun.get_horner_ratfun(mlp, data->f[i], alg.fmap[i],
                               session.main_context()->ww);
    }

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid) || inputnodes.size() != 1) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);

    alg.init(nparsin, nfunctions, *data);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_coeff_ratfun_eval(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args, nfunctions, ncoeffs;
    int nparsin=0;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    typedef FunctionFromCoeffsData Data;
    std::unique_ptr<FunctionFromCoeffs> algptr(new FunctionFromCoeffs());
    std::unique_ptr<Data> data(new Data());
    FunctionFromCoeffs & alg = *algptr;
    MLGetInteger32(mlp, &nparsin);
    MLGetInteger32(mlp, &ncoeffs);
    MLTestHead( mlp, "List", &nfunctions);
    data->f.reset(new HornerRatFunPtr[nfunctions]);
    alg.fmap.reset(new CoeffHornerRatFunMap[nfunctions]);
    for (int i=0; i<nfunctions; ++i) {
      MathGetHornerFun getfun;
      getfun.get_horner_ratfun(mlp, data->f[i], alg.fmap[i],
                               session.main_context()->ww);
    }

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid) || inputnodes.size() != 2) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);

    alg.init(ncoeffs, nparsin, nfunctions, *data);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_ratnum_eval(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args, nnumbers;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    typedef EvalRationalNumbersData Data;
    std::unique_ptr<EvalRationalNumbers> algptr(new EvalRationalNumbers());
    std::unique_ptr<Data> data(new Data());
    EvalRationalNumbers & alg = *algptr;
    MLTestHead( mlp, "List", &nnumbers);
    std::vector<MPRational> numbers(nnumbers);
    for (int i=0; i<nnumbers; ++i)
      get_rational(mlp, numbers[i]);

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid) || inputnodes.size() != 0) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    alg.init(std::move(numbers), *data);

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_reset_ratnum_eval(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args, nnumbers;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid, nodeid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    MLGetInteger32(mlp, &nodeid);

    bool okay = true;

    Graph * graph  = session.graph(graphid);
    Algorithm * alg = session.algorithm(graphid, nodeid);
    if (!alg || !graph->is_mutable())
      okay = false;

    if (okay && dynamic_cast<EvalRationalNumbers *>(alg)) {

      auto & aa = *static_cast<EvalRationalNumbers *>(alg);
      MPRational * numbers = aa.numbers();

      MLTestHead( mlp, "List", &nnumbers);
      if (nnumbers != aa.nparsout) {
        okay = false;
      } else {
        for (int i=0; i<nnumbers; ++i)
          get_rational(mlp, numbers[i]);
        graph->invalidate_reconstruction_cache();
        static_cast<EvalRationalNumbersData*>(session.alg_data(graphid, nodeid))->invalidate();
        session.invalidate_subctxt_alg_data(graphid, nodeid);
      }
    } else {
      okay = false;
    }

    MLNewPacket(mlp);

    if (!okay)
      MLPutSymbol(mlp, "$Failed");
    else
      MLPutSymbol(mlp, "Null");

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_chain(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    std::unique_ptr<Chain> algptr(new Chain());
    Chain & alg = *algptr;
    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }
    alg.init(nparsin.data(), nparsin.size());

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_take(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    std::unique_ptr<Take> algptr(new Take());
    Take & alg = *algptr;
    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    std::vector<Take::InputEl> elems;
    int n_elems;
    MLTestHead(mlp, "List", &n_elems);
    elems.resize(n_elems/2);
    for (int j=0; j<n_elems/2; ++j) {
      int tmp;
      MLGetInteger32(mlp, &tmp);
      elems[j].list = tmp;
      MLGetInteger32(mlp, &tmp);
      elems[j].el = tmp;
    }

    Ret ret = alg.init(nparsin.data(), nparsin.size(), std::move(elems));

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid) || ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_slice(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args, start, end;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &start);
    MLGetInteger32(mlp, &end);

    MLNewPacket(mlp);

    Node * node = nullptr;
    if (!inputnodes.size() || !(node = session.node(graphid, inputnodes[0]))) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    unsigned nparsin = node->algorithm()->nparsout;
    std::unique_ptr<Slice> algptr(new Slice());
    Slice & alg = *algptr;
    if (end == -1)
      end = nparsin;
    Ret ret = alg.init(nparsin, start, end);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_add(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    MLNewPacket(mlp);

    std::unique_ptr<Add> algptr(new Add());
    Add & alg = *algptr;
    unsigned list_len = 0;
    if (inputnodes.size()) {
      const Node * n = session.node(graphid, inputnodes[0]);
      if (n)
        list_len = n->algorithm()->nparsout;
    }
    alg.init(inputnodes.size(), list_len);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_mul(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    MLNewPacket(mlp);

    std::unique_ptr<Mul> algptr(new Mul());
    Mul & alg = *algptr;
    unsigned list_len = 0;
    if (inputnodes.size()) {
      const Node * n = session.node(graphid, inputnodes[0]);
      if (n)
        list_len = n->algorithm()->nparsout;
    }
    alg.init(inputnodes.size(), list_len);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_mat_mul(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args, nrows1, ncols1, ncols2;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    std::unique_ptr<MatrixMul> algptr(new MatrixMul());
    auto & alg = *algptr;

    MLGetInteger32(mlp, &nrows1);
    MLGetInteger32(mlp, &ncols1);
    MLGetInteger32(mlp, &ncols2);

    MLNewPacket(mlp);

    alg.init(nrows1, ncols1, ncols2);

    if (!session.graph_exists(graphid) || inputnodes.size() != 2) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_take_and_add(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    std::unique_ptr<TakeAndAdd> algptr(new TakeAndAdd());
    TakeAndAdd & alg = *algptr;
    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    std::vector<std::vector<TakeAndAdd::InputEl>> elems;
    int n_lists, n_elems=0;
    MLTestHead(mlp, "List", &n_lists);
    elems.resize(n_lists);
    for (int k=0; k<n_lists; ++k) {
      MLTestHead(mlp, "List", &n_elems);
      elems[k].resize(n_elems/2);
      for (int j=0; j<n_elems/2; ++j) {
        int tmp;
        MLGetInteger32(mlp, &tmp);
        elems[k][j].list = tmp;
        MLGetInteger32(mlp, &tmp);
        elems[k][j].el = tmp;
      }
    }

    Ret ret = alg.init(nparsin.data(), nparsin.size(), std::move(elems));

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid) || ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


    int fflowml_alg_take_and_add_bl(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    std::unique_ptr<TakeAndAddBL> algptr(new TakeAndAddBL());
    TakeAndAddBL & alg = *algptr;
    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    std::vector<std::vector<TakeAndAddBL::InputEl>> elems;
    int n_lists, n_elems=0;
    MLTestHead(mlp, "List", &n_lists);
    elems.resize(n_lists);
    for (int k=0; k<n_lists; ++k) {
      MLTestHead(mlp, "List", &n_elems);
      elems[k].resize(n_elems/4);
      for (int j=0; j<n_elems/4; ++j) {
        int tmp;
        MLGetInteger32(mlp, &tmp);
        elems[k][j].list1 = tmp;
        MLGetInteger32(mlp, &tmp);
        elems[k][j].el1 = tmp;
        MLGetInteger32(mlp, &tmp);
        elems[k][j].list2 = tmp;
        MLGetInteger32(mlp, &tmp);
        elems[k][j].el2 = tmp;
      }
    }

    Ret ret = alg.init(nparsin.data(), nparsin.size(), std::move(elems));

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid) || ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_linear_fit(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    std::vector<HornerRatFunPtr> cfun;
    std::vector<HornerRatFunPtr> ifun;
    std::vector<HornerRatFunPtr> lexpr;
    std::vector<HornerHornerRatFunMap> cfunmap;
    std::vector<HornerHornerRatFunMap> ifunmap;
    std::vector<HornerHornerRatFunMap> lexprmap;
    std::vector<std::vector<MPHornerMap>> cfunmpmap;
    std::vector<std::vector<MPHornerMap>> ifunmpmap;
    std::vector<std::vector<MPHornerMap>> lexprmpmap;
    unsigned nvars, sub_nvars, tau_vars;
    int extra_eqs, args_size;
    bool okay = true;

    // input
    MLTestHead( mlp, "List", &args_size);

    // graph and node arguments
    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    auto & ww = session.main_context()->ww;

    // get system data
    MLGetInteger32(mlp, (int*)&tau_vars);
    get_cut_solve_horner_data(mlp, cfun, cfunmap, cfunmpmap,
                              nvars, sub_nvars, ww);
    get_cut_solve_horner_data(mlp, ifun, ifunmap, ifunmpmap,
                              nvars, sub_nvars, ww);
    get_cut_solve_horner_data(mlp, lexpr, lexprmap, lexprmpmap,
                              tau_vars, sub_nvars, ww);

    // extra_eqs options
    MLGetInteger32(mlp, &extra_eqs);

    // define cut object
    typedef AnalyticFitData Data;
    std::unique_ptr<AnalyticFit> algptr(new AnalyticFit());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;

    data->c.swap(cfun);
    data->integr.swap(ifun);
    data->lexpr.swap(lexpr);
    data->cmap.swap(cfunmap);
    data->integrmap.swap(ifunmap);
    data->lexprmap.swap(lexprmap);
    alg.cmpmap.swap(cfunmpmap);
    alg.integrmpmap.swap(ifunmpmap);
    alg.lexprmpmap.swap(lexprmpmap);

    // get subvar sizes
    unsigned max_subvar_size = 0;
    {
      int * lvlist;
      int lvlistsize;
      int ilistsize;

      MLTestHead( mlp, "List", &ilistsize);
      alg.deltavar.reset(new std::unique_ptr<unsigned[]>[ilistsize]);
      for (unsigned i=0; i<(unsigned)ilistsize; ++i) {
        MLGetInteger32List(mlp, &lvlist, &lvlistsize);
        alg.deltavar[i].reset(new unsigned[lvlistsize+1]);
        max_subvar_size = std::max(max_subvar_size,(unsigned)lvlistsize);
        alg.deltavar[i][0] = lvlistsize;
        std::copy(lvlist, lvlist+lvlistsize, alg.deltavar[i].get()+1);
        MLReleaseInteger32List(mlp, lvlist, lvlistsize);
      }

      MLTestHead( mlp, "List", &ilistsize);
      alg.integrvar.reset(new std::unique_ptr<unsigned[]>[ilistsize]);
      for (unsigned i=0; i<(unsigned)ilistsize; ++i) {
        MLGetInteger32List(mlp, &lvlist, &lvlistsize);
        alg.integrvar[i].reset(new unsigned[lvlistsize+1]);
        max_subvar_size = std::max(max_subvar_size,(unsigned)lvlistsize);
        alg.integrvar[i][0] = lvlistsize;
        std::copy(lvlist, lvlist+lvlistsize, alg.integrvar[i].get()+1);
        MLReleaseInteger32List(mlp, lvlist, lvlistsize);
      }

      if (max_subvar_size) {
        data->lsubvar.reset(new UInt[max_subvar_size]);
        alg.max_subvar_size = max_subvar_size;
      }
    }

    // weights
    std::vector<unsigned> w_len;
    w_len.resize(inputnodes.size()-1);
    for (unsigned i=1; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      w_len[i-1] = node->algorithm()->nparsout;
    }

    // needed vars and init
    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);
    Ret ret = alg.init(data->c.size()-1, sub_nvars,
                       data->lexpr.size(), tau_vars,
                       (unsigned*)neededv, needed_size,
                       w_len.data(), w_len.size(), *data,
                       extra_eqs);
    if (ret != SUCCESS) {
      MLNewPacket(mlp);
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    MLReleaseInteger32List(mlp, neededv, needed_size);

    // finished reading input
    MLNewPacket(mlp);

    if (!okay || !session.graph_exists(graphid) || inputnodes.size() < 1) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_numeric_fit(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    std::vector<HornerRatFunPtr> cfun;
    std::vector<HornerRatFunPtr> ifun;
    std::vector<HornerRatFunPtr> lexpr;
    std::vector<MPHornerRatFunMap> cfunmap;
    std::vector<MPHornerRatFunMap> ifunmap;
    std::vector<MPHornerRatFunMap> lexprmap;
    unsigned nvars, tau_vars;
    int extra_eqs, args_size;
    bool okay = true;

    // input
    MLTestHead( mlp, "List", &args_size);

    // graph and node arguments
    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    auto & ww = session.main_context()->ww;

    // get system data
    MLGetInteger32(mlp, (int*)&tau_vars);
    get_cut_solve_numeric_data(mlp, cfun, cfunmap, nvars, ww);
    get_cut_solve_numeric_data(mlp, ifun, ifunmap, nvars, ww);
    get_cut_solve_numeric_data(mlp, lexpr, lexprmap, tau_vars, ww);

    // extra_eqs options
    MLGetInteger32(mlp, &extra_eqs);

    // define cut object
    typedef NumericFitData Data;
    std::unique_ptr<NumericFit> algptr(new NumericFit());
    std::unique_ptr<Data> data(new Data());
    auto & alg = *algptr;

    data->c.swap(cfun);
    data->integr.swap(ifun);
    data->lexpr.swap(lexpr);
    alg.cmap.swap(cfunmap);
    alg.integrmap.swap(ifunmap);
    alg.lexprmap.swap(lexprmap);

    // get subvar sizes
    unsigned max_subvar_size = 0;
    {
      int * lvlist;
      int lvlistsize;
      int ilistsize;

      MLTestHead( mlp, "List", &ilistsize);
      alg.deltavar.reset(new std::unique_ptr<unsigned[]>[ilistsize]);
      for (unsigned i=0; i<(unsigned)ilistsize; ++i) {
        MLGetInteger32List(mlp, &lvlist, &lvlistsize);
        alg.deltavar[i].reset(new unsigned[lvlistsize+1]);
        max_subvar_size = std::max(max_subvar_size,(unsigned)lvlistsize);
        alg.deltavar[i][0] = lvlistsize;
        std::copy(lvlist, lvlist+lvlistsize, alg.deltavar[i].get()+1);
        MLReleaseInteger32List(mlp, lvlist, lvlistsize);
      }

      MLTestHead( mlp, "List", &ilistsize);
      alg.integrvar.reset(new std::unique_ptr<unsigned[]>[ilistsize]);
      for (unsigned i=0; i<(unsigned)ilistsize; ++i) {
        MLGetInteger32List(mlp, &lvlist, &lvlistsize);
        alg.integrvar[i].reset(new unsigned[lvlistsize+1]);
        max_subvar_size = std::max(max_subvar_size,(unsigned)lvlistsize);
        alg.integrvar[i][0] = lvlistsize;
        std::copy(lvlist, lvlist+lvlistsize, alg.integrvar[i].get()+1);
        MLReleaseInteger32List(mlp, lvlist, lvlistsize);
      }

      if (max_subvar_size) {
        data->lsubvar.reset(new UInt[max_subvar_size]);
        alg.max_subvar_size = max_subvar_size;
      }
    }

    // weights
    std::vector<unsigned> w_len;
    w_len.resize(inputnodes.size()-1);
    for (unsigned i=1; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      w_len[i-1] = node->algorithm()->nparsout;
    }

    // needed vars and init
    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);
    Ret ret = alg.init(data->c.size()-1,
                       data->lexpr.size(), tau_vars,
                       (unsigned*)neededv, needed_size,
                       w_len.data(), w_len.size(), *data,
                       extra_eqs);
    if (ret != SUCCESS) {
      MLNewPacket(mlp);
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    MLReleaseInteger32List(mlp, neededv, needed_size);

    // finished reading input
    MLNewPacket(mlp);

    if (!okay || !session.graph_exists(graphid) || inputnodes.size() < 1) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_laurent(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);
    int graphid, subgraphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &subgraphid);

    typedef LaurentExpansionData Data;
    std::unique_ptr<LaurentExpansion> algptr(new LaurentExpansion());
    std::unique_ptr<Data> data(new Data());

    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    Ret ret = algptr->init(session, subgraphid, *data,
                           nparsin.data(), nparsin.size());

    int * order;
    int order_size;
    bool okay = true;
    int sub_nparsout = ret == SUCCESS ? algptr->subgraph()->nparsout : 0;
    MLGetInteger32List(mlp, &order, &order_size);
    if (ret == SUCCESS) {
      if (order_size == sub_nparsout)
        std::copy(order, order+order_size, algptr->order());
      else if (order_size == 1)
        std::fill(algptr->order(), algptr->order()+sub_nparsout, order[0]);
      else
        okay = false;
    } else {
      okay = false;
    }
    MLReleaseInteger32List(mlp, order, order_size);

    int max_deg;
    MLGetInteger32(mlp, &max_deg);
    if (max_deg >= 0)
      algptr->max_degree = max_deg;

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid) || !okay) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_subgraph_fit(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);
    int graphid, subgraphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &subgraphid);

    typedef SubgraphFitData Data;
    std::unique_ptr<SubgraphFit> algptr(new SubgraphFit());
    std::unique_ptr<Data> data(new Data());

    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    int ncoeffs, nsamplevars;
    MLGetInteger32(mlp, &ncoeffs);
    MLGetInteger32(mlp, &nsamplevars);

    int * neededv;
    int needed_size;
    MLGetInteger32List(mlp, &neededv, &needed_size);

    int extra_eqs;
    MLGetInteger32(mlp, &extra_eqs);

    Ret ret = algptr->init(session, subgraphid, *data,
                           nparsin.data(), nparsin.size(),
                           ncoeffs, nsamplevars,
                           (unsigned*)neededv, needed_size,
                           extra_eqs);

    MLReleaseInteger32List(mlp, neededv, needed_size);

    if (ret != SUCCESS) {
      MLNewPacket(mlp);
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_subgraph_multifit(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);
    int graphid, subgraphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &subgraphid);

    typedef SubgraphMultiFitData Data;
    std::unique_ptr<SubgraphMultiFit> algptr(new SubgraphMultiFit());
    std::unique_ptr<Data> data(new Data());

    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    int nsamplevars;
    MLGetInteger32(mlp, &nsamplevars);

    std::vector<std::vector<unsigned>> take;
    MLTestHead(mlp, "List", &nargs);
    take.resize(nargs);
    for (int j=0; j<nargs; ++j) {
      int * takev;
      int take_size;
      MLGetInteger32List(mlp, &takev, &take_size);
      take[j].resize(take_size);
      std::copy(takev, takev+take_size, take[j].data());
      MLReleaseInteger32List(mlp, takev, take_size);
    }

    std::vector<std::vector<unsigned>> needed;
    MLTestHead(mlp, "List", &nargs);
    needed.resize(nargs);
    for (int j=0; j<nargs; ++j) {
      int * neededv;
      int needed_size;
      MLGetInteger32List(mlp, &neededv, &needed_size);
      needed[j].resize(needed_size);
      std::copy(neededv, neededv+needed_size, needed[j].data());
      MLReleaseInteger32List(mlp, neededv, needed_size);
    }

    int extra_eqs;
    MLGetInteger32(mlp, &extra_eqs);

    Ret ret = algptr->init(session, subgraphid, *data,
                           nparsin.data(), nparsin.size(),
                           std::move(take), nsamplevars,
                           needed, extra_eqs);

    if (ret != SUCCESS) {
      MLNewPacket(mlp);
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_subgraph_rec(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);
    int graphid, subgraphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &subgraphid);

    typedef SubgraphRecData Data;
    std::unique_ptr<SubgraphRec> algptr(new SubgraphRec());
    std::unique_ptr<Data> data(new Data());

    int nrec=0, shiftvars=0;
    MLGetInteger32(mlp, &nrec);
    MLGetInteger32(mlp, &shiftvars);
    //ReconstructionOptions opt = get_rec_opt(mlp);

    Ret ret = FAILED;

    if (inputnodes.size()==1) {
      Node * node = session.node(graphid, inputnodes[0]);
      if (node) {
        unsigned npars = node->algorithm()->nparsout;
        ret = algptr->init(session, subgraphid, *data,
                           npars, nrec, shiftvars);
      }
    }

    if (ret != SUCCESS) {
      MLNewPacket(mlp);
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    //algptr->setOptions(opt);

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_default_maxdeg(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int zero;
    MLTestHead( mlp, "List", &zero);
    MLNewPacket(mlp);
    MLPutInteger32(mlp, RatFunReconstruction::DEFAULT_MAX_DEG);

    return LIBRARY_NO_ERROR;
  }

  int fflowml_dump_degrees(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs, gid;
    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &gid);
    const char * filename;
    MLGetString(mlp, &filename);
    Ret ret = session.dump_degrees(gid, filename);
    MLReleaseString(mlp, filename);
    MLNewPacket(mlp);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutSymbol(mlp, "Null");
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_load_degrees(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs, gid;
    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &gid);
    const char * filename;
    MLGetString(mlp, &filename);
    Ret ret = session.load_degrees(gid, filename);
    MLReleaseString(mlp, filename);
    MLNewPacket(mlp);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutSymbol(mlp, "Null");
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_dump_sample_points(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs, gid;
    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &gid);
    const char * filename;
    MLGetString(mlp, &filename);
    ReconstructionOptions opt = get_rec_opt(mlp);
    Ret ret = session.dump_sample_points(gid, opt, filename);
    MLReleaseString(mlp, filename);
    MLNewPacket(mlp);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutSymbol(mlp, "Null");
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_dump_evaluations(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs, gid;
    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &gid);
    const char * filename;
    MLGetString(mlp, &filename);
    Ret ret = session.dump_evaluations(gid, filename);
    MLReleaseString(mlp, filename);
    MLNewPacket(mlp);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutSymbol(mlp, "Null");
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_load_evaluations(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs, gid, nfiles;
    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &gid);
    std::vector<const char *> filename;
    MLTestHead(mlp, "List", &nfiles);
    filename.resize(nfiles);
    for (int j=0; j<nfiles; ++j)
      MLGetString(mlp, &(filename[j]));

    Ret ret = session.load_evaluations(gid, filename.data(), nfiles);

    for (int j=0; j<nfiles; ++j)
      MLReleaseString(mlp, filename[j]);

    MLNewPacket(mlp);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutSymbol(mlp, "Null");
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_samples_file_size(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLTestHead( mlp, "List", &nargs);

    const char * filename;
    MLGetString(mlp, &filename);
    UInt size = samples_file_size(filename);
    MLReleaseString(mlp, filename);
    MLNewPacket(mlp);

    if (size == FAILED) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger64(mlp, size);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_npars_from_degree_info(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLTestHead( mlp, "List", &nargs);

    const char * filename;
    unsigned out[2];
    MLGetString(mlp, &filename);
    Ret ret = algorithm_npars_from_degree_info(filename, out[0], out[1]);
    MLReleaseString(mlp, filename);
    MLNewPacket(mlp);

    if (ret == FAILED) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32List(mlp, (int*)(out), 2);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_graph_evaluate(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs, gid;
    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &gid);

    mlint64 * xin;
    int n_xin;
    std::unique_ptr<UInt[]> x;
    MLGetInteger64List(mlp, &xin, &n_xin);
    x.reset(new UInt[n_xin]);
    std::copy(xin, xin+n_xin, x.get());
    MLReleaseInteger64List(mlp, xin, n_xin);

    int prime_no;
    MLGetInteger32(mlp, &prime_no);

    MLNewPacket(mlp);

    Graph * g = session.graph(gid);
    if (!session.graph_can_be_evaluated(gid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    if (prime_no < 0)
      prime_no = 0;

    unsigned nparsin = g->nparsin[0];
    if (nparsin != unsigned(n_xin) ||
        unsigned(prime_no) > BIG_UINT_PRIMES_SIZE)  {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Mod mod(BIG_UINT_PRIMES[prime_no]);
    unsigned nparsout = g->nparsout;
    std::unique_ptr<UInt[]> xout;
    xout.reset(new UInt[nparsout]);

    {
      const UInt * xxin = x.get();
      Context * ctxt = session.main_context();
      Ret ret = g->evaluate(ctxt, &xxin, mod,
                            ctxt->graph_data(gid),
                            xout.get());
      if (ret != SUCCESS)  {
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
    }

    MLPutInteger64List(mlp, (mlint64*)(xout.get()), nparsout);

    return LIBRARY_NO_ERROR;
  }





  int fflowml_graph_evaluate_list(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    const UInt max_primes = BIG_UINT_PRIMES_SIZE;

    int nargs, gid;
    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &gid);

    Graph * g = session.graph(gid);

    if (!session.graph_can_be_evaluated(gid)) {
      MLNewPacket(mlp);
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    unsigned nparsin = g->nparsin[0];
    unsigned nparsout = g->nparsout;

    // nthreads
    int nthreads;
    MLGetInteger32(mlp, &nthreads);
    if (nthreads < 0)
      nthreads = 0;

    // default prime
    int prime_no;
    MLGetInteger32(mlp, &prime_no);

    if (prime_no < 0)
      prime_no = 0;
    if (unsigned(prime_no) >= BIG_UINT_PRIMES_SIZE)  {
      MLNewPacket(mlp);
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    mlint64 * xin;
    int n_points=0, n_xin=0;
    std::unique_ptr<UInt[]> x;
    SamplePointsVector points;
    MLTestHead( mlp, "List", &n_points);
    points.reserve(n_points);
    for (int j=0; j<n_points; ++j) {
      MLGetInteger64List(mlp, &xin, &n_xin);
      if ((unsigned(n_xin) != nparsin && unsigned(n_xin) != nparsin+1) ||
          (unsigned(n_xin) == nparsin+1 && UInt(xin[nparsin]) >= max_primes)) {
        MLNewPacket(mlp);
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      x.reset(new UInt[nparsin+1]);
      std::copy(xin, xin+n_xin, x.get());
      if (unsigned(n_xin) == nparsin)
        x[nparsin] = BIG_UINT_PRIMES[prime_no];
      else
        x[nparsin] = BIG_UINT_PRIMES[x[nparsin]];
      points.push_back(std::move(x));
      MLReleaseInteger64List(mlp, xin, n_xin);
    }

    MLNewPacket(mlp);

    SamplePointsVector xout;
    {
      Ret ret = session.evaluate_list(gid, points, xout, nthreads);
      if (ret != SUCCESS)  {
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
    }

    MLPutFunction(mlp, "List", n_points);
    for (int j=0; j<n_points; ++j) {
      if (xout[j][0] == FAILED)
        MLPutSymbol(mlp, "$Failed");
      else
        MLPutInteger64List(mlp, (mlint64*)(xout[j].get()), nparsout);
    }

    return LIBRARY_NO_ERROR;
  }



  int fflowml_alg_sparse_mat_mul(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args, nrows1, ncols1, ncols2;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    std::unique_ptr<SparseMatrixMul> algptr(new SparseMatrixMul());
    auto & alg = *algptr;

    MLGetInteger32(mlp, &nrows1);
    MLGetInteger32(mlp, &ncols1);
    MLGetInteger32(mlp, &ncols2);

    int nrows_l1 = 0;
    MLTestHead(mlp, "List", &nrows_l1);
    alg.row1.resize(nrows_l1);
    for (int j=0; j<nrows_l1; ++j) {
      int * cols;
      int len;
      MLGetInteger32List(mlp, &cols, &len);
      alg.row1[j].size = len;
      alg.row1[j].cols.reset(new unsigned[len]);
      unsigned * rcols = alg.row1[j].cols.get();
      std::copy(cols, cols+len, rcols);
      MLReleaseInteger32List(mlp, cols, len);
    }

    int nrows_l2 = 0;
    MLTestHead(mlp, "List", &nrows_l2);
    alg.row2.resize(nrows_l2);
    for (int j=0; j<nrows_l2; ++j) {
      int * cols;
      int len;
      MLGetInteger32List(mlp, &cols, &len);
      alg.row2[j].size = len;
      alg.row2[j].cols.reset(new unsigned[len]);
      unsigned * rcols = alg.row2[j].cols.get();
      std::copy(cols, cols+len, rcols);
      MLReleaseInteger32List(mlp, cols, len);
    }

    MLNewPacket(mlp);

    alg.init(nrows1, ncols1, ncols2);

    if (!session.graph_exists(graphid) || inputnodes.size() != 2) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), nullptr,
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_count_sample_points(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int two;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &two);

    int id, nthreads;
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nthreads);
    if (nthreads < 0)
      nthreads = 0;
    ReconstructionOptions opt = get_rec_opt(mlp);

    std::string file;
    {
      const char * filename;
      MLGetString(mlp, &filename);
      file = std::string(filename);
      MLReleaseString(mlp, filename);
    }

    MLNewPacket(mlp);

    Graph * graph = session.graph(id);
    if (graph == nullptr) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }
    unsigned nparsout = graph->nparsout;
    std::unique_ptr<unsigned[]> count(new unsigned[nparsout]);
    const char * cfile = file.length() ? file.c_str() : nullptr;
    Ret ret = session.count_sample_points(id, opt, count.get(), nthreads,
                                          cfile);
    if (ret == FAILED) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutFunction(mlp, "List", 2);
      MLPutInteger64(mlp, ret);
      MLPutInteger32List(mlp, (int*)(count.get()), nparsout);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_rat_rec(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int two, nints;
    MLNewPacket(mlp);
    MLTestHead(mlp, "List", &two);

    MLTestHead(mlp, "List", &nints);
    std::vector<MPInt> vec;
    vec.resize(nints);
    for (int j=0; j<nints; ++j)
      get_integer(mlp,vec[j]);

    MPInt p;
    get_integer(mlp,p);

    MLNewPacket(mlp);

    std::vector<MPRational> qvec;
    qvec.resize(nints);
    for (int j=0; j<nints; ++j)
      rat_rec(vec[j], p, qvec[j]);

    void (*gmpfreefunc) (void *, size_t);
    mp_get_memory_functions(0, 0, &gmpfreefunc);

    MLPutFunction(mlp, "List", nints);
    for (int i=0; i<nints; ++i)
      put_mprat(mlp, qvec[i], gmpfreefunc);

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_ratexpr_eval(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int n_args, nfunctions;
    int nparsin=0;
    MLNewPacket(mlp);
    MLTestHead( mlp, "List", &n_args);

    int graphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);

    typedef AnalyticExpressionData Data;
    std::unique_ptr<AnalyticExpression> algptr(new AnalyticExpression());
    std::unique_ptr<Data> data(new Data());
    AnalyticExpression & alg = *algptr;

    MLGetInteger32(mlp, &nparsin);

    MLTestHead(mlp, "List", &nfunctions);
    std::vector<std::vector<Instruction>> bytecode(nfunctions);

    for (int i=0; i<nfunctions; ++i) {
      short * datap;
      int countp;
      MLGetInteger16List(mlp, &datap, &countp);
      bytecode[i].resize(countp);
      std::copy(datap, datap+countp, bytecode[i].data());
      MLReleaseInteger16List(mlp, datap, countp);
    }

    int nnumbers;
    MLTestHead( mlp, "List", &nnumbers);
    std::vector<MPRational> numbers(nnumbers);
    for (int i=0; i<nnumbers; ++i)
      get_rational(mlp, numbers[i]);

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid) || inputnodes.size() != 1) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);

    alg.init(nparsin, std::move(bytecode), std::move(numbers), *data);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }


  int fflowml_alg_cached_subgraph(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);
    int graphid, subgraphid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &subgraphid);

    typedef CachedSubGraphData Data;
    std::unique_ptr<CachedSubGraph> algptr(new CachedSubGraph());
    std::unique_ptr<Data> data(new Data());

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    std::vector<unsigned> nparsin(inputnodes.size());
    for (unsigned i=0; i<inputnodes.size(); ++i) {
      Node * node = session.node(graphid, inputnodes[i]);
      if (!node) {
        MLPutSymbol(mlp, "$Failed");
        return LIBRARY_NO_ERROR;
      }
      nparsin[i] = node->algorithm()->nparsout;
    }

    Ret ret = algptr->init(session, subgraphid, *data,
                           nparsin.data(), nparsin.size());

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_cached_from_subgraph(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int nargs;
    MLNewPacket(mlp);

    MLTestHead(mlp, "List", &nargs);
    int graphid, subgraphid, subnodeid;
    std::vector<unsigned> inputnodes;
    MLGetInteger32(mlp, &graphid);
    get_input_nodes(mlp, inputnodes);
    MLGetInteger32(mlp, &subgraphid);
    MLGetInteger32(mlp, &subnodeid);

    typedef CachedFromSubGraphData Data;
    std::unique_ptr<CachedFromSubGraph> algptr(new CachedFromSubGraph());
    std::unique_ptr<Data> data(new Data());

    MLNewPacket(mlp);

    if (!session.graph_exists(graphid)) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Ret ret = algptr->init(&session, subgraphid, subnodeid);

    if (ret != SUCCESS) {
      MLPutSymbol(mlp, "$Failed");
      return LIBRARY_NO_ERROR;
    }

    Graph * graph = session.graph(graphid);
    unsigned id = graph->new_node(std::move(algptr), std::move(data),
                                  inputnodes.data());

    if (id == ALG_NO_ID) {
      MLPutSymbol(mlp, "$Failed");
    } else {
      MLPutInteger32(mlp, id);
    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_cached_subgraph_merge(WolframLibraryData libData, MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int id, nodeid, nargs;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    Algorithm * alg = session.algorithm(id, nodeid);
    MLNewPacket(mlp);

    if (!alg) {

      MLPutSymbol(mlp, "$Failed");

    } else if (dynamic_cast<CachedSubGraph *>(alg)) {

      unsigned nsubctxt = session.subcontext_size();
      std::vector<std::unique_ptr<UIntCache>> caches;
      caches.reserve(nsubctxt + 1);

      CachedSubGraph & subg = *static_cast<CachedSubGraph *>(alg);

      {
        auto * data = static_cast<CachedSubGraphData*>(session.alg_data(id,nodeid));
        if (data && data->cache())
          caches.push_back(std::move(data->move_cache()));
      }

      for (unsigned j=0; j<nsubctxt; ++j) {
        auto * data = static_cast<CachedSubGraphData*>(session.subctxt_alg_data(j,id,nodeid));
        if (data && data->cache())
          caches.push_back(std::move(data->move_cache()));
      }

      subg.merge_caches(caches.data(), caches.size());

      MLPutSymbol(mlp, "Null");

    } else {

      MLPutSymbol(mlp, "$Failed");

    }

    return LIBRARY_NO_ERROR;
  }

  int fflowml_alg_cached_subgraph_default_subcache_size(WolframLibraryData libData,
                                                        MLINK mlp)
  {
    (void)(libData);
    FFLOWML_SET_DBGPRINT();

    int id, nodeid, nargs, cachesize;
    MLNewPacket(mlp);

    MLTestHead( mlp, "List", &nargs);
    MLGetInteger32(mlp, &id);
    MLGetInteger32(mlp, &nodeid);
    MLGetInteger32(mlp, &cachesize);
    Algorithm * alg = session.algorithm(id, nodeid);
    MLNewPacket(mlp);

    if (!alg) {

      MLPutSymbol(mlp, "$Failed");

    } else if (dynamic_cast<CachedSubGraph *>(alg)) {

      CachedSubGraph & subg = *static_cast<CachedSubGraph *>(alg);
      subg.set_default_subcache_size(cachesize);
      MLPutSymbol(mlp, "Null");

    } else {

      MLPutSymbol(mlp, "$Failed");

    }

    return LIBRARY_NO_ERROR;
  }

} // extern "C"
