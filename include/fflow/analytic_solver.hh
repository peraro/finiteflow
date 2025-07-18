#ifndef FFLOW_ANALYTIC_LS_HH
#define FFLOW_ANALYTIC_LS_HH

#include<fflow/alg_linear_solver.hh>
#include<fflow/mp_functions.hh>

namespace fflow {

  class AnalyticDenseSolver;
  class AnalyticSparseSolver;


  struct AnalyticDenseSolverData : public DenseLinearSolverData {
  public:
    DynamicMatrixT<HornerRatFunPtr> c;
  private:
    friend class AnalyticDenseSolver;
    void reset_mod_(const AnalyticDenseSolver & ls, Mod mod);
  private:
    std::unique_ptr<UInt[]> xp_ = nullptr;
    UInt this_mod_ = 0;
  };

  class AnalyticDenseSolver : public DenseLinearSolver {
  public:
    virtual Ret fill_matrix(Context * ctxt,
                            std::size_t n_rows, const std::size_t rows[],
                            AlgInput xi[], Mod mod,
                            AlgorithmData * data,
                            MatrixView & m) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  public:
    DynamicMatrixT<MPHornerRatFunMap> cmap;
  };


  struct AnalyticSparseSolverData : public SparseLinearSolverData {
  public:
    std::vector<HornerRatFunPtr> c;

  private:
    friend class AnalyticSparseSolver;
    void reset_mod_(const AnalyticSparseSolver & ls, Mod mod);
    void reset_evals_(const AnalyticSparseSolver & ls);

  private:
    std::unique_ptr<UInt[]> xp_ = nullptr;
    std::unique_ptr<UInt[]> evals_ = nullptr;
    UInt this_mod_ = 0;
  };

  class AnalyticSparseSolver : public SparseLinearSolver {
  public:
    void delete_unneeded_eqs(AlgorithmData * data);

    virtual Ret fill_matrix(Context * ctxt,
                            unsigned n_rows, const unsigned rows[],
                            AlgInput xi[], Mod mod,
                            AlgorithmData * data,
                            SparseMatrix & m) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  public:
    struct RowInfo {
      std::size_t size;
      std::unique_ptr<unsigned[]> cols;
      std::unique_ptr<std::size_t[]> idx;
    };

    std::vector<RowInfo> rinfo;
    std::vector<MPHornerRatFunMap> cmap;
  };

  inline
  void AnalyticSparseSolverData::reset_evals_(const AnalyticSparseSolver & ls)
  {
    std::size_t len = ls.cmap.size();
    if (!evals_)
      evals_.reset(new UInt[len]);
    for (unsigned j=0; j<len; ++j)
      evals_[j] = MISSING_SAMPLES;
  }


  class AnalyticSparseSolverEx;
  struct AnalyticSparseSolverExData : public SparseLinearSolverData {
  public:
    std::vector<HornerRatFunPtr> c;

  private:
    friend class AnalyticSparseSolverEx;
    void reset_mod_(const AnalyticSparseSolverEx & ls, Mod mod);
    void reset_evals_(const AnalyticSparseSolverEx & ls);

  private:
    std::unique_ptr<UInt[]> xp_ = nullptr;
    std::unique_ptr<UInt[]> evals_ = nullptr;
    UInt this_mod_ = 0;
  };

  class AnalyticSparseSolverEx : public SparseLinearSolver {
  public:
    void delete_unneeded_eqs(AlgorithmData * data);

    virtual Ret fill_matrix(Context * ctxt,
                            unsigned n_rows, const unsigned rows[],
                            AlgInput xi[], Mod mod,
                            AlgorithmData * data,
                            SparseMatrix & m) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

    virtual Ret learn(Context * ctxt,
                      AlgInput xin[], Mod mod, AlgorithmData * data) override;

  private:

    Ret learn_zero_ccs_(Context * ctxt,
                        AlgInput xin[], Mod mod, AlgorithmData * data);

  public:

    struct Weight {
      unsigned node, el;
      std::size_t idx;
    };

    struct Weights {
      unsigned size;
      std::unique_ptr<Weight[]> w;

      void from(Weights && w2)
      {
        size = w2.size;
        w = std::move(w2.w);
      }
    };

    struct RowInfo {
      std::size_t size;
      std::unique_ptr<unsigned[]> cols;
      std::unique_ptr<Weights[]> w;
    };

    std::vector<RowInfo> rinfo;
    std::vector<MPHornerRatFunMap> cmap;
    bool learned_ccs_ = false;
  };

  inline void
  AnalyticSparseSolverExData::reset_evals_(const AnalyticSparseSolverEx & ls)
  {
    std::size_t len = ls.cmap.size();
    if (!evals_)
      evals_.reset(new UInt[len]);
    for (unsigned j=0; j<len; ++j)
      evals_[j] = MISSING_SAMPLES;
  }


} // namespace fflow


#endif // FFLOW_ANALYTIC_Solver_HH
