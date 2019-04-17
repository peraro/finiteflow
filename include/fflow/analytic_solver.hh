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
    std::vector<std::unique_ptr<HornerRatFunPtr[]>> c;

  private:
    friend class AnalyticSparseSolver;
    void reset_mod_(const AnalyticSparseSolver & ls, Mod mod);

  private:
    std::unique_ptr<UInt[]> xp_ = nullptr;
    UInt this_mod_ = 0;
  };

  class AnalyticSparseSolver : public SparseLinearSolver {
  public:
    void delete_unneeded_eqs(AlgorithmData * data);

    virtual Ret fill_matrix(Context * ctxt,
                            std::size_t n_rows, const std::size_t rows[],
                            AlgInput xi[], Mod mod,
                            AlgorithmData * data,
                            SparseMatrix & m) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  public:
    struct RowInfo {
      std::size_t size;
      std::unique_ptr<unsigned[]> cols;
    };

    std::vector<RowInfo> rinfo;
    std::vector<std::unique_ptr<MPHornerRatFunMap[]>> cmap;
  };

} // namespace fflow


#endif // FFLOW_ANALYTIC_Solver_HH
