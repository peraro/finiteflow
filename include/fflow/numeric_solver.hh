#ifndef FFLOW_ALG_NUMERIC_SOLVER_HH
#define FFLOW_ALG_NUMERIC_SOLVER_HH

#include<fflow/alg_linear_solver.hh>
#include<fflow/mp_functions.hh>

namespace fflow {

  class NumericDenseSolver;
  class NumericSparseSolver;


  struct NumericDenseSolverData : public DenseLinearSolverData {};

  class NumericDenseSolver : public DenseLinearSolver {
  public:
    virtual Ret fill_matrix(Context * ctxt,
                            std::size_t n_rows, const std::size_t rows[],
                            AlgInput xi[], Mod mod,
                            AlgorithmData * data,
                            MatrixView & m) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  public:
    DynamicMatrixT<MPRational> c;
  };


  struct NumericSparseSolverData : public SparseLinearSolverData {};

  class NumericSparseSolver : public SparseLinearSolver {
  public:
    void delete_unneeded_eqs(AlgorithmData * data);

    virtual Ret fill_matrix(Context * ctxt,
                            unsigned n_rows, const unsigned rows[],
                            AlgInput xi[], Mod mod,
                            AlgorithmData * data,
                            SparseMatrix & m) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  private:
    void reset_flags_(AlgorithmData * data);

  public:
    struct RowInfo {
      std::size_t size;
      std::unique_ptr<unsigned[]> cols;
      std::unique_ptr<std::size_t[]> idx;
    };

    std::vector<RowInfo> rinfo;
    std::vector<MPRational> c;
  };

} // namespace fflow


#endif // FFLOW_ALG_NUMERIC_SOLVER_HH
