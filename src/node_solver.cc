#include <fflow/node_solver.hh>
#include <fflow/mp_gcd.hh>

namespace fflow {

  Ret NodeDenseSolver::fill_matrix(Context *,
                                   std::size_t n_rows,
                                   const std::size_t rows[],
                                   AlgInput xin[], Mod,
                                   AlgorithmData *,
                                   MatrixView & m) const
  {
    AlgInput x = xin[0];
    const unsigned n_cols = DenseLinearSolver::nvars()+1;

    for (unsigned i=0; i<n_rows; ++i) {

      std::size_t row_i = rows[i];
      UInt * r = m.row(i);
      const DenseLinearSolver::flag_t * info = xinfo();

      for (unsigned j=0; j<n_cols; ++j, ++r)
        if ((*info) & LSVar::IS_NON_ZERO) {
          *r = x[row_i*n_cols + j];
        }
    }

    return SUCCESS;
  }

  AlgorithmData::Ptr
  NodeDenseSolver::clone_data(const AlgorithmData * datain) const
  {
    std::unique_ptr<NodeDenseSolverData> ptr(new NodeDenseSolverData());
    DenseLinearSolver::copy_data(datain, ptr.get());
    return std::move(ptr);
  }


} // namespace fflow
