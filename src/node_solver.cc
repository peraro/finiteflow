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

      for (unsigned j=0; j<n_cols; ++j, ++r, ++info)
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



  void NodeSparseSolver::init_node_solver(std::size_t neqs, std::size_t nvars,
                                          unsigned * needed_vars,
                                          unsigned needed_size,
                                          NodeSparseSolverData & data)

  {
    nparsin.resize(1);
    unsigned nin = 0;
    for (unsigned j=0; j<neqs; ++j)
      nin += rinfo[j].size;
    nparsin[0] = nin;
    init(neqs, nvars, needed_vars, needed_size, data);
  }

  Ret NodeSparseSolver::fill_matrix(Context *,
                                    std::size_t n_rows,
                                    const std::size_t rows[],
                                    AlgInput xin[], Mod,
                                    AlgorithmData *,
                                    SparseMatrix & m) const
  {
    const UInt * xi = xin[0];

    for (unsigned i=0; i<n_rows; ++i) {

      SparseMatrixRow & r = m.row(i);

      const std::size_t row_start = rinfo[rows[i]].start;
      const std::size_t row_size = rinfo[rows[i]].size;
      const unsigned * cols = rinfo[rows[i]].cols.get();
      const UInt * f = xi + row_start;
      const UInt * fend = f + row_size;

      r.resize(row_size);
      unsigned j=0;

      for (; f<fend; ++f, ++j) {
        UInt res = (*f);
        unsigned col = cols[j];
        if (res == FAILED || res == 0)
          return FAILED;
        r.el(j).col = col;
        r.el(j).val = res;
      }
      r.el(j).col = SparseMatrixRow::END;
    }

    return SUCCESS;
  }

  AlgorithmData::Ptr
  NodeSparseSolver::clone_data(const AlgorithmData * datain) const
  {
    std::unique_ptr<NodeSparseSolverData> ptr(new NodeSparseSolverData());
    SparseLinearSolver::copy_data(datain, ptr.get());
    return std::move(ptr);
  }


} // namespace fflow
