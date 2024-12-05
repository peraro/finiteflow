#include <fflow/numeric_solver.hh>
#include <fflow/mp_gcd.hh>

namespace fflow {

  Ret NumericDenseSolver::fill_matrix(Context *,
                                      std::size_t n_rows,
                                      const std::size_t rows[],
                                      AlgInput[], Mod mod,
                                      AlgorithmData *,
                                      MatrixView & m) const
  {
    const MPInt mpmod(mod.n());
    MPInt mpres;
    const unsigned n_cols = c.columns();

    for (unsigned i=0; i<n_rows; ++i) {

      const MPRational * f = c.row(rows[i]);
      const MPRational * fend = f + n_cols;
      UInt * r = m.row(i);
      const DenseLinearSolver::flag_t * info = xinfo();

      for (; f<fend; ++f, ++r, ++info)
        if ((*info) & LSVar::IS_NON_ZERO) {
          rat_mod((*f), mpmod, mpres);
          *r = mpres.to_uint();
        }
    }

    return SUCCESS;
  }

  AlgorithmData::Ptr
  NumericDenseSolver::clone_data(const AlgorithmData * datain) const
  {
    std::unique_ptr<NumericDenseSolverData> ptr(new NumericDenseSolverData());
    DenseLinearSolver::copy_data(datain, ptr.get());
    return std::move(ptr);
  }


  void NumericSparseSolver::delete_unneeded_eqs(AlgorithmData *)
  {
    const SparseLinearSolver::flag_t * info = xinfo();
    const std::size_t n_rows = rinfo.size();
    const std::size_t c_len = c.size();
    std::unique_ptr<bool[]> needed(new bool[n_rows]());
    std::unique_ptr<bool[]> needed_ccs(new bool[c_len]());

    std::size_t nindepeqs = n_indep_eqs();
    const unsigned * ieq = indep_eqs();

    // Mark all needed functions
    for (unsigned ind=0; ind<nindepeqs; ++ind) {
      std::size_t i = ieq[ind];
      needed[i] = true;
      std::size_t row_size = rinfo[i].size;
      const std::size_t * idx = rinfo[i].idx.get();
      const unsigned * cols = rinfo[i].cols.get();
      const unsigned * cols_end = cols + row_size;
      for (; cols<cols_end; ++idx, ++cols)
        if (info[*cols] & LSVar::IS_NON_ZERO)
          needed_ccs[*idx] = true;
    }

    // Delete what is not needed anymore
    for (unsigned j=0; j<c_len; ++j)
      if (!needed_ccs[j])
        c[j] = MPRational();

    // Clean up row info
    for (unsigned i=0; i<n_rows; ++i)
      if (!needed[i])
        rinfo[i] = RowInfo();
  }

  Ret NumericSparseSolver::fill_matrix(Context *,
                                       unsigned n_rows,
                                       const unsigned rows[],
                                       AlgInput[], Mod mod,
                                       AlgorithmData *,
                                       SparseMatrix & m) const
  {
    const SparseLinearSolver::flag_t * info = xinfo();
    const MPInt mpmod(mod.n());
    MPInt mpres;

    for (unsigned i=0; i<n_rows; ++i) {

      SparseMatrixRow & r = m.row(i);

      const std::size_t row_size = rinfo[rows[i]].size;
      const unsigned * cols = rinfo[rows[i]].cols.get();
      const unsigned * cols_end = cols + row_size;
      const std::size_t * idx = rinfo[rows[i]].idx.get();

      r.resize(row_size);
      unsigned oj=0;

      // NOTE: rat_mod is evaluated many times here for the same c[j]
      // in the (very likely) case there are duplicates.
      for (; cols<cols_end; ++cols, ++idx) {
        unsigned col = *cols;
        if (info[col] & LSVar::IS_NON_ZERO) {
          rat_mod(c[*idx], mpmod, mpres);
          r.el(oj).col = col;
          r.el(oj).val.set(mpres.to_uint());
          ++oj;
        }
      }
      r.el(oj).col = SparseMatrixRow::END;
      r.resize(oj);
    }

    return SUCCESS;
  }

  AlgorithmData::Ptr
  NumericSparseSolver::clone_data(const AlgorithmData * datain) const
  {
    std::unique_ptr<NumericSparseSolverData> ptr(new NumericSparseSolverData());
    copy_data(datain, ptr.get());
    return std::move(ptr);
  }

} // namespace fflow
