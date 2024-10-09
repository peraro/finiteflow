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
    std::size_t n_rows = rinfo.size();
    std::unique_ptr<bool[]> needed(new bool[n_rows]());

    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();

    for (unsigned i=0; i<nindepeqs; ++i)
      needed[ieq[i]] = true;

    for (unsigned i=0; i<n_rows; ++i)
      if (!needed[i]) {
        rinfo[i] = RowInfo();
        c[i].reset(nullptr);
      }
  }

  Ret NumericSparseSolver::fill_matrix(Context *,
                                       std::size_t n_rows,
                                       const std::size_t rows[],
                                       AlgInput[], Mod mod,
                                       AlgorithmData *,
                                       SparseMatrix & m) const
  {
    const MPInt mpmod(mod.n());
    MPInt mpres;

    for (unsigned i=0; i<n_rows; ++i) {

      SparseMatrixRow & r = m.row(i);

      const std::size_t row_size = rinfo[rows[i]].size;
      const unsigned * cols = rinfo[rows[i]].cols.get();
      const MPRational * f = c[rows[i]].get();
      const MPRational * fend = f + row_size;

      r.resize(row_size);
      unsigned j=0;

      for (; f<fend; ++f, ++j) {
        rat_mod(*f, mpmod, mpres);
        unsigned col = cols[j];
        r.el(j).col = col;
        r.el(j).val.set(mpres.to_uint());
      }
      r.el(j).col = SparseMatrixRow::END;
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
