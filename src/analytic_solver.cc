#include <fflow/analytic_solver.hh>
#include <fflow/graph.hh>

namespace fflow {

  void AnalyticDenseSolverData::reset_mod_(const AnalyticDenseSolver & ls,
                                           Mod mod)
  {
    this_mod_ = mod.n();
    if (xp_.get() == nullptr)
      xp_.reset(new UInt[ls.nparsin[0]]);
    const unsigned n_cols = c.columns();
    std::size_t nindepeqs = ls.n_indep_eqs();
    const std::size_t * ieq = ls.indep_eqs();
    for (unsigned i=0; i<nindepeqs; ++i)
      for (unsigned j=0; j<n_cols; ++j)
        horner_ratfunmap_mprational(ls.cmap(ieq[i],j), mod, c(ieq[i],j));
  }

  Ret AnalyticDenseSolver::fill_matrix(Context * ctxt,
                                       std::size_t n_rows,
                                       const std::size_t rows[],
                                       AlgInput xin[], Mod mod,
                                       AlgorithmData * datain,
                                       MatrixView & m) const
  {
    auto & data = *static_cast<AnalyticDenseSolverData*>(datain);
    const unsigned nin = nparsin[0];

    if (mod.n() != data.this_mod_)
      data.reset_mod_(*this, mod);

    const UInt * xi = xin[0];
    UInt * xp = data.xp_.get();
    precomp_array_mul_shoup(xi, nin, mod, xp);

    UInt * ww = ctxt->ww.get();
    const unsigned n_cols = data.c.columns();

    for (unsigned i=0; i<n_rows; ++i) {

      const HornerRatFunPtr * f = data.c.row(rows[i]);
      const HornerRatFunPtr * fend = f + n_cols;
      UInt * r = m.row(i);
      const DenseLinearSolver::flag_t * info = xinfo();

      for (; f<fend; ++f, ++r, ++info)
        if ((*info) & LSVar::IS_NON_ZERO) {
          UInt res = (*f).eval(nin, ww, xi, xp, mod);
          if (res == FAILED)
            return FAILED;
          *r = res;
        }
    }

    return SUCCESS;
  }

  AlgorithmData::Ptr
  AnalyticDenseSolver::clone_data(const AlgorithmData * datain) const
  {
    const auto & data = *static_cast<const AnalyticDenseSolverData*>(datain);

    std::unique_ptr<AnalyticDenseSolverData> ptr(new AnalyticDenseSolverData());

    DenseLinearSolver::copy_data(datain, ptr.get());

    auto & newalg = *ptr;
    newalg.c.resize(data.c.rows(), data.c.columns());

    std::size_t nindepeqs = n_indep_eqs();
    const std::size_t * ieq = indep_eqs();
    for (unsigned i=0; i<nindepeqs; ++i)
      for (unsigned j=0; j<data.c.columns(); ++j)
        horner_ratfun_clone(data.c(ieq[i],j), nparsin[0], newalg.c(ieq[i],j));

    return std::move(ptr);
  }


  void AnalyticSparseSolverData::reset_mod_(const AnalyticSparseSolver & ls,
                                            Mod mod)
  {
    this_mod_ = mod.n();
    if (xp_.get() == nullptr)
      xp_.reset(new UInt[ls.nparsin[0]]);
    reset_evals_(ls);
    const SparseLinearSolver::flag_t * info = ls.xinfo();
    std::size_t nindepeqs = ls.n_indep_eqs();
    const unsigned * ieq = ls.indep_eqs();
    for (unsigned ind=0; ind<nindepeqs; ++ind) {
      std::size_t i = ieq[ind];
      std::size_t row_size = ls.rinfo[i].size;
      const std::size_t * idx = ls.rinfo[i].idx.get();
      const unsigned * cols = ls.rinfo[i].cols.get();
      const unsigned * cols_end = cols + row_size;
      for (; cols<cols_end; ++idx, ++cols)
        if ((info[*cols] & LSVar::IS_NON_ZERO)
            && evals_[*idx] == MISSING_SAMPLES) {
          horner_ratfunmap_mprational(ls.cmap[*idx], mod, c[*idx]);
          evals_[*idx] = 1;
        }
    }
  }

  void AnalyticSparseSolver::delete_unneeded_eqs(AlgorithmData * datain)
  {
    auto & data = *static_cast<AnalyticSparseSolverData*>(datain);

    data.reset_evals_(*this);

    const SparseLinearSolver::flag_t * info = xinfo();
    std::size_t n_rows = rinfo.size();
    std::unique_ptr<bool[]> needed(new bool[n_rows]());

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
          data.evals_[*idx] = 1;
    }

    // Delete what is not needed anymore
    std::size_t len = cmap.size();
    for (unsigned j=0; j<len; ++j)
      if (data.evals_[j] == MISSING_SAMPLES) {
        data.c[j] = HornerRatFunPtr();
        cmap[j] = MPHornerRatFunMap();
      }

    // Remove RowInfo of deleted equations
    for (unsigned i=0; i<n_rows; ++i)
      if (!needed[i])
        rinfo[i] = RowInfo();

    // TODO: we could also shrink c and cmap, not doing it for now
  }

  Ret AnalyticSparseSolver::fill_matrix(Context * ctxt,
                                        unsigned n_rows,
                                        const unsigned rows[],
                                        AlgInput xin[], Mod mod,
                                        AlgorithmData * datain,
                                        SparseMatrix & m) const
  {
    auto & data = *static_cast<AnalyticSparseSolverData*>(datain);

    if (mod.n() != data.this_mod_)
      data.reset_mod_(*this, mod);

    data.reset_evals_(*this);

    const SparseLinearSolver::flag_t * info = xinfo();
    const unsigned nin = nparsin[0];
    const UInt * xi = xin[0];

    UInt * xp = data.xp_.get();
    precomp_array_mul_shoup(xi, nin, mod, xp);

    UInt * ww = ctxt->ww.get();
    const HornerRatFunPtr * c = data.c.data();

    for (unsigned i=0; i<n_rows; ++i) {

      SparseMatrixRow & r = m.row(i);

      const std::size_t row_size = rinfo[rows[i]].size;
      const unsigned * cols = rinfo[rows[i]].cols.get();
      const unsigned * cols_end = cols + row_size;
      const std::size_t * idx = rinfo[rows[i]].idx.get();

      r.resize(row_size);
      unsigned oj=0;

      for (; cols<cols_end; ++cols, ++idx) {
        unsigned col = *cols;
        if (info[col] & LSVar::IS_NON_ZERO) {
          UInt & res = data.evals_[*idx];
          if (res == MISSING_SAMPLES)
            res = c[*idx].eval(nin, ww, xi, xp, mod);
          if (res == FAILED || res == 0)
            return FAILED;
          r.el(oj).col = col;
          r.el(oj).val.set(res);
          ++oj;
        }
      }
      r.el(oj).col = SparseMatrixRow::END;
      r.resize(oj);
    }

    return SUCCESS;
  }

  AlgorithmData::Ptr
  AnalyticSparseSolver::clone_data(const AlgorithmData * datain) const
  {
    typedef AnalyticSparseSolverData Data;
    const auto & data = *static_cast<const Data*>(datain);

    std::unique_ptr<Data> ptr(new Data());
    auto & newalg = *ptr;

    copy_data(datain, ptr.get());

    const SparseLinearSolver::flag_t * info = xinfo();
    std::size_t nindepeqs = n_indep_eqs();
    const unsigned * ieq = indep_eqs();
    const unsigned nin = nparsin[0];

    newalg.c.resize(cmap.size());
    newalg.reset_evals_(*this);

    for (unsigned ind=0; ind<nindepeqs; ++ind) {
      std::size_t i = ieq[ind];
      std::size_t row_size = rinfo[i].size;
      const std::size_t * idx = rinfo[i].idx.get();
      const unsigned * cols = rinfo[i].cols.get();
      const unsigned * cols_end = cols + row_size;
      for (; cols<cols_end; ++cols, ++idx)
        if ((info[*cols] & LSVar::IS_NON_ZERO)
            && newalg.evals_[*idx] == MISSING_SAMPLES) {
          horner_ratfun_clone(data.c[*idx], nin, newalg.c[*idx]);
          newalg.evals_[*idx] = 1;
        }
    }

    return std::move(ptr);
  }

} // namespace fflow
