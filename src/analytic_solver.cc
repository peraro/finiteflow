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

  static UInt bad_sparse_ccs_(UInt res)
  {
    if (res == 0)
      logerr("Sparse linear-system matrix of "
             "coefficients contains zeroes");
    else
      logerr("Evaluation of linear system-matrix "
             "of coefficients failed");
    return FAILED;
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
          if (res == MISSING_SAMPLES) {
            res = c[*idx].eval(nin, ww, xi, xp, mod);
            if (FF_ERRCOND(res == FAILED || res == 0))
              return bad_sparse_ccs_(res);
          }
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



  void
  AnalyticSparseSolverExData::reset_mod_(const AnalyticSparseSolverEx & ls,
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
      const AnalyticSparseSolverEx::Weights * ws = ls.rinfo[i].w.get();
      const unsigned * cols = ls.rinfo[i].cols.get();
      const unsigned * cols_end = cols + row_size;
      for (; cols<cols_end; ++ws, ++cols)
        if ((info[*cols] & LSVar::IS_NON_ZERO)) {
          const AnalyticSparseSolverEx::Weight * w = ws->w.get();
          for (unsigned j=0; j<ws->size; ++j, ++w)
            if (evals_[w->idx] == MISSING_SAMPLES) {
              horner_ratfunmap_mprational(ls.cmap[w->idx], mod, c[w->idx]);
              evals_[w->idx] = 1;
            }
        }
    }
  }

  void AnalyticSparseSolverEx::delete_unneeded_eqs(AlgorithmData * datain)
  {
    auto & data = *static_cast<AnalyticSparseSolverExData*>(datain);

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
      const Weights * weights = rinfo[i].w.get();
      const unsigned * cols = rinfo[i].cols.get();
      const unsigned * cols_end = cols + row_size;
      for (; cols<cols_end; ++weights, ++cols) {
        const Weight * w = weights->w.get();
        const Weight * w_end = w + weights->size;
        for (; w < w_end; ++w) {
          if (info[*cols] & LSVar::IS_NON_ZERO)
            data.evals_[w->idx] = 1;
        }
      }
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

  Ret AnalyticSparseSolverEx::learn(Context * ctxt,
                                    AlgInput xin[], Mod mod,
                                    AlgorithmData * data)
  {
    if (!learned_ccs_) {
      reset(data);
      if (learn_zero_ccs_(ctxt, xin, mod, data) != SUCCESS)
        return FAILED;
    }
    return SparseLinearSolver::learn(ctxt, xin, mod, data);
  }

  Ret AnalyticSparseSolverEx::learn_zero_ccs_(Context * ctxt,
                                              AlgInput xin[], Mod mod,
                                              AlgorithmData * datain)
  {
    auto & data = *static_cast<AnalyticSparseSolverExData*>(datain);

    if (mod.n() != data.this_mod_)
      data.reset_mod_(*this, mod);

    data.reset_evals_(*this);

    //const SparseLinearSolver::flag_t * info = xinfo();
    const unsigned nin = nparsin[0];
    const UInt * xi = xin[0];

    UInt * xp = data.xp_.get();
    precomp_array_mul_shoup(xi, nin, mod, xp);

    UInt * ww = ctxt->ww.get();
    const HornerRatFunPtr * c = data.c.data();

    const unsigned n_rows = neqs();
    std::vector<UInt> row_res;

    for (unsigned i=0; i<n_rows; ++i) {

      const std::size_t row_size = rinfo[i].size;
      const unsigned n_cols = row_size;
      const Weights * weights = rinfo[i].w.get();

      row_res.resize(row_size);
      unsigned non_zero_col = 0;

      for (unsigned oj=0; oj<n_cols; ++oj, ++weights) {
        { // NOTE: zero vars have not been checked yet
          const Weight * wgt = weights->w.get();
          const Weight * wgt_end = wgt + weights->size;
          UInt tot_res = 0;
          for (; wgt < wgt_end; ++wgt) {
            UInt & res = data.evals_[wgt->idx];
            if (res == MISSING_SAMPLES) {
              res = c[wgt->idx].eval(nin, ww, xi, xp, mod);
              if (FF_ERRCOND(res == FAILED || res == 0))
                return bad_sparse_ccs_(res);
            }
            addmul_mod(tot_res, res, xin[wgt->node][wgt->el], mod);
          }
          if (tot_res != 0)
            ++non_zero_col;
          row_res[oj] = tot_res;
        }
      }

      if (non_zero_col != row_size) {

        RowInfo & rinf = rinfo[i];
        std::unique_ptr<Weights[]> new_w(new Weights[non_zero_col]);
        std::unique_ptr<unsigned[]> new_cols(new unsigned[non_zero_col]);
        unsigned oj = 0;
        for (unsigned j=0; j<row_size; ++j) {
          if (row_res[j] != 0) {
            new_cols[oj] = rinf.cols[j];
            new_w[oj].from(std::move(rinf.w[j]));
            ++oj;
          }
        }
        rinf.cols = std::move(new_cols);
        rinf.w = std::move(new_w);
        rinf.size = non_zero_col;

      }

    }

    learned_ccs_ = true;

    return SUCCESS;
  }

  Ret AnalyticSparseSolverEx::fill_matrix(Context * ctxt,
                                          unsigned n_rows,
                                          const unsigned rows[],
                                          AlgInput xin[], Mod mod,
                                          AlgorithmData * datain,
                                          SparseMatrix & m) const
  {
    auto & data = *static_cast<AnalyticSparseSolverExData*>(datain);

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
      const Weights * weights = rinfo[rows[i]].w.get();

      r.resize(row_size);
      unsigned oj=0;

      for (; cols<cols_end; ++cols, ++weights) {
        unsigned col = *cols;
        if (info[col] & LSVar::IS_NON_ZERO) {
          const Weight * wgt = weights->w.get();
          const Weight * wgt_end = wgt + weights->size;
          UInt tot_res = 0;
          for (; wgt < wgt_end; ++wgt) {
            UInt & res = data.evals_[wgt->idx];
            if (res == MISSING_SAMPLES) {
              res = c[wgt->idx].eval(nin, ww, xi, xp, mod);
              if (FF_ERRCOND(res == FAILED || res == 0))
                return bad_sparse_ccs_(res);
            }
            addmul_mod(tot_res, res, xin[wgt->node][wgt->el], mod);
          }
          if (tot_res == 0)
            return bad_sparse_ccs_(tot_res);
          r.el(oj).col = col;
          r.el(oj).val.set(tot_res);
          ++oj;
        }
      }
      r.el(oj).col = SparseMatrixRow::END;
      r.resize(oj);
    }

    return SUCCESS;
  }

  AlgorithmData::Ptr
  AnalyticSparseSolverEx::clone_data(const AlgorithmData * datain) const
  {
    typedef AnalyticSparseSolverExData Data;
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
      const Weights * weights = rinfo[i].w.get();
      const unsigned * cols = rinfo[i].cols.get();
      const unsigned * cols_end = cols + row_size;
      for (; cols<cols_end; ++cols, ++weights) {
        const Weight * w = weights->w.get();
        const Weight * w_end = w + weights->size;
        for (; w < w_end; ++w) {
          std::size_t idx = w->idx;
          if ((info[*cols] & LSVar::IS_NON_ZERO)
              && newalg.evals_[idx] == MISSING_SAMPLES) {
            horner_ratfun_clone(data.c[idx], nin, newalg.c[idx]);
            newalg.evals_[idx] = 1;
          }
        }
      }
    }

    return std::move(ptr);
  }

} // namespace fflow
