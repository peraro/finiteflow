#include <fstream>
#include <fflow/alg_linear_solver.hh>

// In the learning phase, out of all the needed vars, we select the
// needed ones which are non-zero, and the subset of independent vars
// which appear in their solution.  Thus we define the corresponding
// coefficients.

namespace fflow {

  void DenseLinearSolver::copy_data(const AlgorithmData * data_in,
                                    AlgorithmData * oth_in) const
  {
    const auto * data = static_cast<const DenseLinearSolverData*>(data_in);
    auto * oth = static_cast<DenseLinearSolverData*>(oth_in);
    oth->mat_.reset(data->mat_.nrows(), data->mat_.ncolumns());
    oth->depv_ = data->depv_;
  }

  void DenseLinearSolver::reset(AlgorithmData * data)
  {
    adata_(data).depv_.clear();
    indepv_.clear();

    std::fill(xinfo_.get(), xinfo_.get()+nvars()+1, LSVar::IS_NON_ZERO);
    for (unsigned i=0; i<nneeded_ext_; ++i)
      xinfo_[needed_ext_[i]] |= LSVar::IS_NEEDED;

    nnindepeqs_ = adata_(data).mat_.nrows();
    for (unsigned i=0; i<nnindepeqs_; ++i)
      indepeqs_[i] = i;
  }

  Ret DenseLinearSolver::reset_needed(AlgorithmData * data,
                                      const unsigned * needed_vars,
                                      unsigned needed_size)
  {
    for (unsigned i=0; i<needed_size; ++i)
      if (needed_vars[i] >= nvars() ||
          !(xinfo_[needed_vars[i]] & flag_t(LSVar::IS_NEEDED)))
        return FAILED;

    nneeded_ext_ = needed_size;
    needed_ext_.reset(new std::size_t[needed_size]);
    for (unsigned i=0; i<needed_size; ++i)
      needed_ext_[i] = needed_vars[i];

    unsigned nv = nvars();
    for (unsigned i=0; i<nv; ++i)
      xinfo_[i] &= ~flag_t(LSVar::IS_NEEDED);
    for (unsigned i=0; i<nneeded_ext_; ++i)
      xinfo_[needed_ext_[i]] |= LSVar::IS_NEEDED;

    if (stage_ >= SECOND_) {
      MatrixView mv = mat_(data).getMatrixView(nnindepeqs_, nvars()+1);
      learn_needed_(data, mv);
      nparsout = nndeps_ * (nnindeps_ + unsigned(!homog_));
    }

    return SUCCESS;
  }

  void DenseLinearSolver::number_eqs_(AlgorithmData * data) const
  {
    for (unsigned i=0; i<nnindepeqs_; ++i)
      mat_(data)(i,mat_(data).ncolumns()-1) = i;
  }

  void DenseLinearSolver::get_dependent_variables_(AlgorithmData * data,
                                                   MatrixView & mv)
  {
    indepv_.clear();
    mv.dependent_vars(adata_(data).depv_);
    for (auto dep : adata_(data).depv_)
      xinfo_[dep] |= LSVar::IS_DEP;
    for (unsigned i=0; i<nvars(); ++i)
      if (!(xinfo_[i] & LSVar::IS_DEP))
        indepv_.push_back(i);
  }

  Ret DenseLinearSolver::check_dependent_variables_(AlgorithmData * data,
                                                    MatrixView & mv) const
  {
    std::size_t ndeps = adata_(data).depv_.size();
    mv.dependent_vars(adata_(data).depv_);
    if (ndeps != adata_(data).depv_.size()) {
      adata_(data).depv_.resize(ndeps);
      return FAILED;
    }
    for (auto dep : adata_(data).depv_)
      if (!(xinfo_[dep] & LSVar::IS_DEP))
        return FAILED;
    return SUCCESS;
  }

  // must be called after get_dependent_variables_ or
  // check_dependent_variables_
  void DenseLinearSolver::get_independent_eqs_(AlgorithmData * data,
                                               MatrixView & mv)
  {
    nnindepeqs_ = adata_(data).depv_.size();
    for (unsigned eq = 0; eq<nnindepeqs_; ++eq)
      indepeqs_[eq] = mv(eq, mat_(data).ncolumns()-1);
  }

  void DenseLinearSolver::replace_zeroes_(MatrixView & mv) const
  {
    unsigned nv = nvars();
    for (unsigned i=0; i<mv.nrows(); ++i)
      for (unsigned j=0; j<nv; ++j)
        if (!(xinfo_[j] & LSVar::IS_NON_ZERO))
          mv(i,j) = 0;
  }

  // must be called after get_dependent_variables_ or
  // check_dependent_variables_
  void DenseLinearSolver::learn_zeroes_(AlgorithmData * data, MatrixView & mv)
  {
    if (only_non_homog_)
      return learn_zeroes_onlynonhomog_(data, mv);

    unsigned eq = 0;
    zero_vars_ = 0;

    for (auto dep : adata_(data).depv_) {

      bool non_zero = false;
      for (unsigned j=dep+1; j<nvars()+1; ++j)
        if (mv(eq,j)) {
          non_zero = true;
          break;
        }

      if (!non_zero) {
        xinfo_[dep] &= ~flag_t(LSVar::IS_NON_ZERO);
        ++zero_vars_;
      }

      ++eq;
    }
  }

  void DenseLinearSolver::learn_zeroes_onlynonhomog_(AlgorithmData * data,
                                                     MatrixView & mv)
  {
    unsigned eq = 0;
    zero_vars_ = 0;

    for (auto dep : adata_(data).depv_) {

      bool non_zero = mv(eq,nvars()) != 0;

      if (!non_zero) {
        xinfo_[dep] &= ~flag_t(LSVar::IS_NON_ZERO);
        ++zero_vars_;
      }

      ++eq;
    }

    for (auto indep : indepv_) {
      xinfo_[indep] |= LSVar::IS_DEP;
      xinfo_[indep] &= ~flag_t(LSVar::IS_NON_ZERO);
    }
    zero_vars_ += indepv_.size();
    indepv_.clear();
  }

  // must be called after get_dependent_variables_ or
  // check_dependent_variables_, on a reduced systems without zero
  // variables
  void DenseLinearSolver::learn_needed_(AlgorithmData * data, MatrixView & mv)
  {
    unsigned eq = 0;
    nndeps_ = 0;
    for (auto dep : adata_(data).depv_) {
      if (xinfo_[dep] & LSVar::IS_NEEDED) {

        for (unsigned j=dep+1; j<nvars(); ++j)
          if (mv(eq,j))
            xinfo_[j] |= LSVar::IS_NEEDED;
        ++nndeps_;

      }
      ++eq;
    }

    nnindeps_ = 0;
    for (unsigned j=0; j<nvars(); ++j)
      if ((xinfo_[j] & LSVar::IS_NEEDED) && !(xinfo_[j] & LSVar::IS_DEP))
        ++nnindeps_;

    needed_dep_.reset(new std::size_t[nndeps_]);
    needed_indep_.reset(new std::size_t[nnindeps_]);
    unsigned ndc=0, nic=0;
    for (unsigned j=0; j<nvars(); ++j)
      if (xinfo_[j] & LSVar::IS_NEEDED && xinfo_[j] & LSVar::IS_NON_ZERO) {
        if (xinfo_[j] & LSVar::IS_DEP)
          needed_dep_[ndc++] = j;
        else
          needed_indep_[nic++] = j;
      }
  }

  Ret DenseLinearSolver::learn_1_(Context * ctxt,
                                  AlgInput xin[], Mod mod,
                                  AlgorithmData * data)
  {
    reset(data);

    auto mv = mat_(data).getMatrixView(nnindepeqs_, nvars()+1);
    Ret ret = fill_matrix(ctxt, nnindepeqs_, indepeqs_.get(),
                          xin, mod, data, mv);
    if (ret == FAILED)
      return FAILED;

    DynamicMatrix mat2(mv.nrows(), mv.ncolumns());
    for (unsigned i=0; i<mv.nrows(); ++i)
      for (unsigned j=0; j<mv.ncolumns(); ++j)
        mat2(i,j) = mv(i,j);

    mv.toReducedRowEcholon(mod);
    if (mv.isImpossibleSystem()) {
      stage_ = SECOND_;
      return SUCCESS;
    }

    get_dependent_variables_(data, mv);
    learn_zeroes_(data, mv);

    number_eqs_(data);
    for (unsigned i=0; i<mv.nrows(); ++i)
      for (unsigned j=0; j<mv.ncolumns(); ++j)
        mv(i,j) = mat2(i,j);

    replace_zeroes_(mv);
    mv.toReducedRowEcholon(mod);

    get_dependent_variables_(data, mv);
    get_independent_eqs_(data, mv);
    learn_needed_(data, mv);

    stage_ = SECOND_;

    return SUCCESS;
  }

  Ret DenseLinearSolver::learn_2_(Context * ctxt, AlgInput xin[], Mod mod,
                                  AlgorithmData * data)
  {
    auto mv = mat_(data).getMatrixView(nnindepeqs_, nvars()+1);
    Ret ret = fill_matrix(ctxt, nnindepeqs_, indepeqs_.get(),
                          xin, mod, data, mv);
    if (ret == FAILED)
      return FAILED;

    replace_zeroes_(mv);
    mv.toReducedRowEcholon(mod);

    if (mv.isImpossibleSystem()) {
      if (is_learning_impossible_(data))
        return SUCCESS;
      invalidate();
      return FAILED;
    }

    if (check_dependent_variables_(data, mv) != SUCCESS) {
      invalidate();
      return FAILED;
    }

    nparsout = nndeps_ * (nnindeps_ + unsigned(!homog_));

    return SUCCESS;
  }

  Ret DenseLinearSolver::learn(Context * ctxt,
                               AlgInput xin[], Mod mod, AlgorithmData * data)
  {
    if (stage_ == FIRST_)
      return learn_1_(ctxt, xin, mod, data);
    return learn_2_(ctxt, xin, mod, data);
  }


  Ret DenseLinearSolver::evaluate(Context * ctxt, AlgInput xin[], Mod mod,
                                  AlgorithmData * data,
                                  UInt xout[]) const
  {
    auto mv = mat_(data).getMatrixView(nnindepeqs_, nvars()+1);
    Ret ret = fill_matrix(ctxt, nnindepeqs_, indepeqs_.get(),
                          xin, mod, data, mv);
    if (ret == FAILED)
      return FAILED;

    replace_zeroes_(mv);
    mv.toReducedRowEcholon(mod);

    if (mv.isImpossibleSystem())
      return FAILED;

    if (check_dependent_variables_(data, mv) != SUCCESS)
      return FAILED;

    unsigned pos = 0, ndeps = adata_(data).depv_.size();
    std::size_t * depv = adata_(data).depv_.data();
    unsigned nv = nvars();
    for (unsigned i=0; i<ndeps; ++i) {

      if (xinfo_[depv[i]] & LSVar::IS_NEEDED) {
        for (unsigned j=0; j<nnindeps_; ++j)
          xout[pos++] = neg_mod(mv(i,needed_indep_[j]), mod);
        if (!homog_)
          xout[pos++] = mv(i, nv);
      }

    }

    return SUCCESS;
  }

  Ret DenseLinearSolver::only_homogeneous(bool flag)
  {
    if (!is_mutable())
      return FAILED;
    invalidate();
    homog_ = flag;
    return SUCCESS;
  }

  Ret DenseLinearSolver::only_non_homogeneous(bool flag)
  {
    if (!is_mutable())
      return FAILED;
    invalidate();
    only_non_homog_ = flag;
    return SUCCESS;
  }


  // Sparse

  void SparseLinearSolver::copy_data(const AlgorithmData * datain,
                                     AlgorithmData * dataout) const
  {
    const auto & data = *static_cast<const SparseLinearSolverData*>(datain);
    auto & oth = *static_cast<SparseLinearSolverData*>(dataout);
    oth.mat_.resize(data.mat_.nrows(), data.mat_.ncolumns());
  }

  Ret SparseLinearSolver::dump_info(const AlgorithmData * datain,
                                    const char * filename) const
  {
    const auto & data = *static_cast<const SparseLinearSolverData*>(datain);
    std::ofstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto dump_int32 = [&file](unsigned z)
      {
        file.write(reinterpret_cast<char*>(&z), sizeof(unsigned));
      };

    auto dump_flag = [&file](flag_t z)
      {
        file.write(reinterpret_cast<char*>(&z), sizeof(flag_t));
      };

    dump_int32(neqs_);
    dump_int32(nvars_);
    dump_int32(data.mat_.nrows());
    dump_int32(data.mat_.ncolumns());
    dump_int32(nndeps_);
    dump_int32(nnindeps_);
    dump_int32(nnindepeqs_);
    dump_int32(nparsin[0]);
    dump_int32(nparsout);
    dump_int32(maxcol_);
    dump_int32(zerodeps_.size());

#define dump_uptr(field, size)              \
    for (unsigned j=0; j<size; ++j)         \
      dump_int32(field[j]);

    dump_uptr(needed_dep_, nndeps_);
    dump_uptr(needed_indep_, nnindeps_);
    dump_uptr(indepeqs_, nnindepeqs_);
    dump_uptr(outeq_pos_, nndeps_-zerodeps_.size());
    dump_uptr(zerodeps_, zerodeps_.size());

    for (unsigned j=0; j<nvars()+1; ++j)
      dump_flag(xinfo_[j]);

    dump_flag(stage_);
    dump_flag(flag_);

    // sparse-output data
    if (!output_is_sparse()) {

      dump_flag(0);

    }  else {

      dump_flag(1);
      const auto & spoutd = *sparseout_data_.get();
      dump_int32(spoutd.size());
      for (const auto & eqdata : spoutd) {
        dump_int32(eqdata.size());
        dump_uptr(eqdata.data(), eqdata.size());
      }

    }

    // eq. dependencies
    dump_int32(eqdeps_.size());
    for (const auto & eqdep : eqdeps_) {
      dump_int32(eqdep.size());
      dump_uptr(eqdep, eqdep.size());
    }


#undef dump_uptr

    return SUCCESS;
  }

  Ret SparseLinearSolver::load_info(AlgorithmData * datain,
                                    const char * filename)
  {
    auto & data = *static_cast<SparseLinearSolverData*>(datain);
    std::ifstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto load_int32 = [&file]() -> unsigned
      {
        unsigned z;
        file.read(reinterpret_cast<char*>(&z), sizeof(unsigned));
        return z;
      };

    auto load_flag = [&file]() -> flag_t
      {
        flag_t z;
        file.read(reinterpret_cast<char*>(&z), sizeof(flag_t));
        return z;
      };

#define LOAD_INT32(z)   \
    do {                \
      z = load_int32(); \
      if (file.fail())  \
        return FAILED;  \
    } while (0)
#define LOAD_FLAG(z)    \
    do {                \
      z = load_flag();  \
      if (file.fail())  \
        return FAILED;  \
    } while (0)

    LOAD_INT32(neqs_);
    LOAD_INT32(nvars_);
    unsigned nr, nc;
    LOAD_INT32(nr);
    LOAD_INT32(nc);
    data.mat_.resize(nr, nc);
    //
    LOAD_INT32(nndeps_);
    needed_dep_.reset(new unsigned[nndeps_]);
    LOAD_INT32(nnindeps_);
    needed_indep_.reset(new unsigned[nnindeps_]);
    LOAD_INT32(nnindepeqs_);
    indepeqs_.reset(new unsigned[nnindepeqs_]);
    nparsin.resize(1);
    LOAD_INT32(nparsin[0]);
    LOAD_INT32(nparsout);
    LOAD_INT32(maxcol_);
    unsigned zerod_size;
    LOAD_INT32(zerod_size);
    zerodeps_.resize(zerod_size);
    outeq_pos_.reset(new unsigned[nndeps_-zerod_size]);

#define load_uptr(field, size)              \
    for (unsigned j=0; j<size; ++j)         \
      LOAD_INT32(field[j]);

    load_uptr(needed_dep_, nndeps_);
    load_uptr(needed_indep_, nnindeps_);
    load_uptr(indepeqs_, nnindepeqs_);
    load_uptr(outeq_pos_, nndeps_-zerod_size);
    load_uptr(zerodeps_, zerod_size);

    xinfo_.reset(new flag_t[nvars()+1]);
    for (unsigned j=0; j<nvars()+1; ++j)
      LOAD_FLAG(xinfo_[j]);

    LOAD_FLAG(stage_);
    LOAD_FLAG(flag_);

    // sparse-output data
    flag_t is_sparse;
    LOAD_FLAG(is_sparse);
    if (is_sparse) {

      sparseout_data_.reset(new std::vector<std::vector<unsigned>>());
      auto & spoutd = *sparseout_data_.get();
      unsigned spoutd_size;
      LOAD_INT32(spoutd_size);
      spoutd.resize(spoutd_size);
      for (unsigned j=0; j<spoutd_size; ++j) {
        unsigned eqdata_size;
        LOAD_INT32(eqdata_size);
        spoutd[j].resize(eqdata_size);
        load_uptr(spoutd[j], eqdata_size);
      }

    }

    // eq. dependencies
    unsigned eqdep_size;
    LOAD_INT32(eqdep_size);
    eqdeps_.resize(eqdep_size);
    for (unsigned j=0; j<eqdep_size; ++j) {
      unsigned thiseqdep_size;
      LOAD_INT32(thiseqdep_size);
      eqdeps_[j].resize(thiseqdep_size);
      load_uptr(eqdeps_[j], thiseqdep_size);
    }

#undef load_uptr
#undef LOAD_INT32
#undef LOAD_FLAG

    return SUCCESS;
  }


  Ret SparseLinearSolver::reset(AlgorithmData * data)
  {
    nnindepeqs_ = adata_(data).mat_.nrows();
    for (unsigned i=0; i<nnindepeqs_; ++i)
      indepeqs_[i] = i;

    zerodeps_.clear();

    return SUCCESS;
  }

  void SparseLinearSolver::set_ext_needed_(const unsigned * needed_vars,
                                           unsigned needed_size)
  {
    std::fill(xinfo_.get(), xinfo_.get()+nvars()+1, LSVar::IS_NON_ZERO);
    for (unsigned j=0; j<needed_size; ++j)
      xinfo_[needed_vars[j]] |= (LSVar::IS_NEEDED | LSVar::IS_NEEDED_EXT);
  }

  void SparseLinearSolver::get_needed_indep_()
  {
    nnindeps_ = 0;
    for (unsigned j=0; j<nvars(); ++j)
      if ((xinfo_[j] & LSVar::IS_NEEDED) && !(xinfo_[j] & LSVar::IS_DEP))
        ++nnindeps_;

    unsigned nic=0;
    needed_indep_.reset(new unsigned[nnindeps_]);
    for (unsigned j=0; j<nvars(); ++j)
      if (xinfo_[j] & LSVar::IS_NEEDED && !(xinfo_[j] & LSVar::IS_DEP))
        needed_indep_[nic++] = j;
  }

  void SparseLinearSolver::relearn_needed_(AlgorithmData * data)
  {
    const unsigned old_nz = zerodeps_.size();

    relearn_zero_needed_();
    std::copy(zerodeps_.begin(), zerodeps_.end(), needed_dep_.get());

    const unsigned nz = zerodeps_.size();
    const auto & mat = mat_(data);
    unsigned ieq=0;
    for (unsigned j=old_nz; j<nndeps_; ++j) {
      const unsigned eqpos = outeq_pos_[j-old_nz];
      const auto & row = mat.row(eqpos);
      const unsigned depv = row.first_nonzero_column();
      if (xinfo_[depv] & LSVar::IS_NEEDED) {
        const unsigned rsize = row.size();
        for (unsigned k=0; k<rsize; ++k)
          xinfo_[row.el(k).col] |= LSVar::IS_NEEDED;
        outeq_pos_[ieq] = eqpos;
        needed_dep_[ieq+nz] = depv;
        ++ieq;
      }
    }
    nndeps_ = ieq;

    if (!output_is_sparse())
      get_needed_indep_();
  }

  void SparseLinearSolver::optimize_nonneeded_indeps_()
  {
    // Optimization: remove non-needed independent variables
    for (unsigned j=0; j<nvars_; ++j)
      if (!(xinfo_[j] & LSVar::IS_DEP) && !(xinfo_[j] & LSVar::IS_NEEDED))
        xinfo_[j] &= ~flag_t(LSVar::IS_NON_ZERO);
  }

  Ret SparseLinearSolver::reset_needed(AlgorithmData * data,
                                       const unsigned * needed_vars,
                                       unsigned needed_size)
  {
    const unsigned nv = nvars();

    for (unsigned i=0; i<needed_size; ++i)
      if (needed_vars[i] >= nv ||
          !(xinfo_[needed_vars[i]] & flag_t(LSVar::IS_NEEDED)))
        return FAILED;

    eqdeps_.clear();

    for (unsigned i=0; i<nv; ++i)
      xinfo_[i] &= ~flag_t(LSVar::IS_NEEDED | LSVar::IS_NEEDED_EXT);
    for (unsigned i=0; i<needed_size; ++i)
      xinfo_[needed_vars[i]] |= (LSVar::IS_NEEDED | LSVar::IS_NEEDED_EXT);

    if (stage_ >= SECOND_) {
      relearn_needed_(data);
      optimize_nonneeded_indeps_();
      if (output_is_sparse())
        set_sparseout_data_(data);
      else
        nparsout = nndeps_ * (nnindeps_ + unsigned(!(flag_ & HOMOG_)));
    }

    return SUCCESS;
  }

  void SparseLinearSolver::number_eqs_(AlgorithmData * data)
  {
    for (unsigned i=0; i<nnindepeqs_; ++i)
      adata_(data).mat_.row(i).id() = i;
  }

  Ret SparseLinearSolver::check_dependent_variables_(AlgorithmData * data) const
  {
    const auto & mat = adata_(data).mat_;
    const unsigned neqs = mat.nrows();

    unsigned check_nneeded_deps = zerodeps_.size();

    for (unsigned j=0; j<neqs; ++j) {

      auto & row = mat.row(j);
      const unsigned col = row.first_nonzero_column();

      if (col == SparseMatrixRow::END)
        break;

      if (!(xinfo_[col] & LSVar::IS_DEP) || (xinfo_[col] & LSVar::IS_SWEEPED))
        return FAILED;

      if (xinfo_[col] & LSVar::IS_NEEDED) {
        ++check_nneeded_deps;

        if (row.size() != 1) {
          if (!(xinfo_[col] & LSVar::IS_NON_ZERO))
            return FAILED;
        } else {
          if ((xinfo_[col] & LSVar::IS_NON_ZERO))
            return FAILED;
        }

      }
    }

    if (check_nneeded_deps != nndeps_)
      return FAILED;

    return SUCCESS;
  }

  void SparseLinearSolver::relearn_zero_needed_()
  {
    unsigned nz = std::remove_if(zerodeps_.begin(),
                                 zerodeps_.end(),
                                 [&](unsigned col)
                                 {
                                   return !(xinfo_[col] & LSVar::IS_NEEDED);
                                 }) - zerodeps_.begin();
    zerodeps_.resize(nz);
  }

  bool SparseLinearSolver::check_zeroes_(AlgorithmData * data)
  {
    if (flag_ & ONLY_NON_HOMOG_)
      return check_zeroes_onlynonhomog_(data);

    const auto & mat = adata_(data).mat_;
    const unsigned neqs = mat.nrows();
    bool any_zero = false;

    if (!zerodeps_.size()) {

      for (unsigned j=0; j<neqs; ++j) {

        auto & row = mat.row(j);
        const unsigned col = row.first_nonzero_column();

        if (col == SparseMatrixRow::END)
          break;

        if (row.size() == 1) {
          any_zero = true;
          xinfo_[col] |= LSVar::IS_DEP;
          xinfo_[col] &= ~flag_t(LSVar::IS_NON_ZERO);
          if (xinfo_[col] & LSVar::IS_NEEDED)
            zerodeps_.push_back(col);
        }

      }

    }

    return any_zero;
  }

  bool SparseLinearSolver::check_zeroes_onlynonhomog_(AlgorithmData * data)
  {
    const auto & mat = adata_(data).mat_;
    const unsigned neqs = mat.nrows();
    bool any_zero = false;

    if (!zerodeps_.size()) {

      for (unsigned j=0; j<neqs; ++j) {

        auto & row = mat.row(j);
        const unsigned col = row.first_nonzero_column();

        if (col == SparseMatrixRow::END)
          break;

        xinfo_[col] |= LSVar::IS_DEP;

        bool non_zero = false;
        const unsigned rsize = row.size();
        if (row.el(rsize-1).col == nvars())
          non_zero = true;

        if (!non_zero) {
          any_zero = true;
          xinfo_[col] &= ~flag_t(LSVar::IS_NON_ZERO);
          if (xinfo_[col] & LSVar::IS_NEEDED)
            zerodeps_.push_back(col);
        }

      }

      for (unsigned j=0; j<nvars(); ++j)
        if (!(xinfo_[j] & LSVar::IS_DEP)) {
          any_zero = true;
          xinfo_[j] |= LSVar::IS_DEP;
          xinfo_[j] &= ~flag_t(LSVar::IS_NON_ZERO);
          if (xinfo_[j] & LSVar::IS_NEEDED)
            zerodeps_.push_back(j);
        }
    }

    return any_zero;
  }

  Ret SparseLinearSolver::learn_1_(Context * ctxt, AlgInput xin[], Mod mod,
                                   AlgorithmData * data)
  {
    reset(data);

    Ret ret = fill_matrix(ctxt, nnindepeqs_, indepeqs_.get(),
                          xin, mod, data, mat_(data));
    if (ret == FAILED)
      return FAILED;

    number_eqs_(data);
    mat_(data).sortRows();

    mat_(data).toReducedRowEcholon(mod, maxcol_, !(flag_ & NO_BACKSUBST_),
                                   !(flag_ & KEEP_ALL_OUTS_));

    if (mat_(data).isImpossibleSystem()) {
      flag_ |= IMPOSSIBLE_;
      stage_ = SECOND_;
      return SUCCESS;
    }

    // Deal with zeroes first
    bool any_zero = check_zeroes_(data);
    if (any_zero) {
      Ret ret = fill_matrix(ctxt, nnindepeqs_, indepeqs_.get(),
                            xin, mod, data, mat_(data));
      if (ret == FAILED)
        return FAILED;
      number_eqs_(data);
      mat_(data).sortRows();
      mat_(data).toReducedRowEcholon(mod, maxcol_, !(flag_ & NO_BACKSUBST_),
                                     !(flag_ & KEEP_ALL_OUTS_));;
    }

    // Filter equations
    nnindepeqs_ = adata_(data).mat_.removeZeroDeps();

    // get dependent variables and needed
    {
      const auto & mat = adata_(data).mat_;
      const unsigned neqs = mat.nrows();
      nndeps_ = zerodeps_.size();

      for (unsigned j=0; j<neqs; ++j) {

        auto & row = mat.row(j);
        const unsigned col = row.first_nonzero_column();

        if (col == SparseMatrixRow::END)
          break;

        xinfo_[col] |= LSVar::IS_DEP;

        if (xinfo_[col] & LSVar::IS_NEEDED) {

          const unsigned rsize = row.size();
          for (unsigned k=1; k<rsize; ++k)
            xinfo_[row.el(k).col] |= LSVar::IS_NEEDED;

          ++nndeps_;
        }

      }

      if (!output_is_sparse())
        get_needed_indep_();

      outeq_pos_.reset(new unsigned[nndeps_ - zerodeps_.size()]);
      needed_dep_.reset(new unsigned[nndeps_]);
      std::copy(zerodeps_.begin(), zerodeps_.end(), needed_dep_.get());
      unsigned ceqpos = 0;
      unsigned ndc = zerodeps_.size();
      for (unsigned j=0; j<neqs; ++j) {

        auto & row = mat.row(j);
        const unsigned col = row.first_nonzero_column();

        if (col == SparseMatrixRow::END)
          break;

        if (xinfo_[col] & LSVar::IS_NEEDED) {
          needed_dep_[ndc++] = row.first_nonzero_column();
          outeq_pos_[ceqpos++] = j;
        }
      }
    }

    // get independent eqs
    {
      const auto & mat = adata_(data).mat_;
      for (unsigned eq = 0; eq<nnindepeqs_; ++eq)
        indepeqs_[eq] = mat.row(eq).id();
    }

    optimize_nonneeded_indeps_();

    stage_ = SECOND_;

    return SUCCESS;
  }

  Ret SparseLinearSolver::learn_2_(Context * ctxt, AlgInput xin[], Mod mod,
                                   AlgorithmData * data)
  {
    Ret ret = fill_matrix(ctxt, nnindepeqs_, indepeqs_.get(),
                          xin, mod, data, mat_(data));
    if (ret == FAILED)
      return FAILED;

    eqdeps_.resize(nnindepeqs_);
    mat_(data).toReducedRowEcholon(mod, maxcol_, !(flag_ & NO_BACKSUBST_),
                                   !(flag_ & KEEP_ALL_OUTS_),
                                   xinfo_.get(), eqdeps_.data());

    if (mat_(data).isImpossibleSystem()) {
      if (is_learning_impossible_(data))
        return SUCCESS;
      invalidate();
      return FAILED;
    }

    if (check_dependent_variables_(data) != SUCCESS) {
      invalidate();
      return FAILED;
    }

    if (output_is_sparse())
      set_sparseout_data_(data);
    else
      nparsout = nndeps_ * (nnindeps_ + unsigned(!(flag_ & HOMOG_)));

    return SUCCESS;
  }

  Ret SparseLinearSolver::learn(Context * ctxt, AlgInput xin[], Mod mod,
                                AlgorithmData * data)
  {
    if (!is_mutable())
      return MUTABILITY_ERROR;
    if (stage_ == FIRST_)
      return learn_1_(ctxt, xin, mod, data);
    return learn_2_(ctxt, xin, mod, data);
  }


  Ret SparseLinearSolver::evaluate(Context * ctxt, AlgInput xin[], Mod mod,
                                   AlgorithmData * data,
                                   UInt xout[]) const
  {
    Ret ret = fill_matrix(ctxt, nnindepeqs_, indepeqs_.get(),
                          xin, mod, data, mat_(data));
    if (ret == FAILED)
      return FAILED;

    SparseMatrix & mat = mat_(data);

    mat.toReducedRowEcholon(mod, maxcol_, !(flag_ & NO_BACKSUBST_),
                            !(flag_ & KEEP_ALL_OUTS_), xinfo_.get());

    if (mat.isImpossibleSystem())
      return FAILED;

    if (check_dependent_variables_(data) != SUCCESS)
      return FAILED;

    unsigned pos = 0;
    unsigned ndeps = nndeps_ - zerodeps_.size();

    if (!output_is_sparse()) {

      if (zerodeps_.size()) {
        const unsigned n0 = zerodeps_.size() * (nnindeps_ + unsigned(!(flag_ & HOMOG_)));
        std::fill(xout, xout + n0, 0);
        pos += n0;
      }

      for (unsigned i=0; i<ndeps; ++i) {

        const unsigned eq = outeq_pos_[i];
        for (unsigned j=0; j<nnindeps_; ++j)
          xout[pos++] = neg_mod(mat(eq,needed_indep_[j]), mod);
        if (!(flag_ & HOMOG_))
          xout[pos++] = mat(eq, nvars());

      }

    } else {

      bool homog = (flag_ & HOMOG_);
      for (unsigned i=0; i<ndeps; ++i) {

        const auto & elems = (*sparseout_data_)[i+zerodeps_.size()];
        const unsigned elsize = elems.size();
        const auto & row = mat.row(outeq_pos_[i]);
        if (row.size() < elsize + 1)
          return FAILED;

        for (unsigned j=0; j<elsize; ++j) {
          if (homog || elems[j] != nvars())
            xout[pos++] = neg_mod(row.el(j+1).val.get(), mod);
          else
            xout[pos++] = row.el(j+1).val.get();
        }

      }

    }

    return SUCCESS;
  }


  void SparseLinearSolver::mark_eq_(const SparseMatrix::EqDeps * eqdeps,
                                    unsigned eq,
                                    bool * marked)
  {
    if (marked[eq])
      return;

    marked[eq] = true;
    for (unsigned deq : eqdeps[eq])
      mark_eq_(eqdeps, deq, marked);
  }

  void SparseLinearSolver::mark_and_sweep_eqs(AlgorithmData * data)
  {
    if (eqdeps_.empty())
      return;

    std::unique_ptr<bool[]> marked(new bool[nnindepeqs_]());

    auto & mat = adata_(data).mat_;
    for (unsigned j=0; j<nnindepeqs_; ++j) {
      auto & row = mat.row(j);
      const unsigned col = row.first_nonzero_column();
      if (col == SparseMatrixRow::END)
        break;
      if (xinfo_[col] & LSVar::IS_NEEDED)
        mark_eq_(eqdeps_.data(), j, marked.get());
    }

    unsigned new_nnindepeqs=0;
    unsigned neqs_old = mat.nrows();
    unsigned outeq=0;
    for (unsigned i=0; i<neqs_old; ++i) {
      const auto & row = mat.row(i);
      unsigned depv = row.first_nonzero_column();
      if (marked[i]) {
        indepeqs_[new_nnindepeqs] = row.id();
        if (xinfo_[depv] & LSVar::IS_NEEDED)
          outeq_pos_[outeq++] = new_nnindepeqs;
        ++new_nnindepeqs;
      } else {
        xinfo_[depv] |= LSVar::IS_SWEEPED;
      }
    }
    nnindepeqs_ = new_nnindepeqs;
    mat.restrict_rows(nnindepeqs_);

    eqdeps_.clear();
    eqdeps_.shrink_to_fit();

    flag_ |= MARKED_AND_SWEEPED_;
  }

  Ret SparseLinearSolver::only_homogeneous(bool flag)
  {
    if (!is_mutable())
      return FAILED;
    invalidate();
    if (flag)
      flag_ |= HOMOG_;
    else
      flag_ &= ~flag_t(HOMOG_);
    return SUCCESS;
  }

  Ret SparseLinearSolver::only_non_homogeneous(bool flag)
  {
    if (!is_mutable() || stage_ != FIRST_)
      return FAILED;
    invalidate();
    if (flag)
      flag_ |= ONLY_NON_HOMOG_;
    else
      flag_ &= ~flag_t(ONLY_NON_HOMOG_);
    return SUCCESS;
  }

  Ret SparseLinearSolver::sparse_output(bool flag)
  {
    if (!is_mutable())
      return FAILED;
    invalidate();

    if (!flag) {
      sparseout_data_.reset(nullptr);
      return SUCCESS;
    }

    if (!output_is_sparse())
      sparseout_data_.reset(new std::vector<std::vector<unsigned>>());

    return SUCCESS;
  }

  Ret SparseLinearSolver::sparse_output_with_maxcol(unsigned maxcol,
                                                    bool backsubst,
                                                    bool keep_all_outs)
  {
    Ret ret = sparse_output();
    if (ret != SUCCESS)
      return ret;
    maxcol_ = maxcol;
    if (!backsubst)
      flag_ |= NO_BACKSUBST_;
    else
      flag_ &= ~flag_t(NO_BACKSUBST_);
    if (!keep_all_outs)
      flag_ |= KEEP_ALL_OUTS_;
    else
      flag_ &= ~flag_t(KEEP_ALL_OUTS_);
    return ret;
  }

  void SparseLinearSolver::set_sparseout_data_(const AlgorithmData * data)
  {
    if (!output_is_sparse())
      sparseout_data_.reset(new std::vector<std::vector<unsigned>>());

    sparseout_data_->resize(nndeps_);

    const SparseMatrix & mat = mat_(data);
    const unsigned zero_size = zerodeps_.size();

    for (unsigned i=0; i<zero_size; ++i)
      (*sparseout_data_)[i].resize(0);

    unsigned nout = 0;
    const bool homog = (flag_ & HOMOG_);
    for (unsigned i=zero_size, eqidx=0; i<nndeps_; ++i,++eqidx) {
      auto & elems = (*sparseout_data_)[i];
      const SparseMatrixRow & row = mat.row(outeq_pos_[eqidx]);
      const unsigned elsize = row.size()-1;
      elems.resize(elsize);
      for (unsigned j=0; j<elsize; ++j)
        elems[j] = row.el(j+1).col;
      if (homog && elsize>=1 && elems[elsize-1] == nvars())
        elems.pop_back();
      nout += elems.size();
    }

    nparsout = nout;
  }

} // namespace fflow
