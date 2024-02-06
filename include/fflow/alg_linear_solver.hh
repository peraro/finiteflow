#ifndef FFLOW_ALG_LINEAR_SOLVER_HH
#define FFLOW_ALG_LINEAR_SOLVER_HH

#include <memory>
#include <fflow/algorithm.hh>
#include <fflow/matrix.hh>

namespace fflow {

  class DenseLinearSolver;
  class SparseLinearSolver;


  // Dense solver

  class DenseLinearSolverData : public AlgorithmData {
  private:
    friend class DenseLinearSolver;
    std::vector<std::size_t> depv_;
    DynamicMatrix mat_;
  };

  class DenseLinearSolver : public Algorithm {
  public:

    typedef LSVar::flag_t flag_t;

    DenseLinearSolver() = default;

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    virtual UInt min_learn_times() override
    {
      return 3;
    }

    virtual Ret learn(Context * ctxt, AlgInput xin[], Mod mod,
                      AlgorithmData * data) override;

    void copy_data(const AlgorithmData * data, AlgorithmData * dataout) const;

    virtual Ret fill_matrix(Context * ctxt,
                            std::size_t n_rows, const std::size_t rows[],
                            AlgInput xi[], Mod mod,
                            AlgorithmData * data,
                            MatrixView & m) const = 0;

    void reset(AlgorithmData * data);

    void init(std::size_t neqs, std::size_t nvars,
              unsigned * needed_vars, unsigned needed_size,
              DenseLinearSolverData & data)
    {
      data.mat_.reset(neqs, nvars+2);
      xinfo_.reset(new flag_t[nvars+1]());
      indepeqs_.reset(new std::size_t[neqs]);
      nneeded_ext_ = needed_size;
      needed_ext_.reset(new std::size_t[needed_size]);
      for (unsigned i=0; i<needed_size; ++i)
        needed_ext_[i] = needed_vars[i];
      neqs_ = neqs;
      nvars_ = nvars;
      stage_ = FIRST_;
    }

    Ret reset_needed(AlgorithmData * data,
                     const unsigned * needed_vars, unsigned needed_size);

    void invalidate()
    {
      stage_ = FIRST_;
    }

    std::size_t neqs() const
    {
      return neqs_;
    }

    std::size_t nvars() const
    {
      return nvars_;
    }

    const flag_t * xinfo() const
    {
      return xinfo_.get();
    }

    std::size_t n_needed_depvars() const
    {
      return nndeps_;
    }

    std::size_t n_needed_indepvars() const
    {
      return nnindeps_;
    }

    const std::size_t * needed_depvars() const
    {
      return needed_dep_.get();
    }

    const std::size_t * needed_indepvars() const
    {
      return needed_indep_.get();
    }

    std::size_t n_indep_eqs() const
    {
      return nnindepeqs_;
    }

    const std::size_t * indep_eqs() const
    {
      return indepeqs_.get();
    }

    Ret only_homogeneous(bool flag = true);

  private:

    static DenseLinearSolverData & adata_(AlgorithmData * data)
    {
      return *static_cast<DenseLinearSolverData *>(data);
    }

    static const DenseLinearSolverData & adata_(const AlgorithmData * data)
    {
      return *static_cast<const DenseLinearSolverData *>(data);
    }

    static DynamicMatrix & mat_(AlgorithmData * data)
    {
      return adata_(data).mat_;
    }

    static const DynamicMatrix & mat_(const AlgorithmData * data)
    {
      return adata_(data).mat_;
    }

    void get_dependent_variables_(AlgorithmData * data, MatrixView & mv);
    Ret check_dependent_variables_(AlgorithmData * data,
                                   MatrixView & mv) const;

    bool is_learning_impossible_(AlgorithmData * data) const
    {
      return adata_(data).depv_.empty() && indepv_.empty();
    }

    void learn_zeroes_(AlgorithmData * data, MatrixView & mv);
    void learn_needed_(AlgorithmData * data, MatrixView & mv);

    void replace_zeroes_(MatrixView & mv) const;

    void get_independent_eqs_(AlgorithmData * data, MatrixView & mv);

    void number_eqs_(AlgorithmData * data) const;

    Ret learn_1_(Context * ctxt, AlgInput xin[], Mod mod, AlgorithmData * data);
    Ret learn_2_(Context * ctxt, AlgInput xin[], Mod mod, AlgorithmData * data);

  private:

    // Learning stages:
    //
    // FIRST_ = spurious eqs, zero vars, needed vars
    //
    // SECOND_ = needed eqs after replacing zeroes
    //
    // LEARNED_ = previous stages already completed
    enum {
      FIRST_ = 0, SECOND_ = 1, LEARNED_ = 2
    };

  private:
    //std::vector<std::size_t> depv_;
    std::vector<std::size_t> indepv_;
    std::unique_ptr<std::size_t[]> needed_ext_;
    std::unique_ptr<std::size_t[]> needed_dep_;
    std::unique_ptr<std::size_t[]> needed_indep_;
    std::unique_ptr<std::size_t[]> indepeqs_;
    std::unique_ptr<flag_t[]> xinfo_;
    std::size_t neqs_, nvars_;
    std::size_t nneeded_ext_;
    std::size_t zero_vars_ = 0, nndeps_ = 0, nnindeps_ = 0, nnindepeqs_ = 0;
    flag_t stage_ = FIRST_;
    bool homog_ = false;
  };


  // Sparse solver

  class SparseLinearSolverData : public AlgorithmData {
  private:
    friend class SparseLinearSolver;
    SparseMatrix mat_;
    std::vector<std::size_t> depv_;
  };

  class SparseLinearSolver : public Algorithm {
  public:

    typedef LSVar::flag_t flag_t;

    SparseLinearSolver() = default;

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    virtual UInt min_learn_times() override
    {
      return 3;
    }

    virtual Ret learn(Context * ctxt,
                      AlgInput xin[], Mod mod, AlgorithmData * data) override;

    virtual Ret fill_matrix(Context * ctxt,
                            std::size_t n_rows, const std::size_t rows[],
                            AlgInput xi[], Mod mod,
                            AlgorithmData * data,
                            SparseMatrix & m) const = 0;

    Ret dump_info(const AlgorithmData * data, const char * filename) const;
    Ret load_info(AlgorithmData * data, const char * filename);

    void copy_data(const AlgorithmData * data, AlgorithmData * dataout) const;

    Ret reset(AlgorithmData * data);

    void init(std::size_t neqs, std::size_t nvars,
              const unsigned * needed_vars, unsigned needed_size,
              SparseLinearSolverData & data)
    {
      data.mat_.resize(neqs, nvars+1);
      xinfo_.reset(new flag_t[nvars+1]());
      indepeqs_.reset(new std::size_t[neqs]);
      nneeded_ext_ = needed_size;
      needed_ext_.reset(new std::size_t[needed_size]);
      for (unsigned i=0; i<needed_size; ++i)
        needed_ext_[i] = needed_vars[i];
      neqs_ = neqs;
      nvars_ = nvars;
      stage_ = FIRST_;
    }

    Ret reset_needed(AlgorithmData * data,
                     const unsigned * needed_vars, unsigned needed_size);

    void invalidate()
    {
      stage_ = FIRST_;
    }

    std::size_t neqs() const
    {
      return neqs_;
    }

    std::size_t nvars() const
    {
      return nvars_;
    }

    const flag_t * xinfo() const
    {
      return xinfo_.get();
    }

    std::size_t n_needed_depvars() const
    {
      return nndeps_;
    }

    std::size_t n_needed_indepvars() const
    {
      return nnindeps_;
    }

    const std::size_t * needed_depvars() const
    {
      return needed_dep_.get();
    }

    const std::size_t * needed_indepvars() const
    {
      return needed_indep_.get();
    }

    std::size_t n_indep_eqs() const
    {
      return nnindepeqs_;
    }

    const std::size_t * indep_eqs() const
    {
      return indepeqs_.get();
    }

    void mark_and_sweep_eqs(AlgorithmData * data);

    Ret only_homogeneous(bool flag = true);

    bool marked_and_sweeped() const
    {
      return marked_and_sweeped_;
    }

    bool output_is_sparse() const
    {
      return sparseout_data_.get() != nullptr;
    }

    Ret sparse_output(bool flag = true);
    Ret sparse_output_with_maxrow(unsigned maxrow, bool back_subst = true);

    // the returned pointer may be null
    const std::vector<std::vector<unsigned>> * sparse_output_data() const
    {
      return sparseout_data_.get();
    }

  private:

    static SparseLinearSolverData & adata_(AlgorithmData * data)
    {
      return *static_cast<SparseLinearSolverData *>(data);
    }

    static const SparseLinearSolverData & adata_(const AlgorithmData * data)
    {
      return *static_cast<const SparseLinearSolverData *>(data);
    }

    static SparseMatrix & mat_(AlgorithmData * data)
    {
      return adata_(data).mat_;
    }

    static const SparseMatrix & mat_(const AlgorithmData * data)
    {
      return adata_(data).mat_;
    }

    void get_dependent_variables_(AlgorithmData * data);
    Ret check_dependent_variables_(AlgorithmData * data) const;

    bool is_learning_impossible_(const AlgorithmData * data) const
    {
      return adata_(data).depv_.empty() && indepv_.empty();
    }

    bool has_max_row_() const
    {
      return output_is_sparse() && (maxrow_ < SparseMatrixRow::END);
    }

    void learn_needed_(AlgorithmData * data);

    void get_independent_eqs_(AlgorithmData * data);

    void number_eqs_(AlgorithmData * data);

    void set_sparseout_data_(const AlgorithmData * data);

    Ret learn_1_(Context * ctxt, AlgInput xin[], Mod mod, AlgorithmData * data);
    Ret learn_2_(Context * ctxt, AlgInput xin[], Mod mod, AlgorithmData * data);

    static void mark_eq_(const SparseMatrix::EqDeps * eqdeps, unsigned eq,
                         bool * marked);

  private:

    // Learning stages:
    //
    // FIRST_ = spurious eqs, zero vars, needed vars
    //
    // SECOND_ = needed eqs after replacing zeroes
    //
    // LEARNED_ = previous stages already completed
    enum {
      FIRST_ = 0, SECOND_ = 1, LEARNED_ = 2
    };

  private:
    // std::vector<std::size_t> depv_; // <-- Moved into AlgorithmData
    std::vector<std::size_t> indepv_;
    std::unique_ptr<std::size_t[]> needed_ext_;
    std::unique_ptr<std::size_t[]> needed_dep_;
    std::unique_ptr<std::size_t[]> needed_indep_;
    std::unique_ptr<std::size_t[]> indepeqs_;
    std::unique_ptr<flag_t[]> xinfo_;
    std::size_t neqs_, nvars_;
    std::size_t nneeded_ext_;
    std::size_t nndeps_ = 0, nnindeps_ = 0, nnindepeqs_ = 0;
    std::vector<SparseMatrix::EqDeps> eqdeps_;
    std::unique_ptr<std::vector<std::vector<unsigned>>> sparseout_data_;
    UInt maxrow_ =  SparseMatrixRow::END;
    flag_t stage_ = FIRST_;
    bool backsubst_ = true, homog_ = false, marked_and_sweeped_ = false;
  };

} // namespace ampf

#endif // FFLOW_ALG_LINEAR_SOLVER_HH
