#ifndef FFLOW_ALG_LINEAR_FIT_HH
#define FFLOW_ALG_LINEAR_FIT_HH

#include <fflow/alg_linear_solver.hh>

namespace fflow {


  class LinearFit;


  class LinearFitData : public DenseLinearSolverData {
  private:
    friend class LinearFit;
    std::unique_ptr<UInt[]> tau_samples_;
    std::unique_ptr<UInt[]> tau_;
  };


  class LinearFit : public Algorithm {
  public:

    typedef DenseLinearSolver::flag_t flag_t;

    LinearFit() = default;

    virtual UInt min_learn_times() override
    {
      return ls_.min_learn_times();
    }

    virtual Ret learn(Context * ctxt, AlgInput x[], Mod mod,
                      AlgorithmData * data) override;

    virtual Ret evaluate(Context * ctxt, AlgInput xin[], Mod mod,
                         AlgorithmData * data, UInt xout[]) const override
    {
      return ls_.evaluate(ctxt, xin, mod, data, xout);
    }

    virtual Ret new_xpoint(Context *, AlgInput[], Mod, AlgorithmData *) const
    {
      return SUCCESS;
    }

    virtual Ret fill_equation(Context * ctxt,
                              AlgInput xi[], const UInt tau[], Mod mod,
                              AlgorithmData * data,
                              UInt res[]) const = 0;

    void copy_data(const AlgorithmData * data, AlgorithmData * dataout) const;

    void reset(AlgorithmData * data)
    {
      ls_.reset(data);
    }

    void init(std::size_t ncoeffs, std::size_t nsamplevars,
              unsigned * needed_coeffs, unsigned needed_size,
              LinearFitData & data, unsigned extra_eqs=2);

    Ret reset_needed(AlgorithmData * data,
                     unsigned * needed_vars, unsigned needed_size)
    {
      return ls_.reset_needed(data, needed_vars, needed_size);
    }

    void invalidate()
    {
      ls_.invalidate();
    }

    unsigned n_sample_vars() const
    {
      return ls_.ntauvars_;
    }

    std::size_t neqs() const
    {
      return ls_.neqs();
    }

    std::size_t ncoeffs() const
    {
      return ls_.nvars();
    }

    unsigned nextra_eqs()
    {
      return neqs() - ncoeffs();
    }

    const flag_t * xinfo() const
    {
      return ls_.xinfo();
    }

    std::size_t n_needed_depvars() const
    {
      return ls_.n_needed_depvars();
    }

    std::size_t n_needed_indepvars() const
    {
      return ls_.n_needed_indepvars();
    }

    const std::size_t * needed_depvars() const
    {
      return ls_.needed_depvars();
    }

    const std::size_t * needed_indepvars() const
    {
      return ls_.needed_indepvars();
    }

    std::size_t n_indep_eqs() const
    {
      return ls_.n_indep_eqs();
    }

    const std::size_t * indep_eqs() const
    {
      return ls_.indep_eqs();
    }

  private:

    struct Solver : public DenseLinearSolver {

      virtual Ret fill_matrix(Context * ctxt,
                              std::size_t n_rows, const std::size_t rows[],
                              AlgInput xi[], Mod mod,
                              AlgorithmData * data,
                              MatrixView & m) const override;

      void data_init_(std::size_t ncoeffs, std::size_t nsamplevars,
                      unsigned extra_eqs, LinearFitData & data) const;

      LinearFit * lf_;
      unsigned ntauvars_, nextra_eqs_;
    };

  private:
    Solver ls_;
  };

} // namespace fflow

#endif // FFLOW_ALG_LINEAR_FIT_HH
