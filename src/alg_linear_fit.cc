#include <fflow/alg_linear_fit.hh>
#include <fflow/primes.hh>
#include <fflow/function_cache.hh>

namespace fflow {


  Ret LinearFit::learn(Context * ctxt, AlgInput x[], Mod mod,
                       AlgorithmData * data)
  {
    if (!is_mutable())
      return MUTABILITY_ERROR;
    Ret ret = ls_.learn(ctxt, x, mod, data);
    nparsout = ls_.nparsout;
    return ret;
  }


  void LinearFit::Solver::data_init_(std::size_t ncoeffs,
                                     std::size_t nsamplevars,
                                     unsigned extra_eqs,
                                     LinearFitData & data) const
  {
    data.tau_samples_.reset(new UInt[nsamplevars+ncoeffs+extra_eqs]);
    data.tau_.reset(new UInt[nsamplevars]);
    Mod mod = Mod(BIG_UINT_PRIMES[BIG_UINT_PRIMES_SIZE-1]);
    for (unsigned i=0; i<nsamplevars+ncoeffs+extra_eqs; ++i)
      data.tau_samples_[i] = sample_uint(OFFSET_2, i, mod);
  }


  void LinearFit::init(std::size_t ncoeffs, std::size_t nsamplevars,
                       unsigned * needed_coeffs, unsigned needed_size,
                       LinearFitData & data, unsigned extra_eqs)
  {
    ls_.data_init_(ncoeffs, nsamplevars, extra_eqs, data);
    ls_.ntauvars_ = nsamplevars;
    ls_.nextra_eqs_ = extra_eqs;
    ls_.lf_ = this;
    ls_.init(ncoeffs+extra_eqs, ncoeffs, needed_coeffs, needed_size, data);
  }


  void LinearFit::copy_data(const AlgorithmData * datain,
                            AlgorithmData * dataout) const
  {
    auto & oth = *static_cast<LinearFitData *>(dataout);
    ls_.data_init_(ncoeffs(), ls_.ntauvars_, ls_.nextra_eqs_, oth);
    ls_.copy_data(datain, dataout);
  }


  Ret LinearFit::Solver::fill_matrix(Context * ctxt,
                                     std::size_t n_rows,
                                     const std::size_t rows[],
                                     AlgInput xi[], Mod mod,
                                     AlgorithmData * datain,
                                     MatrixView & m) const
  {
    auto & data = *static_cast<LinearFitData *>(datain);
    UInt * tau = data.tau_.get();
    UInt * tau_samples = data.tau_samples_.get();

    if (lf_->new_xpoint(ctxt, xi, mod, datain) == FAILED)
      return FAILED;

    for (unsigned sol=0; sol<n_rows; ++sol) {

      for (unsigned i=0; i<ntauvars_; ++i)
        tau[i] = tau_samples[rows[sol] + i];

      if (lf_->fill_equation(ctxt, xi, tau, mod, datain, m.row(sol)) == FAILED)
        return FAILED;

    }

    return SUCCESS;
  }

} // namespace fflow
