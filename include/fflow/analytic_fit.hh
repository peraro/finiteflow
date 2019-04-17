#ifndef FFLOW_ALG_ANALYTIC_FIT_HH
#define FFLOW_ALG_ANALYTIC_FIT_HH

#include <fflow/alg_linear_fit.hh>
#include <fflow/mp_functions.hh>

namespace fflow {

  class AnalyticFit;
  class AnalyticFitData;

  class AnalyticFitData : public LinearFitData  {
  public:

    // tau-dependent data
    std::vector<HornerRatFunPtr> c;
    std::vector<HornerRatFunPtr> integr;
    std::vector<HornerRatFunPtr> lexpr;

    // x-dependent data (including maps of tau-dependent data)
    std::vector<HornerHornerRatFunMap> cmap;
    std::vector<HornerHornerRatFunMap> integrmap;
    std::vector<HornerHornerRatFunMap> lexprmap;

    // other variables
    std::unique_ptr<UInt[]> lvar;
    std::unique_ptr<UInt[]> lsubvar;
    std::unique_ptr<UInt[]> shoup, lvar_shoup;

  private:
    friend class AnalyticFit;

    void reset_mod_(const AnalyticFit & alg, Mod mod);

  private:
    UInt this_mod_ = 0;
  };

  class AnalyticFit : public LinearFit {
  public:

    Ret init(unsigned ncoeffs, unsigned nparams,
             unsigned ntauexpr, unsigned ntauvars,
             unsigned needed_coeffs[], unsigned needed_size,
             const unsigned weights_lists_len[],
             unsigned n_weight_lists,
             AnalyticFitData & data,
             unsigned extra_eqs=2);

    UInt get_rhs(Context * ctxt, AlgInput xi[], const UInt tau[],
                 Mod mod,
                 AnalyticFitData & data) const;

    Ret get_lexpr(Context * ctxt, AlgInput[],
                  const UInt tau_in[], Mod mod,
                  AlgorithmData * datain) const;

    virtual Ret new_xpoint(Context * ctxt, AlgInput xi[], Mod mod,
                           AlgorithmData * data) const override;
    virtual Ret fill_equation(Context * ctxt,
                              AlgInput xi[], const UInt tau[], Mod mod,
                              AlgorithmData * data,
                              UInt res[]) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  private:
    struct InputIntegrand {
      unsigned list;
      unsigned el;
    };

  public:
    // maps from multi-precision data
    std::vector<std::vector<MPHornerMap>> cmpmap;
    std::vector<std::vector<MPHornerMap>> integrmpmap;
    std::vector<std::vector<MPHornerMap>> lexprmpmap;

    // other data
    std::unique_ptr<std::unique_ptr<unsigned[]>[]> deltavar;
    std::unique_ptr<std::unique_ptr<unsigned[]>[]> integrvar;

  private:
    std::unique_ptr<InputIntegrand[]> weights_;

  public:
    UInt ntau, nlvar, max_subvar_size=0;
  };

} // namespace fflow


#endif // FFLOW_ALG_ANALYTIC_FIT_HH
