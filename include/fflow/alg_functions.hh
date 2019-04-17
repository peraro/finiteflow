#ifndef FFLOW_ALG_FUNCTIONS_HH
#define FFLOW_ALG_FUNCTIONS_HH

#include <fflow/algorithm.hh>
#include <fflow/mp_functions.hh>

namespace fflow {

  class AnalyticFunction;

  struct AnalyticFunctionData : public AlgorithmData {
    std::unique_ptr<HornerRatFunPtr[]> f;

  private:
    friend class AnalyticFunction;
    void reset_mod_(Mod mod, const AnalyticFunction & alg);

  private:
    std::unique_ptr<UInt[]> xp_ = nullptr;
    UInt this_mod_ = 0;
  };

  class AnalyticFunction : public Algorithm {
  public:

    void init(unsigned npars, unsigned nfunctions,
              AnalyticFunctionData & data);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  public:
    std::unique_ptr<MPHornerRatFunMap[]> fmap;
  };


  class EvalRationalNumbersData : public AlgorithmData {
  public:
    void resize(std::size_t size);
  private:
    friend class EvalRationalNumbers;
    SmallVector<UInt,4> val_;
    UInt mod_ = 0;
  };

  class EvalRationalNumbers : public Algorithm {
  public:

    void init(std::vector<MPRational> && vec, EvalRationalNumbersData & data);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  private:

    void init_mod_(EvalRationalNumbersData & data, Mod mod) const;

  private:
    std::vector<MPRational> val_;
  };

} // namespace fflow


#endif // FFLOW_ALG_FUNCTIONS_HH
