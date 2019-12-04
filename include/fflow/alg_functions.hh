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


  class FunctionFromCoeffs;

  struct CoeffHornerMap {

    std::unique_ptr<unsigned[]> coeff = nullptr;
    std::unique_ptr<unsigned[]> pos = nullptr;
    std::size_t size = 0;

    void clear()
    {
      coeff.reset(nullptr);
      pos.reset(nullptr);
      size = 0;
    }

    void resize(std::size_t n)
    {
      coeff.reset(new unsigned[n]);
      pos.reset(new unsigned[n]);
      size = n;
    }
  };

  struct CoeffHornerRatFunMap {
    CoeffHornerMap num_map;
    CoeffHornerMap den_map;

    void clear()
    {
      num_map.clear();
      den_map.clear();
    }

    void resize(std::size_t num_size, std::size_t den_size)
    {
      num_map.resize(num_size);
      den_map.resize(den_size);
    }
  };

  struct FunctionFromCoeffsData : public AlgorithmData {
    std::unique_ptr<HornerRatFunPtr[]> f;

  private:
    friend class FunctionFromCoeffs;

  private:
    std::unique_ptr<UInt[]> xp_ = nullptr;
  };

  class FunctionFromCoeffs : public Algorithm {
  public:

    void init(unsigned ncoeffs, unsigned npars, unsigned nfunctions,
              FunctionFromCoeffsData & data);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  private:

    static void coeff_polymap_(const UInt * coeff,
                               const CoeffHornerMap & fmap,
                               UInt * poly);

    static void coeff_ratfunmap_(const UInt * coeff,
                                 const CoeffHornerRatFunMap & fmap,
                                 HornerRatFunPtr & f);

  public:
    std::unique_ptr<CoeffHornerRatFunMap[]> fmap;
  };


  class AnalyticExpression;
  typedef unsigned char Instruction;

  struct AnalyticExpressionData : public AlgorithmData {

  private:
    friend class AnalyticExpression;

  private:
    std::unique_ptr<UInt[]> bignumbers_ = nullptr;
    std::unique_ptr<UInt[]> stack_ = nullptr;
    UInt this_mod_ = 0;
  };

  class AnalyticExpression : public Algorithm {
  public:

    // bytecode instructions
    enum InstrType {
          ADD,
          MUL,
          NEG,
          POW,
          VAR,
          NEGPOW,
          SMALLNUM,
          MEDNUM,
          BIGNUM,
          END
    };

  public:

    void init(unsigned npars,
              std::vector<std::vector<Instruction>> && bytecode,
              std::vector<MPRational> && bignums,
              AnalyticExpressionData & data);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

  private:

    void reset_mod_(Mod mod, AnalyticExpressionData & data) const;
    unsigned compute_max_stack_size_() const;

  private:
    std::vector<std::vector<Instruction>> bytecode_;
    std::vector<MPRational> bignumbers_;
    unsigned max_stack_size_ = 0;
  };

} // namespace fflow


#endif // FFLOW_ALG_FUNCTIONS_HH
