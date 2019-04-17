#ifndef FFLOW_ALGORITHM_HH
#define FFLOW_ALGORITHM_HH

#include <fflow/common.hh>
#include <fflow/small_vector.hh>

namespace fflow {


  class Context;
  class Algorithm;
  class AlgorithmData;


  class AlgorithmData {
  public:

    typedef std::unique_ptr<AlgorithmData> Ptr;

    virtual ~AlgorithmData() {}
  };


  typedef const UInt * const AlgInput;
  typedef UInt * MutAlgInput;


  class Algorithm {
  public:

    typedef AlgInput Input;

    virtual Ret evaluate(Context * ctxt,
                         Input xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const = 0;

    virtual UInt min_learn_times() { return 0; }

    virtual Ret learn(Context *, Input *, Mod, AlgorithmData *)
    {
      return SUCCESS;
    }

    Ret learn_all(Context * ctxt, Mod mod,  AlgorithmData * data,
                  unsigned n_singular, unsigned n_learning);

    Ret learn_all(Context * ctxt, Mod mod,  AlgorithmData * data,
                  unsigned n_singular)
    {
      return learn_all(ctxt, mod, data, n_singular, min_learn_times());
    }

    virtual ~Algorithm() {}

    virtual AlgorithmData::Ptr clone_data(const AlgorithmData *) const
    {
      return nullptr;
    }

    void set_immutable()
    {
      flags_ |= IS_IMMUTABLE_;
    }

    // This is dangerous, so use with caution
    void set_mutable()
    {
      flags_ &= ~IS_IMMUTABLE_;
    }

    bool is_mutable() const
    {
      return !(flags_ & IS_IMMUTABLE_);
    }

    bool has_learned() const
    {
      return flags_ & HAS_LEARNED_;
    }

    void set_learned(bool val = true)
    {
      if (val)
        flags_ |= HAS_LEARNED_;
      else
        flags_ &= ~HAS_LEARNED_;
    }

    // Input algorithm don't have a valid evaluate method
    void set_as_input_alg()
    {
      flags_ |= IS_INPUT_ALG_;
      flags_ |= HAS_LEARNED_;
    }

    bool is_input_alg() const
    {
      return flags_ & IS_INPUT_ALG_;
    }

  private:

    enum {
      IS_IMMUTABLE_ = 1,
      HAS_LEARNED_ = 1 << 1,
      IS_INPUT_ALG_ = 1 << 2
    };

  public:

    SmallVector<unsigned,5> nparsin;
    unsigned nparsout = 0;

  private:
    unsigned flags_ = 0;
  };


  class NoAlgorithm : public Algorithm {
  public:

    // Calling this will always fail
    virtual Ret evaluate(Context *,
                         Input[], Mod, AlgorithmData *,
                         UInt[]) const override;
  };


  class InputAlgorithm : public Algorithm {
  public:

    InputAlgorithm() : Algorithm()
    {
      set_as_input_alg();
    }

    virtual Ret evaluate(Context *,
                         Input[], Mod, AlgorithmData *,
                         UInt[]) const override;
  };


  typedef std::unique_ptr<Algorithm> AlgorithmPtr;

} // namespace fflow

#endif // FFLOW_ALGORITHM_HH
