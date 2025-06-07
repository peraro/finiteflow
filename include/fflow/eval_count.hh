#include <atomic>
#include <fflow/algorithm.hh>

namespace fflow {

  class EvalCount : public Algorithm {

  public:

    void init(unsigned in_len);

    std::size_t getCount() const
    {
      return count_;
    }

    // returns current count and resets it
    std::size_t resetCount(std::size_t val = 0)
    {
      std::size_t ret = count_;
      count_ = val;
      return ret;
    }

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

  private:

    mutable std::atomic_size_t count_;
  };

}
