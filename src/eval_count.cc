#include <algorithm>
#include <fflow/eval_count.hh>

namespace fflow {

  void EvalCount::init(unsigned in_len)
  {
    count_ = 0;
    nparsin.resize(1);
    nparsin[0] = nparsout = in_len;
  }

  Ret EvalCount::evaluate(Context *,
                          AlgInput xin[], Mod, AlgorithmData *,
                          UInt xout[]) const
  {
    //++count_;
    count_.fetch_add(1, std::memory_order_relaxed);
    std::copy(xin[0], xin[0] + nparsin[0], xout);
    return SUCCESS;
  }

}
