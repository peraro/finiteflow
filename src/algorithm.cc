#include <fflow/algorithm.hh>
#include <fflow/function_cache.hh>

namespace fflow {

  Ret Algorithm::learn_all(Context * ctxt, Mod mod,  AlgorithmData * data,
                           unsigned n_singular, unsigned n_learning)
  {
    if (n_learning == 0) {
      set_learned();
      return SUCCESS;
    }

    set_learned(false);

    unsigned npin = nparsin[0];

    std::unique_ptr<UInt[]> xptr(new UInt[npin]);
    UInt * x = xptr.get();

    unsigned learn_ok = 0;
    unsigned singular = 0;
    unsigned try_no = 0;
    const std::size_t max_tries = n_learning + n_singular;

    for (try_no=0; try_no<max_tries; ++try_no) {

      for (unsigned i=0; i<npin; ++i)
        x[i] = sample_uint(OFFSET_5, try_no + i, mod);

      Ret ret = learn(ctxt, &x, mod, data);
      if (ret == SUCCESS) {
        ++learn_ok;
        if (learn_ok >= n_learning) {
          set_learned();
          return SUCCESS;
        }
      } else {
        learn_ok = 0;
        ++singular;
        if (singular > n_singular)
          return ret;
      }

    }

    return FAILED;
  }


  Ret NoAlgorithm::evaluate(Context *,
                            Input[], Mod, AlgorithmData *,
                            UInt[]) const
  {
    logerr("Error: NoAlgorithm::evaluate should never be called.");
    return FAILED;
  }


  Ret InputAlgorithm::evaluate(Context *,
                               Input[], Mod, AlgorithmData *,
                               UInt[]) const
  {
    logerr("Internal error: Graph::GraphInputVars_::evaluate "
           "should never be called.");
    return FAILED;
  }

} // namespace fflow


