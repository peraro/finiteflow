#include <fflow/capi_custom_threads.h>
#include <fflow/graph.hh>
using namespace fflow;

#define session (*global_session)

struct FFMTGraphEvaluator {
  GraphExtParallelEvaluator eval;
};

namespace {

  typedef GraphExtParallelEvaluator MTEval;

}

extern "C" {

  FFMTGraphEvaluator * ffMTEvalBegin(FFGraph graph, unsigned n_threads)
  {
    std::unique_ptr<MTEval> ev(new MTEval(session, graph, n_threads));
    if (!ev->valid())
      return 0;
    return (FFMTGraphEvaluator*)ev.release();
  }

  void ffMTEvalEnd(FFMTGraphEvaluator * thr_eval)
  {
    auto * ev = (MTEval*)thr_eval;
    delete ev;
  }

  FFStatus ffMTEvalWithId(FFMTGraphEvaluator * thr_eval, unsigned thread_id,
                          const FFUInt x[], FFMod mod, FFUInt * out)
  {
    Mod p;
    p.mod_ = mod;
    auto * ev = (MTEval*)thr_eval;
    Ret ret = ev->evaluate(x, p, thread_id, out);
    if (ret != SUCCESS)
      return FF_ERROR;
    return FF_SUCCESS;
  }

} // extern "C"
