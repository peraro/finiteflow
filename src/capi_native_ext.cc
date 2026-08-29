#include "fflow/graph.hh"
#include "fflow/mp_common.hh"
#include "fflow/capi_native_ext.h"
#include <fflow/mp_gcd.hh>
using namespace fflow;

#define session (*global_session)

extern "C" {

  struct FFContext {
    Context ctxt;
  };

}

namespace  {

  class NativeExtensionData : public AlgorithmData {
  public:
    ~NativeExtensionData();

  public:
    mutable void * mutable_data = 0;
    FFDestroyMutDataFun destructor = 0;
  };


  class NativeExtension : public Algorithm {
  public:

    virtual Ret evaluate(Context * ctxt,
                         Input xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    virtual Ret learn(Context * ctxt,
                      Input xin[], Mod mod,
                      AlgorithmData * data) override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;

    virtual UInt min_learn_times() override { return min_learn_times_val; }

    virtual ~NativeExtension();

    FFInOutInfo get_ioinfo() const
    {
      FFInOutInfo info{nparsin.data(), unsigned(nparsin.size()), nparsout};
      return info;
    }

  public:
    mutable void * shared_data = 0;
    FFAlgFunctions funs = {0,0,0,0,0};
    unsigned min_learn_times_val = 0;
    FFNativeTypeID native_id = 0;
  };


  Ret NativeExtension::evaluate(Context * ctxt,  Input xin[], Mod mod,
                                AlgorithmData * data, UInt xout[]) const
  {
    void * mut_data = static_cast<NativeExtensionData*>(data)->mutable_data;
    FFStatus ret = funs.evaluate((FFContext*)ctxt, get_ioinfo(),
                                 shared_data, mut_data,
                                 xin, mod.mod_, xout);
    if (ret != FF_SUCCESS)
      return FAILED;
    return SUCCESS;
  }


  Ret NativeExtension::learn(Context * ctxt,  Input xin[], Mod mod,
                             AlgorithmData * data)
  {
    if (!funs.learn) {
      logerr("Native algorithm did not provide a learn function.");
      return FAILED;
    }

    void * mut_data = static_cast<NativeExtensionData*>(data)->mutable_data;
    FFInInfo info = {nparsin.data(), unsigned(nparsin.size())};
    FFStatus ret = funs.learn((FFContext*)ctxt, info,
                              shared_data, mut_data,
                              &nparsout, xin, mod.mod_);
    if (ret != FF_SUCCESS)
      return FAILED;
    return SUCCESS;
  }


  AlgorithmData::Ptr
  NativeExtension::clone_data(const AlgorithmData * data) const
  {
    const auto * native_data = static_cast<const NativeExtensionData*>(data);
    void * mdata = native_data->mutable_data;
    std::unique_ptr<NativeExtensionData> ptr(new NativeExtensionData());
    auto & newalg = *ptr;

    void * new_data = 0;
    if (funs.clone_mutable_data)
      new_data = funs.clone_mutable_data(get_ioinfo(), shared_data, mdata);

    newalg.mutable_data = new_data;
    newalg.destructor = native_data->destructor;

    return std::move(ptr);
  }


  NativeExtension::~NativeExtension()
  {
    if (funs.destroy_shared_data)
      funs.destroy_shared_data(get_ioinfo(), shared_data);
  }


  NativeExtensionData::~NativeExtensionData()
  {
    if (destructor)
      destructor(mutable_data);
  }


}


extern "C" {


  void ffDbgPrint(const char * str)
  {
    dbgprint(str);
  }

  void ffLogErr(const char * str)
  {
    logerr(str);
  }



  FFNode ffAlgNative(FFGraph graph,
                     const FFNode * input_nodes, unsigned n_input_nodes,
                     unsigned nparsout, unsigned min_learn_times,
                     void * shared_data, void * mutable_data,
                     FFAlgFunctions funcs)
  {
    if (!funcs.evaluate) {
      logerr("The function evaluate cannot be NULL.");
      return FF_ERROR;
    }

    if (mutable_data && !funcs.clone_mutable_data) {
      logerr("The function clone_mutable_data must be provided "
             "if mutable data is not NULL.");
      return FF_ERROR;
    }

    if (min_learn_times && !funcs.learn) {
      logerr("The function learn must be provided "
             "if min_learn_times is not zero.");
      return FF_ERROR;
    }

    typedef NativeExtensionData Data;
    std::unique_ptr<NativeExtension> algptr(new NativeExtension());
    std::unique_ptr<Data> data(new Data());
    NativeExtension & alg = *algptr;

    if (!session.graph_exists(graph)) {
      logerr("Graph does not exist");
      return FF_ERROR;
    }

    alg.nparsin.resize(n_input_nodes);
    for (unsigned i=0; i<n_input_nodes; ++i) {
      unsigned nin = ffNodeNParsOut(graph, input_nodes[i]);
      if (nin == FF_ERROR) {
        logerr("The selected input nodes are not valid.");
        return nin;
      }
      alg.nparsin[i] = nin;
    }
    alg.nparsout = nparsout;

    Graph * g = session.graph(graph);

    alg.shared_data = shared_data;
    alg.funs = funcs;
    alg.min_learn_times_val = min_learn_times;

    data->mutable_data = mutable_data;
    data->destructor = funcs.destroy_mutable_data;

    unsigned id = g->new_node(std::move(algptr), std::move(data), input_nodes);

    if (id == ALG_NO_ID)
      return FF_ERROR;

    return id;
  }

  bool ffAlgIsNative(FFGraph graph, FFNode node)
  {
    if (!session.node_exists(graph,node))
      return false;
    const Algorithm * alg = session.node(graph,node)->algorithm();
    if (dynamic_cast<const NativeExtension*>(alg))
      return true;
    return false;
  }


  void ffEnsureWorkspaceSize(size_t n)
  {
    session.main_context()->ww.ensure_size(n);
  }

  FFUInt * ffContextWorkspace(FFContext * ctx)
  {
    return ((Context*)ctx)->ww.get();
  }

  size_t ffContextWorkspaceSize(void)
  {
    return session.main_context()->ww.size();
  }


  void ffInvalidateMutableData(FFGraph graph, FFNode node)
  {
    session.invalidate_subctxt_alg_data(graph, node);
  }

  void * ffNativeSharedData(FFGraph graph, FFNode node)
  {
    if (!session.node_exists(graph,node)) {
      logerr("The node does not exist.");
      return 0;
    }

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<NativeExtension*>(alg)) {
      NativeExtension & na = *static_cast<NativeExtension*>(alg);
      return na.shared_data;
    }

    logerr("The node is not a native extension.");
    return 0;
  }

  void * ffNativeMutableData(FFGraph graph, FFNode node)
  {
    if (!session.node_exists(graph,node)) {
      logerr("The node does not exist.");
      return 0;
    }

    Algorithm * alg = session.node(graph,node)->algorithm();

    if (dynamic_cast<NativeExtension*>(alg)) {
      AlgorithmData * alg_data = session.alg_data(graph, node);
      if (!alg_data)
        return 0;
      return static_cast<NativeExtensionData*>(alg_data)->mutable_data;
    }

    logerr("The node is not a native extension.");
    return 0;
  }


  // Native ids

  static FFNativeTypeID native_id = 0;

  FFNativeTypeID ffNewNativeTypeID(void)
  {
    return ++native_id;
  }

  FFStatus ffAssignNativeTypeID(FFGraph graph, FFNode node, FFNativeTypeID id)
  {
    if (!session.node_exists(graph,node)) {
      logerr("The node does not exist.");
      return FF_ERROR;
    }

    Algorithm * alg = session.node(graph,node)->algorithm();
    if (!alg->is_mutable()) {
      logerr("You cannot change the type ID of an immutable node.");
      return FF_ERROR;
    }

    if (dynamic_cast<NativeExtension*>(alg)) {
      NativeExtension & na = *static_cast<NativeExtension*>(alg);
      if (na.native_id != 0) {
        logerr("Re-assigning the type ID of a native node is not allowed.");
        return FF_ERROR;
      }
      na.native_id = id;
      return FF_SUCCESS;
    }

    logerr("The node is not a native extension.");
    return FF_ERROR;
  }

  bool ffHasNativeTypeID(FFGraph graph, FFNode node, FFNativeTypeID type_id)
  {
    if (type_id == 0 || !session.node_exists(graph,node))
      return false;

    const Algorithm * alg = session.node(graph,node)->algorithm();
    if (dynamic_cast<const NativeExtension*>(alg)) {
      const auto & na = *static_cast<const NativeExtension*>(alg);
      return na.native_id == type_id;
    }
    return false;
  }


  // Rational numbers

  struct FFRatNumList {
    std::unique_ptr<MPRational[]> ratnums;
    std::size_t len;
  };

  FFRatNumList * ffRatNumList(const FFCStr * ratnums, unsigned length)
  {
    std::unique_ptr<FFRatNumList> ret(new FFRatNumList());
    ret->len = length;
    ret->ratnums.reset(new MPRational[length]);
    for (unsigned j=0; j<length; ++j) {
      if (ret->ratnums[j].set(ratnums[j]) != 0)
        return 0;
    }
    return ret.release();
  }

  FFStatus ffRatNumListMod(const FFRatNumList * ratnums, FFMod mod,
                           FFUInt output[])
  {
    if (!ratnums) {
      logerr("NULL pointer to FFRatNumList in ffRatNumListMod");
      return FF_ERROR;
    }

    MPInt mmod(mod.n);
    MPInt z;
    std::size_t len = ratnums->len;
    const MPRational * r = ratnums->ratnums.get();
    for (unsigned j=0; j<len; ++j) {
      if (!rat_mod(r[j], mmod, z))
        return FF_ERROR;
      output[j] = z.to_uint();
    }

    return FF_SUCCESS;
  }

  size_t ffRatNumListLen(const FFRatNumList * ratnums)
  {
    if (!ratnums) {
      logerr("NULL pointer to FFRatNumList in ffRatNumListLen");
      return 0;
    }

    return ratnums->len;
  }

  void ffFreeRatNumList(FFRatNumList * ratnums)
  {
    delete ratnums;
  }


}
