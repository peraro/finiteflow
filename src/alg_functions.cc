#include <fflow/alg_functions.hh>
#include <fflow/mp_gcd.hh>
#include <fflow/graph.hh>

namespace fflow {

  void AnalyticFunctionData::reset_mod_(Mod mod, const AnalyticFunction & alg)
  {
    unsigned nparsin = alg.nparsin[0];
    unsigned nparsout = alg.nparsout;
    if (xp_.get() == nullptr)
      xp_.reset(new UInt[nparsin]);
    this_mod_ = mod.n();
    for (unsigned i=0; i<nparsout; ++i)
      horner_ratfunmap_mprational(alg.fmap[i], mod, f[i]);
  }

  void AnalyticFunction::init(unsigned npars, unsigned nfunctions,
                              AnalyticFunctionData & data)
  {
    nparsin.resize(1);
    nparsin[0] = npars;
    nparsout = nfunctions;

    if (!data.f.get()) {
      data.f.reset(new HornerRatFunPtr[nfunctions]);
      fmap.reset(new MPHornerRatFunMap[nfunctions]);
    }
    data.xp_.reset(new UInt[npars]);
  }

  Ret AnalyticFunction::evaluate(Context * ctxt,
                                 AlgInput xinin[], Mod mod,
                                 AlgorithmData * datain,
                                 UInt xout[]) const
  {
    auto & data = *static_cast<AnalyticFunctionData*>(datain);

    if (mod.n() != data.this_mod_)
      data.reset_mod_(mod,*this);

    const UInt * xin = xinin[0];
    UInt * xp = data.xp_.get();
    const unsigned nin = nparsin[0];
    auto * f = data.f.get();

    precomp_array_mul_shoup(xin, nin, mod, xp);

    UInt * ww = ctxt->ww.get();
    for (unsigned i=0; i<nparsout; ++i) {
      UInt res = f[i].eval(nin, ww, xin, xp, mod);
      if (res == FAILED)
        return FAILED;
      xout[i] = res;
    }

    return SUCCESS;
  }

  AlgorithmData::Ptr
  AnalyticFunction::clone_data(const AlgorithmData * datain) const
  {
    auto & data = *static_cast<const AnalyticFunctionData *>(datain);
    std::unique_ptr<AnalyticFunctionData> ptr(new AnalyticFunctionData());
    auto & newalg = *ptr;
    newalg.f.reset(new HornerRatFunPtr[nparsout]);
    for (unsigned i=0; i<nparsout; ++i)
      horner_ratfun_clone(data.f[i], nparsin[0], newalg.f[i]);
    return std::move(ptr);
  }


  void EvalRationalNumbersData::resize(std::size_t size)
  {
    val_.resize(size);
  }

  void EvalRationalNumbers::init(std::vector<MPRational> && vec,
                                 EvalRationalNumbersData & data)
  {
    val_ = std::move(vec);
    nparsout = val_.size();
    data.resize(val_.size());
  }

  void EvalRationalNumbers::init_mod_(EvalRationalNumbersData & data,
                                      Mod mod) const
  {
    MPInt mmod(mod.n());
    MPInt z;
    UInt * val = data.val_.data();
    unsigned idx=0;
    for (const auto & q : val_) {
      rat_mod(q, mmod, z);
      val[idx++] = z.to_uint();
    }
    data.mod_ = mod.n();
  }

  Ret EvalRationalNumbers::evaluate(Context*,
                                    AlgInput*, Mod mod,
                                    AlgorithmData * datain,
                                    UInt xout[]) const
  {
    auto & data = *static_cast<EvalRationalNumbersData*>(datain);
    if (mod.n() != data.mod_)
      init_mod_(data, mod);

    const unsigned n = nparsout;
    const UInt * val = data.val_.data();
    std::copy(val, val+n, xout);

    return SUCCESS;
  }

  AlgorithmData::Ptr
  EvalRationalNumbers::clone_data(const AlgorithmData *) const
  {
    typedef std::unique_ptr<EvalRationalNumbersData> Ptr;
    Ptr newptr(new EvalRationalNumbersData());
    newptr->resize(nparsout);
    return std::move(newptr);
  }

} // namespace fflow
