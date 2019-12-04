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


  void FunctionFromCoeffs::coeff_polymap_(const UInt * coeff,
                                          const CoeffHornerMap & fmap,
                                          UInt * poly)
  {
    unsigned size = fmap.size;
    for (std::size_t j=0; j<size; ++j)
      poly[fmap.pos[j]] = coeff[fmap.coeff[j]];
  }

  void FunctionFromCoeffs::coeff_ratfunmap_(const UInt * coeff,
                                            const CoeffHornerRatFunMap & fmap,
                                            HornerRatFunPtr & f)
  {
    coeff_polymap_(coeff, fmap.num_map, f.num_ptr().get());
    coeff_polymap_(coeff, fmap.den_map, f.den_ptr().get());
  }

  void FunctionFromCoeffs::init(unsigned ncoeffs, unsigned npars,
                                unsigned nfunctions,
                                FunctionFromCoeffsData & data)
  {
    nparsin.resize(2);
    nparsin[0] = ncoeffs;
    nparsin[1] = npars;
    nparsout = nfunctions;

    if (!data.f.get()) {
      data.f.reset(new HornerRatFunPtr[nfunctions]);
      fmap.reset(new CoeffHornerRatFunMap[nfunctions]);
    }
    data.xp_.reset(new UInt[npars]);
  }

  Ret FunctionFromCoeffs::evaluate(Context * ctxt,
                                   AlgInput xinin[], Mod mod,
                                   AlgorithmData * datain,
                                   UInt xout[]) const
  {
    auto & data = *static_cast<FunctionFromCoeffsData*>(datain);

    const UInt * coeff = xinin[0];
    const UInt * xin = xinin[1];
    UInt * xp = data.xp_.get();
    const unsigned nin = nparsin[1];
    auto * f = data.f.get();

    for (unsigned i=0; i<nparsout; ++i)
      coeff_ratfunmap_(coeff, fmap[i], f[i]);

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
  FunctionFromCoeffs::clone_data(const AlgorithmData * datain) const
  {
    auto & data = *static_cast<const FunctionFromCoeffsData *>(datain);
    std::unique_ptr<FunctionFromCoeffsData> ptr(new FunctionFromCoeffsData());
    auto & newalg = *ptr;
    newalg.f.reset(new HornerRatFunPtr[nparsout]);
    newalg.xp_.reset(new UInt[nparsin[1]]);
    for (unsigned i=0; i<nparsout; ++i)
      horner_ratfun_clone(data.f[i], nparsin[1], newalg.f[i]);
    return std::move(ptr);
  }



  // Rational expression

  void
  AnalyticExpression::init(unsigned npars,
                           std::vector<std::vector<Instruction>> && bytecode,
                           std::vector<MPRational> && bignums,
                           AnalyticExpressionData & data)
  {
    nparsin.resize(1);
    nparsin[0] = npars;

    bytecode_ = std::move(bytecode);
    nparsout = bytecode_.size();

    bignumbers_ = std::move(bignums);

    max_stack_size_ = compute_max_stack_size_();

    if (bignumbers_.size())
      data.bignumbers_.reset(new UInt[bignumbers_.size()]);

    data.stack_.reset(new UInt[max_stack_size_]);
  }

  AlgorithmData::Ptr
  AnalyticExpression::clone_data(const AlgorithmData *) const
  {
    typedef std::unique_ptr<AnalyticExpressionData> Ptr;
    Ptr newptr(new AnalyticExpressionData());
    if (bignumbers_.size())
      newptr->bignumbers_.reset(new UInt[bignumbers_.size()]);
    newptr->stack_.reset(new UInt[max_stack_size_]);
    return std::move(newptr);
  }

  void
  AnalyticExpression::reset_mod_(Mod mod, AnalyticExpressionData & data) const
  {
    MPInt mmod(mod.n());
    MPInt z;
    UInt * val = data.bignumbers_.get();
    unsigned idx=0;
    for (const auto & q : bignumbers_) {
      rat_mod(q, mmod, z);
      val[idx++] = z.to_uint();
    }
    data.this_mod_ = mod.n();
  }

  Ret AnalyticExpression::evaluate(Context *,
                                   AlgInput xin[], Mod mod,
                                   AlgorithmData * datain,
                                   UInt xout[]) const
  {
    AnalyticExpressionData & data = *static_cast<AnalyticExpressionData*>(datain);

    if (data.this_mod_ != mod.n())
      reset_mod_(mod, data);

    const UInt * bignums = data.bignumbers_.get();
    const UInt * in = xin[0];

#define POP() (*(--stack))
#define PUSH(el) (*(stack++) = el)

    unsigned outidx=0;

    for (const auto & bytecode : bytecode_) {

      const unsigned char * instr = bytecode.data();
      UInt * stack = data.stack_.get();

      while (*instr != END)
        switch (InstrType(*instr)) {

        case ADD:
          {
            std::size_t len = POP();
            UInt tmp = POP();
            for (unsigned j=0; j<len-1; ++j)
              tmp = add_mod(tmp, POP(), mod);
            PUSH(tmp);
            ++instr;
          }
          break;

        case MUL:
          {
            std::size_t len = POP();
            UInt tmp = POP();
            for (unsigned j=0; j<len-1; ++j)
              tmp = mul_mod(tmp, POP(), mod);
            PUSH(tmp);
            ++instr;
          }
          break;

        case NEG:
          {
            const UInt res = neg_mod(POP(), mod);
            PUSH(res);
            ++instr;
          }
          break;

        case POW:
          {
            const UInt a = POP();
            const UInt b = POP();
            const UInt res = power(a,b,mod);
            PUSH(res);
            ++instr;
          }
          break;

        case VAR:
          {
            const std::size_t idx = POP();
            const UInt res = in[idx];
            PUSH(res);
            ++instr;
          }
          break;

        case NEGPOW:
          {
            const UInt a = POP();
            const UInt b = POP();
            if (a == 0)
              return FAILED;
            const UInt res = mul_inv(power(a,b,mod),mod);
            PUSH(res);
            ++instr;
          }
          break;

        case SMALLNUM:
          {
            const UInt res = *(++instr);
            PUSH(res);
            ++instr;
          }
          break;

        case MEDNUM:
          {
            const UInt len = *(++instr);
            UInt tmp = 0;
            for (unsigned j=0; j<len; ++j) {
              tmp = tmp << 8*sizeof(Instruction);
              tmp |= *(++instr);
            }
            PUSH(tmp);
            ++instr;
          }
          break;

        case BIGNUM:
          {
            const std::size_t idx = POP();
            const UInt res = bignums[idx];
            PUSH(res);
            ++instr;
          }
          break;

        case END:
          break;

        }

      xout[outidx] = POP();
      ++outidx;

    }

#undef POP
#undef PUSH

    return SUCCESS;
  }

  unsigned AnalyticExpression::compute_max_stack_size_() const
  {
    unsigned max_stack=0;

#define POP() (--stack);
#define PUSH() (max_stack = std::max(max_stack,++stack));

    for (const auto & bytecode : bytecode_) {

      const unsigned char * instr = bytecode.data();
      unsigned stack=0;

      while (*instr != END)
        switch (InstrType(*instr)) {

        case ADD:
          {
            POP();
            POP();
            POP();
            PUSH();
            ++instr;
          }
          break;

        case MUL:
          {
            POP();
            POP();
            POP();
            PUSH();
            ++instr;
          }
          break;

        case NEG:
          {
            POP();
            PUSH();
            ++instr;
          }
          break;

        case POW:
          {
            POP();
            POP();
            PUSH();
            ++instr;
          }
          break;

        case VAR:
          {
            POP();
            PUSH();
            ++instr;
          }
          break;

        case NEGPOW:
          {
            POP();
            POP();
            PUSH();
            ++instr;
          }
          break;

        case SMALLNUM:
          {
            PUSH();
            ++instr;
            ++instr;
          }
          break;

        case MEDNUM:
          {
            const UInt len = *(++instr);
            for (unsigned j=0; j<len; ++j)
              ++instr;
            ++instr;
            PUSH();
          }
          break;

        case BIGNUM:
          {
            POP();
            PUSH();
            ++instr;
          }
          break;

        case END:
          break;

        }

#undef POP
#undef PUSH

    }

    return max_stack;
  }


} // namespace fflow
