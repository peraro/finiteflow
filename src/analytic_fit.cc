#include <fflow/analytic_fit.hh>
#include <fflow/graph.hh>

namespace fflow {

  namespace  {

    void mpmap_mpdata_(const std::vector<std::vector<MPHornerMap>> & cmpmap,
                       std::vector<HornerHornerRatFunMap> & cmap, Mod mod,
                       const LinearFit::flag_t * xinf = nullptr)
    {
      const auto * cf = cmpmap.data();
      for (auto & cfun : cmap) {
        if (xinf == nullptr || (*xinf++ & LSVar::IS_NON_ZERO))
          horner2_polymap_mprational(*cf, mod, cfun);
        ++cf;
      }
    }

  } // namespace


  void AnalyticFitData::reset_mod_(const AnalyticFit & alg, Mod mod)
  {
    this_mod_ = mod.n();
    MPInt tmp, mmod(mod.n());
    const LinearFit::flag_t * xinf = alg.xinfo();
    mpmap_mpdata_(alg.cmpmap, cmap, mod, xinf);
    mpmap_mpdata_(alg.integrmpmap, integrmap, mod);
    mpmap_mpdata_(alg.lexprmpmap, lexprmap, mod);
  }

  Ret AnalyticFit::init(unsigned ncoeffs, unsigned nparams,
                        unsigned ntauexpr, unsigned ntauvars,
                        unsigned needed_coeffs[], unsigned needed_size,
                        const unsigned weights_lists_len[],
                        unsigned n_weight_lists,
                        AnalyticFitData & data,
                        unsigned extra_eqs)
  {
    if (n_weight_lists) {

      unsigned tot_integr = std::accumulate(weights_lists_len,
                                            weights_lists_len + n_weight_lists,
                                            0);
      if (tot_integr != data.integr.size())
        return FAILED;

      nparsin.resize(n_weight_lists+1);
      std::copy(weights_lists_len, weights_lists_len+n_weight_lists,
                nparsin.data()+1);

      weights_.reset(new InputIntegrand[tot_integr]);
      unsigned idx=0;
      for (unsigned i=0; i<n_weight_lists; ++i)
        for (unsigned j=0; j<weights_lists_len[i]; ++j)
          weights_[idx++] = InputIntegrand{i,j};
    } else {
      nparsin.resize(1);
    }

    nparsin[0] = nparams;
    LinearFit::init(ncoeffs, ntauvars, needed_coeffs, needed_size,
                    data, extra_eqs);

    ntau = ntauvars;
    nlvar = ntauexpr ? ntauexpr : ntauvars;

    data.lvar.reset(new UInt[nlvar]);
    std::size_t shoupsize = std::max<std::size_t>({nparsin[0], nlvar, ntau});
    data.shoup.reset(new UInt[shoupsize]);
    data.lvar_shoup.reset(new UInt[nlvar]);

    return SUCCESS;
  }

  Ret AnalyticFit::new_xpoint(Context * ctxt, AlgInput xi[], Mod mod,
                              AlgorithmData * datain) const
  {
    auto & data = *static_cast<AnalyticFitData*>(datain);

    if (data.this_mod_ != mod.n())
      data.reset_mod_(*this, mod);

    UInt * ww = ctxt->ww.get();

    const UInt * x = xi[0];
    const unsigned nin = nparsin[0];
    UInt * shoup = data.shoup.get();

    const flag_t * xinf = xinfo();
    precomp_array_mul_shoup(x, nin, mod, shoup);

    HornerRatFunPtr * dest = nullptr;

    dest = data.c.data();
    for (auto & cf : data.cmap) {
      if ((*xinf++ & LSVar::IS_NON_ZERO))
        horner_ratfunmap_horner(cf, nin, ww, x, shoup, mod, *dest);
      ++dest;
    }

    dest = data.integr.data();
    for (auto & cf : data.integrmap)
      horner_ratfunmap_horner(cf, nin, ww, x, shoup, mod, *(dest++));

    dest = data.lexpr.data();
    for (auto & cf : data.lexprmap)
      horner_ratfunmap_horner(cf, nin, ww, x, shoup, mod, *(dest++));

    return SUCCESS;
  }

  Ret AnalyticFit::get_lexpr(Context * ctxt, AlgInput[],
                             const UInt tau_in[], Mod mod,
                             AlgorithmData * datain) const
  {
    auto & data = *static_cast<AnalyticFitData*>(datain);
    UInt * ww = ctxt->ww.get();

    auto & lexpr = data.lexpr;
    UInt * shoup = data.shoup.get();
    UInt * lvar = data.lvar.get();

    // get cut solutions

    if (lexpr.empty()) {

      std::copy(tau_in, tau_in + ntau, lvar);
      precomp_array_mul_shoup(lvar, ntau, mod, shoup);

    } else {

      precomp_array_mul_shoup(tau_in, ntau, mod, shoup);
      unsigned nlexpr = lexpr.size();
      for (unsigned i=0; i<nlexpr; ++i) {
        UInt tmp = lexpr[i].eval(ntau, ww, tau_in, shoup, mod);
        if (tmp == FAILED)
          return FAILED;
        lvar[i] = tmp;
      }
      precomp_array_mul_shoup(lvar, nlexpr, mod, shoup);

    }

    return SUCCESS;
  }

  Ret AnalyticFit::fill_equation(Context * ctxt,
                                 AlgInput xi[], const UInt tau[], Mod mod,
                                 AlgorithmData * datain,
                                 UInt term[]) const
  {
    auto & data = *static_cast<AnalyticFitData*>(datain);

    UInt rhs = get_lexpr(ctxt, xi, tau, mod, &data);
    if (rhs == FAILED)
      return FAILED;

    rhs = get_rhs(ctxt, xi, tau, mod, data);
    if (rhs == FAILED)
      return FAILED;

    term[ncoeffs()] = rhs;

    UInt * ww = ctxt->ww.get();
    UInt * lsubvar = data.lsubvar.get();
    UInt * lvar = data.lvar.get();
    UInt * lvar_shoup = data.lvar_shoup.get();
    UInt * shoup = data.shoup.get();
    auto * c = data.c.data();

    for (unsigned i=0; i<ncoeffs(); ++i)
      if (xinfo()[i] & LSVar::IS_NON_ZERO) {
        const unsigned * deltavarlv = deltavar[i].get();
        unsigned lsubvarsize = deltavarlv[0];
        for (unsigned j=0; j<lsubvarsize; ++j) {
          lsubvar[j] = lvar[deltavarlv[1+j]];
          lvar_shoup[j] = shoup[deltavarlv[1+j]];
        }
        UInt res = c[i].eval(lsubvarsize, ww, lsubvar, lvar_shoup, mod);
        if (res == FAILED)
          return FAILED;
        term[i] = res;
      }

    return SUCCESS;
  }

  UInt AnalyticFit::get_rhs(Context * ctxt,
                            AlgInput xin[], const UInt[],
                            Mod mod,
                            AnalyticFitData & data) const
  {
    UInt * ww = ctxt->ww.get();

    UInt tmp;
    UInt res = 0;
    AlgInput * weightsin = 0;
    const InputIntegrand * wid = weights_.get();
    if (weights_.get())
      weightsin = xin+1;

    UInt * lsubvar = data.lsubvar.get();
    UInt * lvar_shoup = data.lvar_shoup.get();
    const UInt * lvar = data.lvar.get();
    const UInt * shoup = data.shoup.get();
    const auto * c = data.c.data();

    const unsigned * deltavarlv = deltavar[ncoeffs()].get();
    unsigned lsubvarsize = deltavarlv[0];
    for (unsigned i=0; i<lsubvarsize; ++i) {
      lsubvar[i] = lvar[deltavarlv[1+i]];
      lvar_shoup[i] = shoup[deltavarlv[1+i]];
    }
    tmp = c[ncoeffs()].eval(lsubvarsize, ww, lsubvar, lvar_shoup, mod);
    if (tmp == FAILED)
      return FAILED;
    res = add_mod(res, tmp, mod);

    auto & integr = data.integr;
    const std::unique_ptr<unsigned[]> * integrlv = integrvar.get();
    for (const auto & cf : integr) {
      lsubvarsize = (*integrlv)[0];
      for (unsigned i=0; i<lsubvarsize; ++i) {
        lsubvar[i] = lvar[(*integrlv)[1+i]];
        lvar_shoup[i] = shoup[(*integrlv)[1+i]];
      }
      tmp = cf.eval(lsubvarsize, ww, lsubvar, lvar_shoup, mod);
      if (tmp == FAILED)
        return FAILED;
      if (wid) {
        tmp = mul_mod(tmp, weightsin[wid->list][wid->el], mod);
        ++wid;
      }
      res = add_mod(res, tmp, mod);
      ++integrlv;
    }

    return res;
  }

  AlgorithmData::Ptr
  AnalyticFit::clone_data(const AlgorithmData * datain) const
  {
    auto & data = *static_cast<const AnalyticFitData*>(datain);

    std::unique_ptr<AnalyticFitData> algptr(new AnalyticFitData());
    auto & alg = *algptr;

    LinearFit::copy_data(datain, algptr.get());

    alg.lvar.reset(new UInt[nlvar]);
    alg.shoup.reset(new UInt[std::max<std::size_t>(nparsin[0], nlvar)]);
    alg.lvar_shoup.reset(new UInt[nlvar]);

    if (max_subvar_size)
      alg.lsubvar.reset(new UInt[max_subvar_size]);

    unsigned csize = data.c.size();
    alg.cmap.resize(csize);
    alg.c.resize(csize);
    for (unsigned i=0; i<csize; ++i)
      horner_horner_ratfun_clone(data.cmap[i], data.c[i],
                                 deltavar[i][0], nparsin[0],
                                 alg.cmap[i], alg.c[i]);

    unsigned integrsize = data.integr.size();
    alg.integrmap.resize(integrsize);
    alg.integr.resize(integrsize);
    for (unsigned i=0; i<integrsize; ++i)
      horner_horner_ratfun_clone(data.integrmap[i], data.integr[i],
                                 integrvar[i][0], nparsin[0],
                                 alg.integrmap[i], alg.integr[i]);

    unsigned lexprsize = data.lexpr.size();
    alg.lexprmap.resize(lexprsize);
    alg.lexpr.resize(lexprsize);
    for (unsigned i=0; i<lexprsize; ++i)
      horner_horner_ratfun_clone(data.lexprmap[i], data.lexpr[i],
                                 ntau, nparsin[0],
                                 alg.lexprmap[i], alg.lexpr[i]);

    return std::move(algptr);
  }

} // namespace fflow
