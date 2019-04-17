#ifndef FFLOW_MP_FUNCTIONS_HH
#define FFLOW_MP_FUNCTIONS_HH

#include <fflow/polynomial.hh>
#include <fflow/rational_function.hh>
#include <fflow/mp_common.hh>
#include <fflow/gcd.hh>
#include <fflow/ratfun_parser.hh>
#include <fflow/shared_array.hh>

namespace fflow {

  class MPReconstructedPoly {

  public:
    MPReconstructedPoly() : monomials_(), cq_(), cz_(), mod_(1), op_{0} {}

    explicit MPReconstructedPoly(std::size_t nvars)
      : monomials_(), cq_(), cz_(), mod_(1), op_{nvars} {}

    MPReconstructedPoly(const MPReconstructedPoly & oth)
       : monomials_(), cq_(), cz_(), mod_(1), op_{0}
    {
      copy(oth);
    }

    MPReconstructedPoly(MPReconstructedPoly && oth)
       : monomials_(), cq_(), cz_(), mod_(1), op_{0}
    {
      copy(std::move(oth));
    }

    void copy(const MPReconstructedPoly & oth);

    void copy(MPReconstructedPoly && oth)
    {
      op_ = oth.op_;
      monomials_ = std::move(oth.monomials_);
      cq_ = std::move(oth.cq_);
      cz_ = std::move(oth.cz_);
    }

    // note: these two assume mod are compatible
    void copy(const SparsePoly & oth);
    void copy(SparsePoly && oth)
    {
      op_ = oth.op();
      monomials_ = std::move(oth.monomials_);
      build_from_monomials_();
    }

    void swap(MPReconstructedPoly & oth)
    {
      std::swap(op_,oth.op_);
      monomials_.swap(oth.monomials_);
      cq_.swap(oth.cq_);
      cz_.swap(oth.cz_);
      mod_.swap(oth.mod_);
    }

    MPInt & mod()
    {
      return mod_;
    }

    const MPInt & mod() const
    {
      return mod_;
    }

    std::size_t nvars() const
    {
      return op_.n_vars;
    }

    std::size_t size() const
    {
      return monomials_.size();
    }

    const Monomial & monomial(std::size_t i) const
    {
      return monomials_[i];
    }
    const MPRational & coeff(std::size_t i) const
    {
      return cq_[i];
    }

    void clear()
    {
      monomials_.clear();
      cq_.clear();
      cz_.clear();
    }

    void to_one()
    {
      clear();
      monomials_.push_back(Monomial(nvars()));
      cq_.push_back(MPRational(1));
      cz_.push_back(MPInt(1));
    }

    UInt eval(const UInt x[], Mod mod) const;

    void monomialsMod(Mod mod) const;

    HornerPtr toHornerPtr(Mod mod) const;

    // merge results from another polynomial
    void merge(const SparsePoly & oth,
               const MPInt & c1, const MPInt & c2, const MPInt & mod12);
    void merge(SparsePoly && oth,
               const MPInt & c1, const MPInt & c2, const MPInt & mod12);

    std::string to_str(const std::string vars[]) const;

    void serialize(std::string & str) const;
    void unserialize(const std::string & str);

  private:
    void build_from_monomials_();

    template <typename SparsePolyT>
    void merge_(SparsePolyT && oth,
                const MPInt & c1, const MPInt & c2, const MPInt & mod12);

    friend class MPReconstructedRatFun;

  private:
    mutable std::vector<Monomial> monomials_;
    std::vector<MPRational> cq_;
    std::vector<MPInt> cz_;
    MPInt mod_;
    MonomialOp op_;
  };


  class MPReconstructedRatFun {
  public:

    MPReconstructedRatFun() : num_(), den_() {}

    explicit MPReconstructedRatFun(std::size_t nvars)
      : num_(nvars), den_(nvars) {}

    MPReconstructedRatFun(const MPReconstructedPoly & num,
                          const MPReconstructedPoly & den)
      : num_(num), den_(den) {}

    MPReconstructedRatFun(MPReconstructedPoly && num,
                          MPReconstructedPoly && den)
      : num_(std::move(num)), den_(std::move(den)) {}

    MPReconstructedRatFun(const MPReconstructedRatFun & oth)
      : num_(), den_()
    {
      copy(oth);
    }

    MPReconstructedRatFun(MPReconstructedRatFun && oth)
      : num_(), den_()
    {
      copy(std::move(oth));
    }

    void swap(MPReconstructedRatFun & oth)
    {
      num_.swap(oth.num_);
      den_.swap(oth.den_);
    }

    MPReconstructedRatFun & operator=(const MPReconstructedRatFun & oth)
    {
      MPReconstructedRatFun newfun(oth);
      newfun.swap(*this);
      return *this;
    }

    MPReconstructedRatFun & operator=(MPReconstructedRatFun && oth)
    {
      oth.swap(*this);
      return *this;
    }

    const MPReconstructedPoly & numerator() const
    {
      return num_;
    }

    MPReconstructedPoly & numerator()
    {
      return num_;
    }

    const MPReconstructedPoly & denominator() const
    {
      return den_;
    }

    MPReconstructedPoly & denominator()
    {
      return den_;
    }

    void setMod(const MPInt & mod)
    {
      num_.mod() = mod;
      den_.mod() = mod;
    }

    void copy(const MPReconstructedRatFun & oth)
    {
      num_.copy(oth.num_);
      den_.copy(oth.den_);
    }

    void copy(MPReconstructedRatFun && oth)
    {
      num_.copy(std::move(oth.num_));
      den_.copy(std::move(oth.den_));
    }

    void copy(const SparseRationalFunction & oth)
    {
      num_.copy(oth.numerator());
      den_.copy(oth.denominator());
    }

    void copy(SparseRationalFunction && oth)
    {
      num_.copy(std::move(oth.numerator()));
      den_.copy(std::move(oth.denominator()));
    }

    UInt eval(const UInt x[], Mod mod) const
    {
      UInt num = num_.eval(x,mod);
      UInt den = den_.eval(x,mod);
      return div_mod(num, den, mod);
    }

    void toHornerRatFunPtr(HornerRatFunPtr & hfun, Mod mod) const
    {
      hfun.num_ptr() = num_.toHornerPtr(mod);
      hfun.den_ptr() = den_.toHornerPtr(mod);
    }

    void merge(const SparseRationalFunction & oth,
               const MPInt & c1, const MPInt & c2, const MPInt & mod12)
    {
      num_.merge(oth.numerator(),c1,c2,mod12);
      den_.merge(oth.denominator(),c1,c2,mod12);
    }

    void merge(SparseRationalFunction && oth,
               const MPInt & c1, const MPInt & c2, const MPInt & mod12)
    {
      num_.merge(std::move(oth.numerator()),c1,c2,mod12);
      den_.merge(std::move(oth.denominator()),c1,c2,mod12);
    }

    void to_zero()
    {
      num_.clear();
      den_.to_one();
    }

    std::string to_str(const std::string vars[]) const;

    void serialize(std::string & str) const;
    void unserialize(const std::string & str);

    Ret parse(const std::string *vars, unsigned nvars,
              const char *start, const char *end)
    {
      *this = MPReconstructedRatFun(nvars);
      return parse_rat_fun(vars, nvars, start, end,
                           num_.monomials_, num_.cq_,
                           den_.monomials_, den_.cq_);
    }

  private:
    MPReconstructedPoly num_, den_;
  };


  // Mapping MPRational to UInt in Horner

  void horner_map_mprational(const MPRational * coeff, Mod mod,
                             std::size_t size, const unsigned * mpos,
                             UInt * dest);

  struct MPHornerMap {

    std::unique_ptr<MPRational[]> coeff = nullptr;
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
      coeff.reset(new MPRational[n]);
      pos.reset(new unsigned[n]);
      size = n;
    }
  };

  struct MPHornerRatFunMap {
    MPHornerMap num_map;
    MPHornerMap den_map;

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

  inline void horner_polymap_mprational(const MPHornerMap & map, Mod mod,
                                        UInt * dest)
  {
    horner_map_mprational(map.coeff.get(), mod, map.size, map.pos.get(), dest);
  }

  inline void horner_ratfunmap_mprational(const MPHornerRatFunMap & map,
                                          Mod mod,
                                          HornerRatFunPtr & dest)
  {
    horner_polymap_mprational(map.num_map, mod, dest.num_ptr().get());
    horner_polymap_mprational(map.den_map, mod, dest.den_ptr().get());
  }

  inline void horner_ratfun_clone(const HornerRatFunPtr & fun,
                                  unsigned nvars,
                                  HornerRatFunPtr & destfun)
  {
    destfun.num_ptr() = horner_clone(fun.num(), nvars);
    destfun.den_ptr() = horner_clone(fun.den(), nvars);
  }



  // Mapping HornerPtr to UInt in Horner

  void horner_map_horner(const HornerPtr coeff[],
                         unsigned nvars, UInt * __restrict workspace,
                         const UInt * __restrict x,
                         const UInt * __restrict xp_shoup,
                         Mod mod,
                         std::size_t size,
                         const unsigned * __restrict mpos,
                         UInt * dest);

  struct HornerHornerMap {

    typedef UInt * UIntPtr;

    std::unique_ptr<HornerPtr[]> coeff = nullptr;
    std::unique_ptr<unsigned[]> pos = nullptr;
    std::size_t size = 0;

    void resize(std::size_t n)
    {
      coeff.reset(new HornerPtr[n]);
      pos.reset(new unsigned[n]);
      size = n;
    }
  };

  struct HornerHornerRatFunMap {
    HornerHornerMap num_map;
    HornerHornerMap den_map;

    void resize(std::size_t num_size, std::size_t den_size)
    {
      num_map.resize(num_size);
      den_map.resize(den_size);
    }
  };

  inline void horner_polymap_horner(const HornerHornerMap & map,
                                    unsigned nvars,
                                    UInt * __restrict workspace,
                                    const UInt * __restrict x,
                                    const UInt * __restrict xp_shoup,
                                    Mod mod,
                                    UInt * dest)
  {
    horner_map_horner(map.coeff.get(), nvars, workspace,
                      x, xp_shoup, mod, map.size, map.pos.get(), dest);
  }

  inline void horner_ratfunmap_horner(const HornerHornerRatFunMap & map,
                                      unsigned nvars,
                                      UInt * __restrict workspace,
                                      const UInt * __restrict x,
                                      const UInt * __restrict xp, // shoup
                                      Mod mod,
                                      HornerRatFunPtr & dest)
  {
    horner_polymap_horner(map.num_map, nvars, workspace, x, xp, mod,
                          dest.num_ptr().get());
    horner_polymap_horner(map.den_map, nvars, workspace, x, xp, mod,
                          dest.den_ptr().get());
  }

  void horner_polymap_horner_copy(const HornerHornerMap & map,
                                  unsigned nsubvars,
                                  HornerHornerMap & dest);

  inline void horner_ratfunmap_horner_copy(const HornerHornerRatFunMap & map,
                                           unsigned nsubvars,
                                           HornerHornerRatFunMap & dest)
  {
    horner_polymap_horner_copy(map.num_map, nsubvars, dest.num_map);
    horner_polymap_horner_copy(map.den_map, nsubvars, dest.den_map);
  }

  void horner_horner_ratfun_clone(const HornerHornerRatFunMap & map,
                                  const HornerRatFunPtr & fun,
                                  unsigned nvars,
                                  unsigned nsubvars,
                                  HornerHornerRatFunMap & destmap,
                                  HornerRatFunPtr & destfun);


  // small utility
  inline void map_mp_data(const std::vector<MPHornerMap> & map, Mod mod,
                          HornerPtr dest[])
  {
    for (auto & cmap : map)
      horner_polymap_mprational(cmap, mod, (dest++)->get());
  }

  inline void map_mp_data(const MPHornerMap * map_begin,
                          const MPHornerMap * map_end,
                          Mod mod,
                          HornerPtr dest[])
  {
    for (const MPHornerMap * cmap = map_begin; cmap != map_end; ++cmap)
      horner_polymap_mprational(*cmap, mod, (dest++)->get());
  }

  inline void
  horner2_polymap_mprational(const std::vector<MPHornerMap> & cmpmap,
                             Mod mod,
                             HornerHornerRatFunMap & cfun)
  {
    const auto & cf = cmpmap;
    unsigned num_size = cfun.num_map.size;
    map_mp_data(cf.data(), cf.data()+num_size,
                mod, cfun.num_map.coeff.get());
    unsigned den_size = cfun.den_map.size;
    map_mp_data(cf.data()+num_size, cf.data()+num_size+den_size,
                mod, cfun.den_map.coeff.get());
  }

} // namespace fflow


#endif // FFLOW_MP_FUNCTIONS_HH
