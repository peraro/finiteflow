#include <sstream>
#include <fflow/mp_functions.hh>
#include <fflow/mp_gcd.hh>

namespace fflow {

  void MPReconstructedPoly::copy(const MPReconstructedPoly & oth)
  {
    op_ = oth.op_;
    mod_ = oth.mod_;

    std::size_t newsize = oth.monomials_.size();
    monomials_.resize(newsize);
    cq_.resize(newsize);
    cz_.resize(newsize);

    for (std::size_t i = 0; i<newsize; ++i) {
      monomials_[i] = op_.copy(oth.monomials_[i]);
      cq_[i] = oth.cq_[i];
      cz_[i] = oth.cz_[i];
    }
  }


  void MPReconstructedPoly::copy(const SparsePoly & oth)
  {
    op_ = oth.op();

    std::size_t newsize = oth.size();
    monomials_.resize(newsize);

    std::size_t i = 0;
    for (const auto & m : oth.monomials_)
      monomials_[i++] = op_.copy(m);

    build_from_monomials_();
  }

  void MPReconstructedPoly::build_from_monomials_()
  {
    std::size_t newsize = monomials_.size();
    cq_.resize(newsize);
    cz_.resize(newsize);

    const Monomial * m = monomials_.data();

    for (std::size_t i = 0; i<newsize; ++i) {
      cz_[i] = m[i].coeff();
      cq_[i] = MPRational(rat_rec(m[i].coeff(), mod_.to_uint()));
    }
  }

  void MPReconstructedPoly::build_from_monomials_exactly_()
  {
    std::size_t newsize = monomials_.size();
    cq_.resize(newsize);
    cz_.resize(newsize);

    const Monomial * m = monomials_.data();

    for (std::size_t i = 0; i<newsize; ++i) {
      cz_[i] = m[i].coeff();
      cq_[i] = MPRational(m[i].coeff());
    }
  }

  template <typename SparsePolyT>
  void MPReconstructedPoly::merge_(SparsePolyT && oth,
                                   const MPInt & c1, const MPInt & c2,
                                   const MPInt & mod12)
  {
    const std::size_t oth_size = oth.size();
    std::vector<Monomial> & p1(monomials_);
    const std::vector<Monomial> & p2(oth.monomials_);

    unsigned i=0, j=0;
    for (j=0; j<oth_size; ++j) {

      int cmp=0;

      while (i<p1.size() && ((cmp = op_.compare(p1[i],p2[j])) > 0)) {

        // p2 is missing an entry: assuming its 0
        chinese_remainder_from_coeffs(cz_[i], 0,
                                      c1, c2, mod12,
                                      cz_[i]);
        ++i;

      }

      if (i == p1.size())
        break;

      if (cmp == 0) {

        // normal merging: both coeff.s non-zero
        chinese_remainder_from_coeffs(cz_[i], p2[j].coeff(),
                                      c1, c2, mod12,
                                      cz_[i]);
        ++i;


      } else {

        // if it gets here, p1 is missing one coefficient and we must
        // insert p2[j] before p1[i], being the exact opposite of the
        // previous case
        Monomial m = SparsePoly::monomial_copy_(std::forward<SparsePolyT>(oth),j);
        MPInt c(m.coeff());
        chinese_remainder_from_coeffs(c, 0,
                                      c2, c1, mod12,
                                      c);
        p1.insert(p1.begin()+i,std::move(m));
        cz_.insert(cz_.begin()+i,std::move(c));
        ++i;

      }

    }

    // here either j==p2.size() or i==p1.size()

    // update remaining entries of p1, if any
    for (; i<p1.size(); ++i) {
      chinese_remainder_from_coeffs(cz_[i], 0,
                                    c1, c2, mod12,
                                    cz_[i]);
    }

    // add remaining entries of p2, if any
    for (; j<oth.size(); ++j) {
      Monomial m = SparsePoly::monomial_copy_(std::forward<SparsePolyT>(oth),j);
      MPInt c(m.coeff());
      chinese_remainder_from_coeffs(c, 0,
                                    c2, c1, mod12,
                                    c);
      p1.push_back(std::move(m));
      cz_.push_back(std::move(c));
    }

    // update cq_ and mod_
    mod_ = mod12;
    const std::size_t newsize = monomials_.size();
    cq_.resize(newsize);
    for (unsigned i=0; i<newsize; ++i)
      rat_rec(cz_[i], mod_, cq_[i]);

    // if oth was rv, clear it, because its state might be invalid
    SparsePoly::clear_if_rv_(std::forward<SparsePolyT>(oth));
  }

  void MPReconstructedPoly::merge(const SparsePoly & oth,
                                  const MPInt & c1, const MPInt & c2,
                                  const MPInt & mod12)
  {
    merge_<const SparsePoly &>(oth,c1,c2,mod12);
  }
  void MPReconstructedPoly::merge(SparsePoly && oth,
                                  const MPInt & c1, const MPInt & c2,
                                  const MPInt & mod12)
  {
    merge_<SparsePoly &&>(std::move(oth),c1,c2,mod12);
  }


  std::string MPReconstructedPoly::to_str(const std::string vars[]) const
  {
    if (monomials_.empty())
      return std::string("0");

    MemoryWriter os;

    for (std::size_t i=0; i<size(); ++i) {
      if (i)
        os << "  +  ";
      if (monomials_[i].degree() == 0) {
        cq_[i].print(os);
        continue;
      }
      if (cq_[i].cmp(MPRational(1)) != 0) {
        cq_[i].print(os);
        os << "*";
      }
      op_.print(os,monomials_[i],vars,false);
    }

    return os.str();
  }

  UInt MPReconstructedPoly::eval(const UInt x[], Mod mod) const
  {
    UInt res  = 0;
    std::size_t i=0;
    MPInt c;
    MPInt mpmod(mod.n());
    for (auto & m : monomials_) {
      rat_mod(cq_[i], mpmod, c);
      m.coeff() = c.to_uint();
      res = add_mod(res, op_.eval(m,x,mod), mod);
      ++i;
    }
    return res;
  }

  void MPReconstructedPoly::monomialsMod(Mod mod) const
  {
    std::size_t i = 0;
    MPInt c;
    MPInt mpmod(mod.n());
    for (auto & m : monomials_) {
      rat_mod(cq_[i], mpmod, c);
      m.coeff() = c.to_uint();
      ++i;
    }
  }

  HornerPtr MPReconstructedPoly::toHornerPtr(Mod mod) const
  {
    monomialsMod(mod);
    return hornerptr_from_sparse_poly(monomials_.data(),
                                      monomials_.size(),
                                      nvars());
  }

  void MPReconstructedPoly::serialize(std::string & str) const
  {
    std::ostringstream ss;
    std::size_t n_vars = op_.n_vars, n_terms = monomials_.size();

    ss << n_vars << " " << n_terms;

    for (unsigned i=0; i<n_terms; ++i) {
      ss << " ";
      ss << cq_[i] << "* ";

      const Monomial & m = monomials_[i];

      for (unsigned j=0; j<n_vars; ++j)
        ss << m.exponent(j) << " ";
    }

    str = ss.str();
  }

  void MPReconstructedPoly::unserialize(const std::string & str)
  {
    std::istringstream ss(str);
    std::size_t n_vars, n_terms;

    std::string buff;

    ss >> n_vars >> n_terms;

    op_ = MonomialOp{n_vars};

    monomials_.clear();
    monomials_.reserve(n_terms);

    cq_.clear();
    cq_.reserve(n_terms);

    cz_.clear();

    mod_ = MPInt(1);

    for (unsigned i=0; i<n_terms; ++i) {
      std::getline(ss, buff, '*');
      cq_.push_back(MPRational(buff.c_str()));

      Monomial m(n_vars);
      unsigned degree = 0;
      m.coeff() = 1;

      for (unsigned j=0; j<n_vars; ++j) {
        ss >> m.exponent(j);
        degree += m.exponent(j);
      }
      m.degree() = degree;

      monomials_.push_back(std::move(m));
    }
  }


  void MPReconstructedRatFun::serialize(std::string & str) const
  {
    std::ostringstream ss;
    std::string buff;
    num_.serialize(buff);
    ss << buff << " & ";
    buff.erase();
    den_.serialize(buff);
    ss << buff;
    str = ss.str();
  }

  void MPReconstructedRatFun::unserialize(const std::string & str)
  {
    std::istringstream ss(str);
    std::string buff;
    std::getline(ss, buff, '&');
    num_.unserialize(buff);
    std::getline(ss, buff);
    den_.unserialize(buff);
  }

  std::string MPReconstructedRatFun::to_str(const std::string vars[]) const
  {
    return format("({})/({})",num_.to_str(vars),den_.to_str(vars));
  }


  void horner_map_mprational(const MPRational * coeff, Mod mod,
                             std::size_t size, const unsigned * mpos,
                             UInt * dest)
  {
    MPInt c;
    MPInt mmod(mod.n());

    for (std::size_t j=0; j<size; ++j) {
      rat_mod(coeff[j], mmod, c);
      dest[mpos[j]] = c.to_uint();
    }
  }


  void horner_map_horner(const HornerPtr coeff[],
                         unsigned nvars, UInt * __restrict workspace,
                         const UInt * __restrict x,
                         const UInt * __restrict xp,
                         Mod mod,
                         std::size_t size,
                         const unsigned * __restrict mpos,
                         UInt * dest)
  {
    for (std::size_t j=0; j<size; ++j)
      dest[mpos[j]] = horner_eval(coeff[j].get(), nvars, workspace, x, xp, mod);
  }


  void horner_polymap_horner_copy(const HornerHornerMap & map,
                                  unsigned nsubvars,
                                  HornerHornerMap & dest)
  {
    const unsigned size = map.size;
    dest.resize(size);
    for (unsigned i=0; i<size; ++i)
      dest.coeff[i] = horner_clone(map.coeff[i].get(), nsubvars);
    std::copy(map.pos.get(), map.pos.get()+size, dest.pos.get());
  }

  void horner_horner_ratfun_clone(const HornerHornerRatFunMap & map,
                                  const HornerRatFunPtr & fun,
                                  unsigned nvars,
                                  unsigned nsubvars,
                                  HornerHornerRatFunMap & destmap,
                                  HornerRatFunPtr & destfun)
  {
    destfun.num_ptr() = horner_clone(fun.num(), nvars);
    destfun.den_ptr() = horner_clone(fun.den(), nvars);
    horner_ratfunmap_horner_copy(map, nsubvars, destmap);
  }

} // namespace fflow
