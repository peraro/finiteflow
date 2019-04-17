#include <algorithm>
#include <numeric>
#include <iostream>
#include <fflow/polynomial.hh>
#include <fflow/primes.hh>

namespace fflow {


  // Univariate

  UInt UPoly::eval(UInt x, Mod mod) const
  {
    int deg = getDegree();
    UInt res = coeff_[deg];
    for (int i=deg-1; i>=0; --i)
      res = add_mod(mul_mod(res,x,mod), coeff_[i], mod);
    return res;
  }

  void UPoly::terms(UInt x, Mod mod, UInt res[]) const
  {
    int deg = getDegree();
    res[0] = 1;
    for (int i=1; i<=deg; ++i)
      res[i] = mul_mod(res[i-1], x, mod);
  }

  bool UPoly::is_equal(const UPoly & other) const
  {
    if (getDegree() != other.getDegree())
      return false;
    return std::equal(coeff_.begin(),coeff_.end(),other.coeff_.begin());
  }

  void UPoly::shift_left(std::size_t n)
  {
    for (unsigned i=0; i<coeff_.size()-n; ++i)
      coeff_[i] = coeff_[i+n];
    for (unsigned i=coeff_.size()-n; i<coeff_.size(); ++i)
      coeff_[i] = 0;
  }

  std::string UPoly::to_str(const std::string & var) const
  {
    MemoryWriter os;
    int deg = getDegree();
    for (int i=deg; i>=0; --i) {
      if (i != deg)
        os << " + ";
      os << coeff_[i];
      if (i == 1)
        os << "*" << var;
      else if (i != 0)
        os << "*" << var << "^" << i;
    }
    return os.str();
  }


  // Multivariate

  // Sparse poly

  void SparsePoly::fromUnivariate(const UPoly & up, std::size_t var)
  {
    clear();
    monomials_.reserve(up.getDegree() + 1);
    for (int i=up.getDegree(); i>=0; --i) {
      if (up[i]==0)
        continue;
      Monomial m(nvars());
      m.coeff() = up[i];
      m.degree() = i;
      m.exponent(var) = i;
      monomials_.push_back(std::move(m));
    }
  }

  void SparsePoly::remove_zeroes_()
  {
    auto first = std::find_if(monomials_.begin(), monomials_.end(),
                              [](const Monomial & m)
                              { return m.coeff() == 0;});
    if (first != monomials_.end()) {
      for(auto i = first; ++i != monomials_.end(); )
        if ((*i).coeff() != 0)
          *first++ = std::move(*i);
      monomials_.resize(first-monomials_.begin());
    }
  }


  namespace {

    struct SparseAddOP {
      static UInt binary(UInt x, UInt y, Mod mod)
      {
        return add_mod(x, y, mod);
      }
      static UInt unary(UInt x, Mod mod)
      {
        (void)(mod);
        return x;
      }
    };

    struct SparseSubOP {
      static UInt binary(UInt x, UInt y, Mod mod)
      {
        return sub_mod(x , y, mod);
      }
      static UInt unary(UInt x, Mod mod)
      {
        return neg_mod(x, mod);
      }
    };

  } // namespace


  template <typename SparsePolyT, typename AddOP>
  void SparsePoly::add_(SparsePolyT && oth, Mod mod)
  {
    unsigned i=0, j=0;
    int cmp=0;

    std::vector<Monomial> & p1(monomials_);
    const std::vector<Monomial> & p2(oth.monomials_);
    const std::size_t oth_size = oth.size();

    for (j=0; j<oth_size; ++j) {

      while (i<p1.size() && ((cmp = op_.compare(p1[i],p2[j])) > 0))
        ++i;

      if (i == p1.size())
        break;

      if (cmp == 0) {

        p1[i].coeff() = AddOP::binary(p1[i].coeff(), p2[j].coeff(), mod);

        if (p1[i].coeff() == 0)
          p1.erase(p1.begin()+i);
        else
          ++i;

      } else {

        // if it gets here, we must insert p2[j] before p1[i]
        Monomial m = monomial_copy_(std::forward<SparsePolyT>(oth),j);
        p1.insert(p1.begin()+i,std::move(m));
        p1[i].coeff() = AddOP::unary(p1[i].coeff(),mod);
        ++i;

      }

    }

    // append remaining elements of p2
    for (; j<oth.size(); ++j) {
      Monomial m = monomial_copy_(std::forward<SparsePolyT>(oth),j);
      p1.push_back(std::move(m));
    }

    // if oth was an rv, we must clear it, because monomial_copy_
    // might have left it in an invalid state with null entries
    clear_if_rv_(std::forward<SparsePolyT>(oth));
  }



  void SparsePoly::add(const SparsePoly & oth, Mod mod)
  {
    add_<const SparsePoly &, SparseAddOP>(oth,mod);
  }
  void SparsePoly::add(SparsePoly && oth, Mod mod)
  {
    add_<SparsePoly &&, SparseAddOP>(std::move(oth),mod);
  }
  void SparsePoly::sub(const SparsePoly & oth, Mod mod)
  {
    add_<const SparsePoly &, SparseSubOP>(oth,mod);
  }
  void SparsePoly::sub(SparsePoly && oth, Mod mod)
  {
    add_<SparsePoly &&, SparseSubOP>(std::move(oth),mod);
  }


  void SparsePoly::mulUnivariateLinear(std::size_t var, UInt c, Mod mod)
  {
    // create (c)*p(x_i)
    SparsePoly oth;
    oth.copy(*this);
    oth.mul(c,mod);

    // multiply this by x_var
    for (auto & m : monomials_) {
      m.exponent(var) += 1;
      m.degree() += 1;
    }

    // add the two
    add(std::move(oth),mod);
  }


  void SparsePoly::fromDegree(std::size_t deg)
  {
    clear();
    std::size_t size = binomial(deg+nvars(), deg);
    monomials_.resize(size);

    Monomial m(nvars());
    m.coeff() = 1;
    monomials_[size-1] = std::move(m);

    const Monomial * tlower = monomials_.data() + size - 1;
    Monomial * thigher = monomials_.data() + size - 2;

    for (std::size_t r=1; r<=deg; ++r) {

      std::size_t r0 = 1;
      std::size_t j = 0;

      for (std::size_t varno=0; varno<nvars(); ++varno) {

        for (std::size_t i=0; i<r0; ++i, ++j) {
          m = op_.copy(*(tlower-i));
          m.exponent(nvars()-varno-1) += 1;
          m.degree() += 1;
          *(thigher-j) = std::move(m);
        }

        r0 *= r+varno;
        r0 /= varno + 1;

      }

      tlower = thigher;
      thigher -= j;

    }
  }


  void SparsePoly::shift(const UInt c[], Mod mod)
  {
    for (unsigned i=0; i<op_.n_vars; ++i) {

      if (c[i] == 0)
        continue;

      SparsePoly toadd(op_.n_vars);

      for (auto & t : monomials_) {

        Monomial::VarExponent a = t.exponent(i);
        UInt cc = 1;

        SparsePoly tmp(op_.n_vars);

        // (x+c)^a = \sum_j binomial(a,j) x^(a-j) c^j
        // note: j=0 is already in, hence not added
        for (int j=1; j<=a; ++j) {

          cc = mul_mod(cc, c[i], mod);

          Monomial m = t.copy(op_.n_vars);
          m.degree() -= j;
          m.exponent(i) -= j;
          m.coeff() = mul_mod(m.coeff(),
                              mul_mod(cc, binomial(a,j,mod), mod), mod);

          // no need to use add: this way they are correctly ordered
          if (m.coeff())
            tmp.monomials_.push_back(std::move(m));

        }

        if (!tmp.is_zero())
          toadd.add(std::move(tmp),mod);

      }

      add(std::move(toadd),mod);

    }
  }


  UInt SparsePoly::eval(const UInt x[], Mod mod) const
  {
    UInt res  = 0;
    for (const auto & m : monomials_)
      res = add_mod(res, op_.eval(m,x,mod), mod);
    return res;
  }


  UInt SparsePoly::evalDegree(unsigned deg, const UInt x[], Mod mod) const
  {
    const int d = deg;
    UInt res  = 0;

    std::vector<Monomial>::const_iterator m = monomials_.begin();

    // find first monomial with degree=m
    for (; m != monomials_.end(); ++m)
      if ((*m).degree() == d)
        break;

    // evaluate until degree changes
    for (; m != monomials_.end() && (*m).degree() == d; ++m)
      res = add_mod(res, op_.eval(*m,x,mod), mod);

    return res;
  }


  void SparsePoly::homogenize(std::size_t var, std::size_t degree)
  {
    const Monomial::VarExponent deg = Monomial::VarExponent(degree);
    for (auto & m : monomials_) {
      m.exponent(var) += deg-m.degree();
      m.degree() = deg;
    }
    sort_();
  }


  void MonomialOp::print(MemoryWriter & os, const Monomial & m,
                         const std::string vars[],
                         bool print_coeff) const
  {
    if (m.coeff() == 0) {
      os << "0";
      return;
    }

    bool first_term = false;

    // print coefficient
    if (print_coeff && m.coeff() != 1) {
      os << m.coeff();
      first_term = true;
    }

    // print monomial
    if (m.degree()) {
      const Monomial::VarExponent * exp = m.exponents();
      for (std::size_t i=0; i<n_vars; ++i)
        if (exp[i] != 0) {
          if (first_term)
            os << "*";
          else
            first_term = true;
          os << vars[i];
          if (exp[i] != 1)
            os << "^" << int(exp[i]);
        }
    }
  }

  std::string SparsePoly::to_str(const std::string vars[]) const
  {
    if (is_zero())
      return std::string("0");
    MemoryWriter os;
    op_.print(os,monomials_[0],vars);
    for (std::size_t i=1; i<size(); ++i) {
      os << "  +  ";
      op_.print(os,monomials_[i],vars);
    }
    return os.str();
  }


  // Newton univariate

  UInt NewtonUPoly::eval(UInt x, Mod mod) const
  {
    int deg = getDegree();
    UInt res = ai_[deg];
    for (int i=deg-1; i>=0; --i)
      res = add_mod(mul_mod(res, sub_mod(x, getXi(i), mod), mod),
                    ai_[i], mod);
    return res;
  }


  UInt NewtonUPoly::evalNewtonSystemRow(UInt x, Mod mod, unsigned row,
                                        UInt & eta_row) const
  {
    int deg = row-1;
    UInt res = 0;
    UInt term = 1;
    for (int i=0; i<=deg; ++i) {
      addmul_mod(res, ai_[i], term, mod);
      term = mul_mod(term, sub_mod(x, getXi(i), mod), mod);
    }
    eta_row = term;
    return res;
  }


  void NewtonUPoly::toStdRep(Mod mod, UInt coeffs[]) const
  {
    int deg = getDegree();
    coeffs[deg] = ai_[deg];
    for (int i=deg-1; i>=0; --i) {
      coeffs[i] = ai_[i];
      for (int j=i; j<deg; ++j)
        addmul_mod(coeffs[j], neg_mod(getXi(i),mod), coeffs[j+1], mod);
    }
  }


  void NewtonUPoly::toSparsePolyVar(unsigned first_var,
                                    Mod mod,
                                    SparsePoly & res) const
  {
    UPoly p(getDegree());
    toStdRep(mod,p.coeffs());
    res.fromUnivariate(p,first_var);
  }


  // Newton multivariate

  unsigned NewtonMPoly::nvars() const
  {
    if (ai_.empty()) // this should not happen
      return 0;

    unsigned n = 2;
    const NewtonPoly * p = ai_.front().get();

    while (!p->is_univariate()) {
      ++n;
      p = static_cast<const NewtonMPoly *>(p)->ai_.front().get();
    }

    return n;
  }


  UInt NewtonMPoly::eval(UInt x, UInt xi[], Mod mod) const
  {
    unsigned deg = ai_.size()-1;
    UInt res = (*ai_[deg]).eval(xi[0],xi+1,mod);
    for (int i=deg-1; i>=0; --i)
      res = add_mod(mul_mod(res, sub_mod(x, getXi(i), mod), mod),
                    (*ai_[i]).eval(xi[0],xi+1,mod),
                    mod);
    return res;
  }


  void NewtonMPoly::toSparsePolyVar(unsigned first_var,
                                    Mod mod,
                                    SparsePoly & res) const
  {
    int deg = getDegree();

    // res = a_deg
    writeSparsePolyVar(*ai_[deg],first_var+1,mod,res);

    for (int i=deg-1; i>=0; --i) {

      // tmp = a_i
      SparsePoly tmp(res.nvars());
      writeSparsePolyVar(*ai_[i],first_var+1,mod,tmp);

      // res *= (x-x_i)
      res.mulUnivariateLinear(first_var, neg_mod(getXi(i),mod), mod);

      // res += tmp
      res.add(std::move(tmp),mod);
    }

  }


  // Horner

  UInt horner_eval(const UInt * __restrict c, unsigned nvars,
                   UInt * __restrict workspace,
                   const UInt * __restrict x,
                   const UInt * __restrict xp,
                   Mod mod)
  {
    UInt tmp;
    UInt * ww = workspace;

    // First compute the univariate polynomials from the coefficients
    {
      UInt xv = x[nvars-1];
      UInt xvp = xp[nvars-1];
      while (*c != HORNER_END) {

        std::size_t n = *c;
        ++c;

        tmp = 0;
        for (unsigned i=0; i<n; ++i, ++c)
          tmp = add_mod(mul_mod_shoup(tmp, xv, xvp, mod), *c, mod);
        *ww = tmp;
        ++ww;

      }
    }

    // Then, recursively, the polynomials in all the other variables
    for (int var=int(nvars)-2; var>=0; --var) {

      const UInt * wr = workspace;
      ww = workspace;
      ++c;

      UInt xv = x[var];
      UInt xvp = xp[var];
      while (*c != HORNER_END) {
        std::size_t n = *c;
        ++c;
        tmp = 0;
        if (!n) {
          ++wr;
        } else {
          for (unsigned i=0; i<n; ++i, ++wr)
            tmp = add_mod(mul_mod_shoup(tmp, xv, xvp, mod), *wr, mod);
        }
        *ww = tmp;
        ++ww;
      }

    }

    return *workspace;
  }

  std::size_t horner_size(const UInt * __restrict cin, unsigned nvars)
  {
    const UInt * c = cin;
    {
      while (*c != HORNER_END) {
        std::size_t n = *c;
        c += n+1;
      }
    }

    for (int var=int(nvars)-2; var>=0; --var) {
      ++c;
      while (*c != HORNER_END)
        ++c;
    }

    return c + 1 - cin;
  }

  HornerPtr horner_clone(const UInt * __restrict horner_data, unsigned nvars)
  {
    if (!horner_data)
      return HornerPtr(nullptr);
    std::size_t size = horner_size(horner_data, nvars);
    HornerPtr res(new UInt[size]);
    std::copy(horner_data, horner_data+size, res.get());
    return res;
  }


  // The required workspace size is the number of univariate (sub-)polynomials
  std::size_t horner_required_workspace(const UInt * c)
  {
    std::size_t size = 0;
    while (*c != HORNER_END) {
      std::size_t n = *c;
      ++c;
      ++size;
      c += n;
    }
    return size;
  }


  namespace  {

    std::size_t horner_get_coeffs_implem_(const UInt * c0,
                                          UInt * c,
                                          const Monomial * m,
                                          const std::size_t * pos_a,
                                          const std::size_t * pos_b,
                                          unsigned nvars,
                                          unsigned first_var,
                                          unsigned current_var,
                                          unsigned * map_pos)
    {
      std::size_t data_size = 0;
      bool last_var = first_var == nvars-1;
      bool this_var = first_var == current_var;

      if (!last_var && !this_var) {
        std::size_t this_size;

        // check zero case first
        if (pos_a == pos_b)
          return 0;

        const std::size_t * pos1 = pos_a;
        int exp = m[*pos_a].exponent(first_var);

        while (pos1 != pos_b) {

          const std::size_t * pos2;
          pos2 = std::find_if(pos1, pos_b,
                              [m, exp, first_var] (std::size_t p) -> bool
                              {
                                return m[p].exponent(first_var) != exp;
                              });
          this_size = horner_get_coeffs_implem_(c0, c, m, pos1, pos2,
                                                nvars, first_var+1,
                                                current_var,
                                                map_pos);
          if (c)
            c += this_size;
          data_size += this_size;

          // fill gaps with zeros
          int nextp = (pos2 == pos_b) ? 0 : m[*pos2].exponent(first_var);
          int gaps = (pos2 == pos_b) ? exp : exp-nextp-1;
          for (int i=0; i<gaps; ++i) {
            if (c)
              *(c++) = 0;
            ++data_size;
          }
          exp = nextp;
          pos1 = pos2;
        }

      } else if (!last_var && this_var) {

        // check zero case first
        if (pos_a == pos_b) {
          if (c)
            *c = 0;
          return 1;
        }

        // encode size from degree
        unsigned exp = m[*pos_a].exponent(current_var);
        if (c) {
          *c = exp + 1;
          ++c;
        }
        ++data_size;

      } else {

        // check zero case first
        if (pos_a == pos_b) {
          if (c)
            *c = 0;
          return 1;
        }

        // encode size from degree
        unsigned exp = m[*pos_a].exponent(nvars-1);
        if (c) {
          *c = exp + 1;
          ++c;
        }
        ++data_size;

        // fill the data
        if (c) {
          for (int i=exp; i>=0; --i, ++c) {
            if (pos_a != pos_b && m[*pos_a].exponent(nvars-1) == i) {
              *c = m[*pos_a].coeff();
              if (map_pos)
                map_pos[*pos_a] = c-c0;
              ++pos_a;
            } else {
              *c = 0;
            }
          }
        }
        data_size += exp+1;
      }

      return data_size;
    }

  } // namespace



  std::size_t horner_from_sparse_poly(UInt * c,
                                      const Monomial * m,
                                      std::size_t size,
                                      unsigned nvars,
                                      unsigned first_var,
                                      unsigned * map_pos)
  {
    std::size_t data_size = 0;

    std::vector<std::size_t> sorted(size);
    std::iota(sorted.begin(), sorted.end(), 0);

    std::sort(sorted.begin(), sorted.end(),
              [m,nvars] (std::size_t i, std::size_t j) -> bool
              {
                return MonomialCompareLex{nvars}(m[i],m[j]) > 0;
              });

    const std::size_t * pos_a = sorted.data();
    const std::size_t * pos_b = sorted.data() + size;
    const UInt * c0 = c;

    for (int var = int(nvars)-1; var>=int(first_var); --var) {
      std::size_t this_data_size;
      this_data_size = horner_get_coeffs_implem_(c0, c, m, pos_a, pos_b,
                                                 nvars, first_var, var,
                                                 map_pos);
      if (c) {
        c += this_data_size;
        *(c++) = HORNER_END;
      }
      data_size += this_data_size + 1;
    }

    return data_size;
  }


  HornerPtr hornerptr_from_sparse_poly(const Monomial * m,
                                       std::size_t size,
                                       unsigned nvars,
                                       unsigned first_var,
                                       unsigned * map_pos)
  {
    std::size_t data_size = horner_from_sparse_poly(nullptr, m, size, nvars,
                                                    first_var);
    HornerPtr ptr(new UInt[data_size]);
    horner_from_sparse_poly(ptr.get(), m, size, nvars, first_var, map_pos);
    return ptr;
  }

} // namespace fflow
