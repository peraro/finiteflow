#include <fflow/matrix.hh>
#include <fflow/primes.hh>
#include <fflow/gcd.hh>
#include <algorithm>
#include <numeric>

namespace fflow {

  MatrixView & MatrixView::toRowEcholon(Mod mod)
  {
    unsigned r=0, p;
    Row eli, elr;
    UInt ielrk, elik, elik_mp;
    for (unsigned k=0; k<m_; ++k) {

      for (p=r; p<n_; ++p)
        if (el(p,k) != 0)
          break;

      if (p==n_)
        continue;

      if (r != p)
        swapRows(r,p);

      elr = row(r);
      ielrk = (elr[k] == 1) ? 1 : mul_inv(elr[k],mod);
      elr[k] = 1;
      for (unsigned j=k+1; j<m_; ++j)
        if (elr[j])
          elr[j] = mul_mod(elr[j], ielrk, mod);

      for (unsigned i=r+1; i<n_; ++i) {
        eli = row(i);
        elik = eli[k];
        eli[k] = 0;
        if (elik) {
          elik_mp = precomp_mul_shoup(elik, mod);
          for (unsigned j=k+1; j<m_; ++j)
            if (elr[j]) {
              eli[j] = sub_mod(eli[j],
                               mul_mod_shoup(elr[j], elik, elik_mp, mod), mod);
            }
        }
      }

      ++r;
    }
    return *this;
  }


  MatrixView & MatrixView::toReducedRowEcholon(Mod mod)
  {
    toRowEcholon(mod);

    unsigned p;
    UInt elip, elip_mp;
    Row elr, eli;
    for (int r=n_-1; r>=0; --r) {
      elr = row(r);

      for (p=0; p<m_; ++p)
        if (elr[p] != 0)
          break;

      if (p==m_)
        continue;

      for (int i=0; i<r; ++i) {
        eli = row(i);
        elip = eli[p];
        eli[p] = 0;
        if (elip) {
          elip_mp = precomp_mul_shoup(elip, mod);
          for (unsigned j=p+1; j<m_; ++j)
            if (elr[j])
              eli[j] = sub_mod(eli[j],
                               mul_mod_shoup(elr[j], elip, elip_mp, mod), mod);
        }
      }

    }
    return *this;
  }


  // Algorithm: find last non-zero r.h.s., and if the system is
  // impossible, the corresponding row should be of the form
  // {0,...,0,1}
  bool MatrixView::isImpossibleSystem()
  {
    int i=n_-1;
    for (; i>=0; --i)
      if (el(i,m_-1) != 0)
        break;

    if (i<0)
      return false;

    for (unsigned j=0; j<m_-1; ++j)
      if (el(i,j) != 0)
        return false;

    return true;
  }


  bool MatrixView::isZeroDimSystem()
  {
    const std::size_t n_vars = m_-1;
    for (unsigned j=0; j<n_; ++j)
      for (unsigned i=0; i<n_vars; ++i) {
        if (j!=i && el(j,i)!=0)
          return false;
        else if (j==i && el(j,i)!=1)
          return false;
      }
    return true;
  }


  void MatrixView::independent_vars(std::vector<std::size_t> & vars)
  {
    const std::size_t n_vars = m_-1;
    std::size_t current_col = 0;

    vars.clear();

    for (unsigned j=0; j<n_; ++j) {
      while (el(j,current_col) == 0) {
        vars.push_back(current_col);
        ++current_col;
        if (current_col == n_vars)
          return;
      }
      ++current_col;
      if (current_col == n_vars)
        return;
    }

    while (current_col < n_vars)
      vars.push_back(current_col++);
  }


  void MatrixView::dependent_vars(std::vector<std::size_t> & vars)
  {
    const std::size_t n_vars = m_-1;
    std::size_t current_col = 0;

    vars.clear();

    for (unsigned j=0; j<n_; ++j) {
      for (; current_col < n_vars; ++current_col)
        if (el(j,current_col) != 0) {
          vars.push_back(current_col++);
          break;
        }
    }
  }


  std::ostream & operator << (std::ostream & os,
                              const MatrixView & m)
  {
    for (unsigned i=0; i<m.nrows(); ++i) {
      for (unsigned j=0; j<m.ncolumns(); ++j)
        os << m(i,j) << "    ";
      os << "\n";
    }
    return os;
  }

#if 0
  void SparseMatrixRow::gauss_elimination(SparseMatrixRow & r1,
                                          const SparseMatrixRow & r2,
                                          std::size_t pivot,
                                          Mod mod)
  {
    const UInt r1p = r1.get(pivot);
    unsigned idx=0;
    if (!r1p)
      return;
    auto * r1el = &r1.el(0);
    const auto * r2el = &r2.el(0);
    resize(r1.size()+r2.size());
    while (1) {
      const UInt r1col = r1el->col;
      const UInt r2col = r2el->col;
      if (r1col < r2col) {
        auto & thisel = el(idx);
        thisel.col = r1col;
        thisel.val = r1el->val;
        ++idx;
        ++r1el;
      } else if (r1col > r2col) {
        auto & thisel = el(idx);
        thisel.col = r2col;
        thisel.val = neg_mod(mul_mod(r2el->val, r1p, mod), mod);
        ++idx;
        ++r2el;
      } else if (r1col != END) {
        if (r2col != pivot) {
          UInt res = ambc_mod(r1el->val, r2el->val, r1p, mod);
          if (res) {
            auto & thisel = el(idx);
            thisel.col = r1col;
            thisel.val = res;
            ++idx;
          }
        }
        ++r1el;
        ++r2el;
      } else {
        auto & thisel = el(idx);
        thisel.col = END;
        size() = idx;
        id() = r1.id();
        //r1.swap(*this);
        copy_into(r1);
        return;
      }
    }
  }

#else

  // TODO: add back shoup mul option
  void SparseMatrixRow::gauss_elimination(const SparseMatrixRow & r1,
                                          const SparseMatrixRow & r2,
                                          std::size_t pivot,
                                          Mod mod)
  {
    resize(r1.size()+r2.size());
#if 0
    // should never happen!
    if (!r1p)
      return;
#endif
    __restrict auto * thisel = &el(0);
    __restrict const auto * r1el = &r1.el(0);
    __restrict const auto * r2el = &r2.el(0);

    while (r1el->col < pivot) {
      thisel->col = r1el->col;
      thisel->val = r1el->val;
      ++thisel;
      ++r1el;
    }

    // Assuming r1el->col == pivot now
    const UInt r1p = r1el->val;
    ++r1el;
    ++r2el;

#if 0
    while (1) {
      const UInt r1col = r1el->col;
      const UInt r2col = r2el->col;
      if (r1col == r2col && r1col == END)
        break;
      const bool r1_less = r1col < r2col;
      const bool r2_less = r2col < r1col;
      const UInt r1val = r2_less ? 0 : r1el->val;
      const UInt r2val = r1_less ? 0 : r2el->val;
      UInt res = ambc_mod(r1val, r2val, r1p, mod);
      if (res) {
        thisel->col = r2_less ? r2col : r1col;
        thisel->val = res;
        ++thisel;
      }
      r1el += !r2_less;
      r2el += !r1_less;
    }

    {
      thisel->col = END;
      size() = thisel - el0;
      id() = r1.id();
      copy_into(r1);
      return;
    }
#endif

    while (1) {
      while (r1el->col < r2el->col) {
        thisel->col = r1el->col;
        thisel->val = r1el->val;
        ++thisel;
        ++r1el;
      }

      while(r1el->col > r2el->col) {
        thisel->col = r2el->col;
        thisel->val = neg_mod(mul_mod(r2el->val, r1p, mod), mod);
        ++thisel;
        ++r2el;
      }

      while (r1el->col == r2el->col && r1el->col != END) {
        UInt res = ambc_mod(r1el->val, r2el->val, r1p, mod);
        if (res) {
          thisel->col = r1el->col;
          thisel->val = res;
          ++thisel;
        }
        ++r1el;
        ++r2el;
      }

      if (r1el->col == r2el->col && r1el->col == END) {
        thisel->col = END;
        size() = thisel - (&el(0));
        id() = r1.id();
        return;
      }
    }

  }

#endif

  void SparseMatrixRow::debug_print(std::ostream & os)
  {
    os << "|" << data_[0].cap << "|" << data_[0].id << "|";
    for (SparseMatCell * c = data_ + 1; (*c).col != END; ++c)
      os << c->col << "," << c->val << "|";
    os << std::endl;
  }

  void SparseMatrix::sortRows()
  {
    std::vector<unsigned> ii(n_);
    std::iota(ii.begin(), ii.end(), 0);
    const SparseMatrixRow * r = rows_.get();
    std::sort(ii.begin(), ii.end(),
              [r](unsigned i, unsigned j) -> bool
              {
                int cmp = compare<UInt>(r[i].el(0).col, r[j].el(0).col);
                if (cmp)
                  return cmp > 0;
                cmp = compare(r[i].size(), r[j].size());
                if (cmp)
                  return cmp < 0;
                unsigned size=r[i].size();
                for (unsigned idx=1; idx<size; ++idx)
                  if ((cmp = compare(r[i].el(idx).col, r[j].el(idx).col)))
                    return cmp > 0;
                return false;
              });
    std::unique_ptr<SparseMatrixRow[]> newrows(new SparseMatrixRow[n_+SparseMatrix::N_WORKING_ROWS]);
    for (unsigned j=0; j<n_; ++j)
      newrows[j] = std::move(rows_[ii[j]]);
    std::swap(rows_, newrows);
  }

  void SparseMatrix::toRowEcholon(Mod mod, EqDeps * eqdeps)
  {
    std::fill(var_eq_.get(), var_eq_.get() + m_, NO_EQ_);
    UInt elif;
    for (unsigned i=0; i<n_; ++i) {
      elif = rows_[i].get_first();
      if (!elif)
        continue;

      {
        unsigned eq = NO_EQ_;
        bool first_sub = true;
        SparseMatrixRow * r = &rows_[i];
        while ((eq = r->eq_to_substitute(var_eq_.get())) != NO_EQ_) {
          if (first_sub) {
            rows_[i].copy_into(working_row(1));
            first_sub = false;
          }
          UInt first_col = rows_[eq].first_nonzero_column();
          working_row(0).gauss_elimination(working_row(1),
                                           rows_[eq],first_col,mod);
          swap_working_rows();
          r = &working_row(1);
          if (eqdeps)
            eqdeps[i].push_back(eq);
        }
        if (!first_sub)
          working_row(1).copy_into(rows_[i]);
      }

      elif = rows_[i].get_first();
      if (!elif)
        continue;

      var_eq_[rows_[i].first_nonzero_column()] = i;

      if (elif != 1) {
        UInt ielif = mul_inv(elif, mod);
        rows_[i].mul(ielif, mod);
      }
    }
  }

  void SparseMatrix::toReducedRowEcholon(Mod mod,
                                         flag_t * flags,
                                         EqDeps * eqdeps)
  {
    toRowEcholon(mod, eqdeps);

    std::fill(var_eq_.get(), var_eq_.get() + m_, NO_EQ_);

    unsigned nindep = std::remove_if(rows_.get(), rows_.get()+n_,
                                     [](const SparseMatrixRow & r)
                                     {
                                       return r.is_zero();
                                     }) - rows_.get();

    for (unsigned i=0; i<nindep; ++i)
      var_eq_[rows_[i].first_nonzero_column()] = i;

    if (nindep < 2)
      return;

    for (unsigned i=0; i<nindep-1; ++i) {
      if (flags) {
        unsigned this_var = rows_[i].first_nonzero_column();
        if (!(flags[this_var] & LSVar::IS_NEEDED))
          continue;
      }
      unsigned eq = NO_EQ_;
      bool first_sub = true;
      SparseMatrixRow * r = &rows_[i];
      while ((eq = r->eq_to_back_substitute(var_eq_.get())) != NO_EQ_) {
        if (first_sub) {
          rows_[i].copy_into(working_row(1));
          first_sub = false;
        }
        UInt first_col = rows_[eq].first_nonzero_column();
        working_row(0).gauss_elimination(working_row(1),
                                         rows_[eq],first_col, mod);
        swap_working_rows();
        r = &working_row(1);
        if (eqdeps)
          eqdeps[i].push_back(eq);
      }
      if (!first_sub)
        working_row(1).copy_into(rows_[i]);
    }
  }

  void SparseMatrix::toRowEcholonWithMaxRow(Mod mod, unsigned maxrow,
                                            EqDeps * eqdeps)
  {
    std::fill(var_eq_.get(), var_eq_.get() + m_, NO_EQ_);
    UInt elif;
    for (unsigned i=0; i<n_; ++i) {
      elif = rows_[i].get_first();
      if (!elif)
        continue;

      {
        unsigned eq = NO_EQ_;
        bool first_sub = true;
        SparseMatrixRow * r = &rows_[i];
        while ((eq = r->eq_to_substitute(var_eq_.get(), maxrow)) != NO_EQ_) {
          if (first_sub) {
            rows_[i].copy_into(working_row(1));
            first_sub = false;
          }
          UInt first_col = rows_[eq].first_nonzero_column();
          working_row(0).gauss_elimination(working_row(1),
                                           rows_[eq],first_col,mod);
          swap_working_rows();
          r = &working_row(1);
          if (eqdeps)
            eqdeps[i].push_back(eq);
        }
        if (!first_sub)
          working_row(1).copy_into(rows_[i]);
      }

      elif = rows_[i].get_first();
      if (!elif)
        continue;

      // In this version, two (or more) equations may end up having "solutions"
      // for the same unknown.  In this case, we keep the simplest equation and
      // throw away the other(s).
      //
      // TODO: fix this.  We should (optionally) keep all of them!!!
      if (var_eq_[rows_[i].first_nonzero_column()] == NO_EQ_)
        var_eq_[rows_[i].first_nonzero_column()] = i;
      else
        rows_[i].clear();

      if (elif != 1) {
        UInt ielif = mul_inv(elif, mod);
        rows_[i].mul(ielif, mod);
      }
    }
  }

  void SparseMatrix::toReducedRowEcholonWithMaxRow(Mod mod,
                                                   unsigned maxrow,
                                                   bool reduced,
                                                   flag_t * flags,
                                                   EqDeps * eqdeps)
  {
    toRowEcholonWithMaxRow(mod, maxrow, eqdeps);

    std::fill(var_eq_.get(), var_eq_.get() + m_, NO_EQ_);

    unsigned nindep = std::remove_if(rows_.get(), rows_.get()+n_,
                                     [](const SparseMatrixRow & r)
                                     {
                                       return r.is_zero();
                                     }) - rows_.get();

    for (unsigned i=0; i<nindep; ++i)
      var_eq_[rows_[i].first_nonzero_column()] = i;

    if (!reduced || nindep < 2)
      return;

    for (unsigned i=0; i<nindep-1; ++i) {
      if (flags) {
        unsigned this_var = rows_[i].first_nonzero_column();
        if (!(flags[this_var] & LSVar::IS_NEEDED))
          continue;
      }
      unsigned eq = NO_EQ_;
      bool first_sub = true;
      SparseMatrixRow * r = &rows_[i];
      while ((eq = r->eq_to_back_substitute(var_eq_.get(),maxrow)) != NO_EQ_) {
        if (first_sub) {
          rows_[i].copy_into(working_row(1));
          first_sub = false;
        }
        UInt first_col = rows_[eq].first_nonzero_column();
        working_row(0).gauss_elimination(working_row(1),
                                         rows_[eq],first_col,mod);
        swap_working_rows();
        r = &working_row(1);
        if (eqdeps)
          eqdeps[i].push_back(eq);
      }
      if (!first_sub)
        working_row(1).copy_into(rows_[i]);
    }
  }


  void SparseMatrix::dependent_vars(std::vector<std::size_t> & vars) const
  {
    vars.clear();

    for (unsigned j=0; j<n_; ++j) {
      UInt col = rows_[j].first_nonzero_column();
      if (col == SparseMatrixRow::END)
        break;
      vars.push_back(col);
    }
  }

  bool SparseMatrix::isImpossibleSystem() const
  {
    return !(var_eq_[m_-1] == NO_EQ_);
  }


  std::ostream & operator << (std::ostream & os,
                              const SparseMatrix & m)
  {
    for (unsigned i=0; i<m.nrows(); ++i) {
      for (unsigned j=0; j<m.ncolumns(); ++j)
        os << m(i,j) << "    ";
      os << "\n";
    }
    return os;
  }

  void SparseMatrix::debug_print(std::ostream & os)
  {
    for (unsigned i =0; i<n_; ++i)
      rows_[i].debug_print(os);
  }


} // namespace fflow
