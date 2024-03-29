#ifndef FFLOW_MATRIX_HH
#define FFLOW_MATRIX_HH

#include <fflow/common.hh>
#include <algorithm>
#include <iostream>
#include <initializer_list>
#include <vector>
#include <array>
#include <fflow/small_vector.hh>

namespace fflow {

  // Common

  // collects bolean info on each variable
  struct LSVar {
    enum {
      IS_DEP = 1,
      IS_NEEDED = 1 << 2,
      IS_NON_ZERO = 1 << 3,
      IS_SWEEPED = 1 << 4,
      INIT = 0
    };

    typedef std::uint8_t flag_t;
  };


  // Common view for static and dynamic matrices of integers
  class MatrixView {
  public:
    typedef UInt * Row;

    MatrixView(Row * rows, unsigned nrows, unsigned ncols)
      : rows_(rows), n_(nrows), m_(ncols) {}

    UInt el(unsigned i, unsigned j) const
    {
      return rows_[i][j];
    }

    UInt & el(unsigned i, unsigned j)
    {
      return rows_[i][j];
    }

    UInt operator() (unsigned i, unsigned j) const
    {
      return el(i,j);
    }

    UInt & operator() (unsigned i, unsigned j)
    {
      return el(i,j);
    }

    const UInt * row(unsigned i) const
    {
      return rows_[i];
    }

    UInt * row(unsigned i)
    {
      return rows_[i];
    }

    MatrixView & swapRows(unsigned i, unsigned j)
    {
      cswap(rows_[i],rows_[j]);
      return *this;
    }

    unsigned nrows() const
    {
      return n_;
    }

    unsigned ncolumns() const
    {
      return m_;
    }

    MatrixView & copy(const MatrixView & oth)
    {
      for (unsigned i=0; i<n_; ++i)
        for (unsigned j=0; j<m_; ++j)
          rows_[i][j] = oth(i,j);
      return *this;
    }

    MatrixView & toRowEcholon(Mod mod);
    MatrixView & toReducedRowEcholon(Mod mod);

    // returns whether the matrix corresponds to an impossible system:
    // needs to be called after toReducedRowEcholon
    bool isImpossibleSystem();

    bool isZeroDimSystem();

    // write the independent vars in the vectors and returns 0,
    // assuming the system is not impossible (isImpossibleSystem() can
    // be called first, in order to check).
    void independent_vars(std::vector<std::size_t> & vars);

    // returns dependent variables
    void dependent_vars(std::vector<std::size_t> & vars);

  private:
    Row * rows_;
    unsigned n_, m_;
  };


  template <unsigned N, unsigned M>
  class StaticMatrix {
  public:
    typedef UInt * Row;

    StaticMatrix() : data_()
    {
      reinitialize_rows_();
    }

    explicit StaticMatrix(std::initializer_list<UInt> l) : data_()
    {
      unsigned j=0;
      for (auto i : l)
        data_[j++] = i;
      reinitialize_rows_();
    }

    UInt operator() (unsigned i, unsigned j) const
    {
      return rows_[i][j];
    }

    UInt & operator() (unsigned i, unsigned j)
    {
      return rows_[i][j];
    }

    StaticMatrix & swap_rows(unsigned i, unsigned j)
    {
      cswap(rows_[i],rows_[j]);
      return *this;
    }

    unsigned nrows() const
    {
      return N;
    }

    unsigned ncolumns() const
    {
      return M;
    }

    MatrixView getMatrixView()
    {
      return MatrixView(rows_,N,M);
    }

    MatrixView getSubMatrixView(std::size_t n, std::size_t m)
    {
      return MatrixView(rows_,n,m);
    }

  private:

    void reinitialize_rows_()
    {
      for (unsigned i=0; i<N; ++i)
        rows_[i] = data_.data() + i*M;
    }

  private:
    std::array<UInt,N*M> data_;
    Row rows_[N];
  };


  class DynamicMatrix {
  public:
    typedef UInt * Row;

    DynamicMatrix() : data_(), rows_(), n_(0), m_(0) {}

    DynamicMatrix(unsigned nrows, unsigned ncols, std::initializer_list<UInt> l)
      : data_(l), rows_(nrows), n_(nrows), m_(ncols)
    {
      reinitialize_rows_();
    }

    DynamicMatrix(unsigned nrows, unsigned ncols)
      : data_(nrows*ncols), rows_(nrows), n_(nrows), m_(ncols)
    {
      reinitialize_rows_();
    }

    DynamicMatrix & reset(unsigned nrows, unsigned ncols)
    {
      n_ = nrows;
      m_ = ncols;
      data_.clear();
      data_.resize(nrows*ncols);
      rows_.resize(nrows);
      reinitialize_rows_();
      return *this;
    }

    bool isDimension(unsigned nrows, unsigned ncols) const
    {
      return (n_ == nrows) && (m_ == ncols);
    }

    DynamicMatrix & resetRowOrder()
    {
      reinitialize_rows_();
      return *this;
    }

    UInt operator() (unsigned i, unsigned j) const
    {
      return rows_[i][j];
    }

    UInt & operator() (unsigned i, unsigned j)
    {
      return rows_[i][j];
    }

    UInt el(unsigned i, unsigned j) const
    {
      return rows_[i][j];
    }

    UInt & el(unsigned i, unsigned j)
    {
      return rows_[i][j];
    }

    DynamicMatrix & swap_rows(unsigned i, unsigned j)
    {
      cswap(rows_[i],rows_[j]);
      return *this;
    }

    unsigned nrows() const
    {
      return n_;
    }

    unsigned ncolumns() const
    {
      return m_;
    }

    MatrixView getMatrixView()
    {
      return MatrixView(rows_.data(),n_,m_);
    }

    // a matrix view for a submatrix
    MatrixView getMatrixView(std::size_t nrow, std::size_t ncols)
    {
      return MatrixView(rows_.data(),nrow,ncols);
    }

    // note: as a side effect it brings back rows in original order
    // (should not matter though)
    DynamicMatrix & resizeRows(unsigned n)
    {
      n_ = n;
      rows_.resize(n);
      data_.resize(n*m_);
      reinitialize_rows_();
      return *this;
    }

    DynamicMatrix & resizeColumns(unsigned m)
    {
      std::vector<UInt> newdata(n_*m);
      unsigned mmin = std::min(m,m_);
      for (unsigned i=0; i<n_; ++i)
        for (unsigned j=0; j<mmin; ++j)
          newdata[i*m + j] = el(i,j);
      data_.swap(newdata);
      m_ = m;
      reinitialize_rows_();
      return *this;
    }

  private:

    void reinitialize_rows_()
    {
      for (unsigned i=0; i<n_; ++i)
        rows_[i] = data_.data() + i*m_;
    }

  private:
    std::vector<UInt> data_;
    std::vector<Row> rows_;
    unsigned n_, m_;
  };



  std::ostream & operator << (std::ostream & os,
                              const MatrixView & m);

  union SparseMatCell {
    struct {
      UInt col, val;
    };
    struct {
      unsigned size, cap;
      UInt id;
    };
  };

  class SparseMatrixRow {
  public:

    enum {
      DEF_SIZE = 8,
      EL_SIZE = sizeof(SparseMatCell),
      END = fflow::FAILED
    };

    SparseMatrixRow()
      : data_(static_cast<SparseMatCell*>(std::malloc(EL_SIZE*(DEF_SIZE+2))))
    {
      data_[0].size = 0;
      data_[0].cap = DEF_SIZE;
      data_[0].id = 0;
      data_[1].col = END;
    }

    SparseMatrixRow(SparseMatrixRow && oth) : data_(oth.data_)
    {
      oth.data_ = nullptr;
    }

    void swap(SparseMatrixRow & oth)
    {
      std::swap(data_, oth.data_);
    }

    SparseMatrixRow & operator=(SparseMatrixRow && oth)
    {
      oth.swap(*this);
      return *this;
    }

    void copy_into(SparseMatrixRow & oth) const
    {
      oth.resize(size());
      std::copy(data_+1, data_+size()+2, oth.data_+1);
    }

    UInt capacity() const
    {
      return data_[0].cap;
    }

    UInt id() const
    {
      return data_[0].id;
    }

    UInt & id()
    {
      return data_[0].id;
    }

    unsigned size() const
    {
      return data_[0].size;
    }

    unsigned & size()
    {
      return data_[0].size;
    }

    UInt first_nonzero_column() const
    {
      return el(0).col;
    }

    bool resize_needed(std::size_t n) const
    {
      return n > data_[0].cap;
    }

    void resize(std::size_t n)
    {
      if (resize_needed(n)) {
        std::size_t newsize = std::max(2*data_[0].cap,unsigned(n));
        void * newdata = std::realloc(data_, EL_SIZE*(newsize+2));
        data_ = static_cast<SparseMatCell*>(newdata);
        data_[0].cap = newsize;
      }
      size() = n;
    }

    void clear()
    {
      size() = 0;
      data_[1].col = END;
    }

    void free()
    {
      if (data_) {
        std::free(data_);
        data_ = nullptr;
      }
    }

    ~SparseMatrixRow()
    {
      std::free(data_);
    }

    SparseMatCell & el(std::size_t idx)
    {
      return data_[1+idx];
    }

    const SparseMatCell & el(std::size_t idx) const
    {
      return data_[1+idx];
    }

    void fromDenseMatrixRow(const UInt * r, std::size_t ncols,
                            bool incl_id)
    {
      std::size_t nonzero = ncols - std::count(r, r+ncols, 0);
      if (resize_needed(nonzero))
        resize(nonzero);
      if (incl_id)
        data_[0].id = r[ncols];
      unsigned idx=0;
      for (unsigned i=0; i<ncols; ++i)
        if (r[i]) {
          el(idx).col = i;
          el(idx).val = r[i];
          ++idx;
        }
      el(idx).col = END;
    }

    void toDenseMatrixRow(UInt * r, std::size_t ncols,
                          bool incl_id)
    {
      std::fill(r, r+ncols, 0);
      if (incl_id)
        r[ncols] = data_[0].id;
      for (unsigned idx=0; el(idx).col != END; ++idx)
        r[el(idx).col] = el(idx).val;
    }

    // z must be nonzero
    void mul(UInt z, Mod mod)
    {
      for (unsigned idx=0; el(idx).col != END; ++idx)
        el(idx).val = mul_mod(el(idx).val, z, mod);
    }

    unsigned eq_to_substitute(const unsigned * eqs) const
    {
      const unsigned NO_EQ = ~unsigned(0);
      unsigned eq = NO_EQ;
      for (unsigned idx=0; el(idx).col != END; ++idx)
        if (eqs[el(idx).col] < eq)
          eq = eqs[el(idx).col];
      return eq;
    }

    unsigned eq_to_back_substitute(const unsigned * eqs) const
    {
      const unsigned NO_EQ = ~unsigned(0);
      unsigned eq = NO_EQ;
      for (unsigned idx=1; el(idx).col != END; ++idx)
        if (eqs[el(idx).col] < eq)
          eq = eqs[el(idx).col];
      return eq;
    }

    unsigned any_eq_to_substitute(const unsigned * eqs) const
    {
      const unsigned NO_EQ = ~unsigned(0);
      unsigned eq = NO_EQ;
      for (unsigned idx=0; el(idx).col != END; ++idx)
        if (eqs[el(idx).col] != NO_EQ)
          return eqs[el(idx).col];
      return eq;
    }

    unsigned eq_to_substitute(const unsigned * eqs, unsigned maxrow) const
    {
      const unsigned NO_EQ = ~unsigned(0);
      unsigned eq = NO_EQ;
      for (unsigned idx=0; el(idx).col < maxrow; ++idx)
        if (eqs[el(idx).col] < eq)
          eq = eqs[el(idx).col];
      return eq;
    }

    unsigned eq_to_back_substitute(const unsigned * eqs, unsigned maxrow) const
    {
      const unsigned NO_EQ = ~unsigned(0);
      unsigned eq = NO_EQ;
      for (unsigned idx=1; el(idx).col < maxrow; ++idx)
        if (eqs[el(idx).col] < eq)
          eq = eqs[el(idx).col];
      return eq;
    }

    unsigned any_eq_to_substitute(const unsigned * eqs, unsigned maxrow) const
    {
      const unsigned NO_EQ = ~unsigned(0);
      unsigned eq = NO_EQ;
      for (unsigned idx=0; el(idx).col < maxrow; ++idx)
        if (eqs[el(idx).col] != NO_EQ)
          return eqs[el(idx).col];
      return eq;
    }

    // r1 = r1-r1(pivot)*r2
    void gauss_elimination(SparseMatrixRow & r1,
                           const SparseMatrixRow & r2,
                           std::size_t pivot,
                           Mod mod);

    UInt get(UInt col) const
    {
      unsigned idx=0;
      while (el(idx).col < col)
        ++idx;
      return el(idx).col == col ? el(idx).val : 0;
    }

    // same as get but assumes col is either the first non-vanishing
    // column or not present at all
    UInt get_first(UInt col) const
    {
      return el(0).col == col ? el(0).val : 0;
    }

    UInt get_first() const
    {
      return el(0).col != END ? el(0).val : 0;
    }

    bool is_zero() const
    {
      return el(0).col == END;
    }

    void debug_print(std::ostream & os);

  private:
    SparseMatCell * data_;
  };


  inline void swap(SparseMatrixRow & r1, SparseMatrixRow & r2)
  {
    r1.swap(r2);
  }


  class SparseMatrix {
  public:

    typedef LSVar::flag_t flag_t;

    typedef SmallVector<unsigned,16> EqDeps;

    SparseMatrix() : rows_(nullptr), n_(0), m_(0) {}

    void resize(std::size_t n, std::size_t m)
    {
      if (n_ != n)
        rows_.reset(new SparseMatrixRow[n+1]);
      n_ = n;
      m_ = m;
      var_eq_.reset(new unsigned[m]);
    }

    void restrict_rows(std::size_t n)
    {
      if (n<n_) {
        for (unsigned i=n+1; i<n_; ++i)
          rows_[i].free();
        n_ = n;
      }
    }

    void fromDenseMatrix(const MatrixView & mv, bool incl_id)
    {
      resize(mv.nrows(), mv.ncolumns());
      m_ = mv.ncolumns();
      for (unsigned i=0; i<n_; ++i)
        rows_[i].fromDenseMatrixRow(mv.row(i), mv.ncolumns(), incl_id);
    }

    void toDenseMatrix(MatrixView & mv, bool incl_id)
    {
      for (unsigned i=0; i<n_; ++i)
        rows_[i].toDenseMatrixRow(mv.row(i), m_, incl_id);
    }

    void sortRows();

    void toRowEcholon(Mod mod, EqDeps * eqdeps = nullptr);
    void toReducedRowEcholon(Mod mod,
                             flag_t * flags = nullptr,
                             EqDeps * eqdeps = nullptr);
    void toReducedRowEcholon(Mod mod, EqDeps * eqdeps)
    {
      toReducedRowEcholon(mod, nullptr, eqdeps);
    }

    void toRowEcholonWithMaxRow(Mod mod, unsigned maxrow,
                                EqDeps * eqdeps = nullptr);
    void toReducedRowEcholonWithMaxRow(Mod mod,
                                       unsigned maxrow, bool reduced,
                                       flag_t * flags = nullptr,
                                       EqDeps * eqdeps = nullptr);
    void toReducedRowEcholonWithMaxRow(Mod mod, unsigned maxrow, bool reduced,
                                       EqDeps * eqdeps)
    {
      toReducedRowEcholonWithMaxRow(mod, maxrow, reduced, nullptr, eqdeps);
    }

    UInt el(unsigned i, unsigned j) const
    {
      return rows_[i].get(j);
    }

    UInt operator() (unsigned i, unsigned j) const
    {
      return el(i,j);
    }

    void swapRows(std::size_t i, std::size_t j)
    {
      rows_[i].swap(rows_[j]);
    }

    const SparseMatrixRow & row(std::size_t i) const
    {
      return rows_[i];
    }

    SparseMatrixRow & row(std::size_t i)
    {
      return rows_[i];
    }

    SparseMatrixRow & working_row()
    {
      return rows_[n_];
    }

    unsigned nrows() const
    {
      return n_;
    }

    unsigned ncolumns() const
    {
      return m_;
    }

    void debug_print(std::ostream & os);

    // returns dependent variables
    void dependent_vars(std::vector<std::size_t> & vars) const;

    bool isImpossibleSystem() const;

  private:
    enum {
      NO_EQ_ = ~unsigned(0)
    };

  private:
    std::unique_ptr<SparseMatrixRow[]> rows_;
    std::size_t n_, m_;
    std::unique_ptr<unsigned[]> var_eq_;
  };

  std::ostream & operator << (std::ostream & os,
                              const SparseMatrix & m);

} // namespace fflow


#endif // FFLOW_MATRIX_HH
