#include <fflow/alg_lists.hh>

namespace fflow {

  void Chain::init(const unsigned npars[], unsigned npars_size)
  {
    nparsin.resize(npars_size);
    unsigned tot = 0;
    for (unsigned i=0; i<npars_size; ++i) {
      tot += npars[i];
      nparsin[i] = npars[i];
    }
    nparsout = tot;
  }

  Ret Chain::evaluate(Context *,
                      AlgInput xin[], Mod, AlgorithmData *,
                      UInt xout[]) const
  {
    const unsigned nps = nparsin.size();

    unsigned idx=0;
    for (unsigned i=0; i<nps; ++i) {
      const unsigned nps_i = nparsin[i];
      for (unsigned j=0; j<nps_i; ++j)
        xout[idx++] = xin[i][j];
    }

    return SUCCESS;
  }


  Ret Take::init(const unsigned npars[], unsigned npars_size,
                 std::vector<InputEl> && elems)
  {
    nparsin.resize(npars_size);
    for (unsigned i=0; i<npars_size; ++i) {
      nparsin[i] = npars[i];
    }
    nparsout = elems.size();

    for (const auto & el : elems) {
      if (!(el.list < nparsin.size() && el.el < nparsin[el.list])) {
        logerr("Indexes in take pattern out of bounds");
        return FAILED;
      }
    }

    elems_ = std::move(elems);
    return SUCCESS;
  }

  Ret Take::evaluate(Context *,
                     AlgInput xin[], Mod, AlgorithmData *,
                     UInt xout[]) const
  {
    unsigned idx=0;
    for (const auto & el : elems_)
      xout[idx++] = xin[el.list][el.el];

    return SUCCESS;
  }


  Ret Slice::init(unsigned npars, unsigned start, unsigned end)
  {
    if (end > npars || start > end)
      return FAILED;
    nparsin.resize(1);
    nparsin[0] = npars;
    start_ = start;
    end_ = end;
    nparsout = end-start;
    return SUCCESS;
  }

  Ret Slice::evaluate(Context *,
                      AlgInput xin[], Mod, AlgorithmData *,
                      UInt xout[]) const
  {
    const UInt * x = xin[0];
    std::copy(x + start_, x + end_, xout);
    return SUCCESS;
  }


  void MatrixMul::init(unsigned nrows1, unsigned ncols1, unsigned ncols2)
  {
    nr1_ = nrows1;
    nc1_ = ncols1;
    nc2_ = ncols2;

    nparsin.resize(2);
    nparsin[0] = nrows1*ncols1;
    nparsin[1] = ncols1*ncols2;
    nparsout = nrows1*ncols2;
  }

  Ret MatrixMul::evaluate(Context *,
                          AlgInput xin[], Mod mod, AlgorithmData *,
                          UInt xout[]) const
  {
    std::fill(xout, xout + nparsout, 0);
    const unsigned nrows1 = nr1_;
    const unsigned ncols1 = nc1_;
    const unsigned ncols2 = nc2_;
    const UInt * xin1 = xin[0];
    const UInt * xin2 = xin[1];
    for (unsigned i=0; i<nrows1; ++i)
      for (unsigned k=0; k<ncols1; ++k) {
        const UInt x1 = xin1[i*ncols1 + k];
        if (x1)
          for (unsigned j=0; j<ncols2; ++j) {
            const UInt x2 = xin2[k*ncols2 + j];
            if (x2) {
              UInt & xo = xout[i*ncols2 + j];
              xo = add_mod(xo, mul_mod(x1, x2, mod), mod);
            }
          }
      }
    return SUCCESS;
  }


  void Add::init(unsigned nlists, unsigned list_len)
  {
    nparsin.resize(nlists);
    std::fill(nparsin.begin(), nparsin.end(), list_len);
    nparsout = list_len;
  }

  Ret Add::evaluate(Context *,
                    AlgInput xin[], Mod mod, AlgorithmData *,
                    UInt xout[]) const
  {
    const unsigned nl = nparsin.size();
    const unsigned ne = nparsout;

    std::fill(xout, xout + ne, 0);
    for (unsigned j=0; j<nl; ++j)
      for (unsigned i=0; i<ne; ++i)
        xout[i] = add_mod(xout[i], xin[j][i], mod);

    return SUCCESS;
  }


  void Mul::init(unsigned nlists, unsigned list_len)
  {
    nparsin.resize(nlists);
    std::fill(nparsin.begin(), nparsin.end(), list_len);
    nparsout = list_len;
  }

  Ret Mul::evaluate(Context *,
                    AlgInput xin[], Mod mod, AlgorithmData *,
                    UInt xout[]) const
  {
    const unsigned nl = nparsin.size();
    const unsigned ne = nparsout;

    std::fill(xout, xout + ne, 1);
    for (unsigned j=0; j<nl; ++j)
      for (unsigned i=0; i<ne; ++i) {
        UInt xji = xin[j][i];
        if (xji != 1)
          xout[i] = mul_mod(xout[i], xin[j][i], mod);
      }

    return SUCCESS;
  }


  void NonZeroes::init(unsigned input_len)
  {
    nparsin.resize(1);
    nparsin[0] = input_len;
    nparsout = ~unsigned(0);
  }

  Ret NonZeroes::learn(Context *, AlgInput * xin, Mod, AlgorithmData *)
  {
    const unsigned nout = nparsin[0];
    const UInt * out = xin[0];
    unsigned nnz = std::count(out, out + nout, UInt(0));

    if (nparsout == ~unsigned(0)) {
      nonzero_.reset(new unsigned[nout-nnz]);
      unsigned idx=0;
      for (unsigned j=0; j<nout; ++j)
        if (out[j] != 0)
          nonzero_[idx++] = j;
    } else if (nparsout == nout-nnz) {
      unsigned idx=0;
      for (unsigned j=0; j<nout; ++j)
        if (out[j] != 0 && nonzero_[idx++] != j) {
          nparsout = ~unsigned(0);
          return FAILED;
        }
    } else {
      nparsout = ~unsigned(0);
      return FAILED;
    }

    nparsout = nout-nnz;
    return SUCCESS;
  }

  Ret NonZeroes::evaluate(Context *,
                          AlgInput xin[], Mod, AlgorithmData *,
                          UInt xout[]) const
  {
    const unsigned nout = nparsout;
    const UInt * x = xin[0];
    for (unsigned j=0; j<nout; ++j)
      xout[j] = x[nonzero_[j]];
    return SUCCESS;
  }


  Ret TakeAndAdd::init(const unsigned npars[], unsigned npars_size,
                       std::vector<std::vector<InputEl>> && elems)
  {
    nparsin.resize(npars_size);
    for (unsigned i=0; i<npars_size; ++i) {
      nparsin[i] = npars[i];
    }
    nparsout = elems.size();

    for (const auto & ellist : elems)
      for (const auto & el : ellist) {
        if (!(el.list < nparsin.size() && el.el < nparsin[el.list])) {
          logerr("Indexes in take pattern out of bounds");
          return FAILED;
        }
      }

    elems_ = std::move(elems);
    return SUCCESS;
  }

  Ret TakeAndAdd::evaluate(Context *,
                           AlgInput xin[], Mod mod,
                           AlgorithmData *,
                           UInt xout[]) const
  {
    unsigned idx=0;
    for (const auto & ellist : elems_) {
      xout[idx] = 0;
      for (const auto & el : ellist)
        xout[idx] = add_mod(xout[idx], xin[el.list][el.el], mod);
      ++idx;
    }

    return SUCCESS;
  }


  Ret TakeAndAddBL::init(const unsigned npars[], unsigned npars_size,
                         std::vector<std::vector<InputEl>> && elems)
  {
    nparsin.resize(npars_size);
    for (unsigned i=0; i<npars_size; ++i) {
      nparsin[i] = npars[i];
    }
    nparsout = elems.size();

    for (const auto & ellist : elems)
      for (const auto & el : ellist) {
        if (!(el.list1 < nparsin.size() && el.el1 < nparsin[el.list1] &&
              el.list2 < nparsin.size() && el.el2 < nparsin[el.list2])) {
          logerr("Indexes in take pattern out of bounds");
          return FAILED;
        }
      }

    elems_ = std::move(elems);
    return SUCCESS;
  }

  Ret TakeAndAddBL::evaluate(Context *,
                             AlgInput xin[], Mod mod,
                             AlgorithmData *,
                             UInt xout[]) const
  {
    unsigned idx=0;
    for (const auto & ellist : elems_) {
      xout[idx] = 0;
      for (const auto & el : ellist)
        xout[idx] = apbc_mod(xout[idx],
                             xin[el.list1][el.el1], xin[el.list2][el.el2],
                             mod);
      ++idx;
    }

    return SUCCESS;
  }


  Ret SparseMatrixMul::init(unsigned nrows1, unsigned ncols1, unsigned ncols2)
  {
    nr1_ = nrows1;
    nc1_ = ncols1;
    nc2_ = ncols2;

    nparsin.resize(2);
    nparsout = nrows1*ncols2;

    if (!(row1.size() == nrows1  && row2.size() == ncols1)) {
      logerr("Wrong number of rows in matrix multiplication");
      return FAILED;
    }

    std::size_t r1size = 0;
    for (auto & r : row1) {
      r.start = r1size;
      r1size += r.size;
      for (unsigned j=0; j<r.size; ++j)
        if (!(r.cols[j] < ncols1)) {
          logerr("Index out of bounds in matrix multiplication");
          return FAILED;
        }
    }
    nparsin[0] = r1size;

    std::size_t r2size = 0;
    for (auto & r : row2) {
      r.start = r2size;
      r2size += r.size;
      for (unsigned j=0; j<r.size; ++j)
        if (!(r.cols[j] < ncols2)) {
          logerr("Index out of bounds in matrix multiplication");
          return FAILED;
        }
    }
    nparsin[1] = r2size;

    return SUCCESS;
  }

  Ret SparseMatrixMul::evaluate(Context *,
                                AlgInput xin[], Mod mod, AlgorithmData *,
                                UInt xout[]) const
  {
    std::fill(xout, xout + nparsout, 0);
    const unsigned nrows1 = nr1_;
    //const unsigned ncols1 = nc1_;
    const unsigned ncols2 = nc2_;
    const UInt * xin1 = nullptr;
    const UInt * xin2 = nullptr;

    for (unsigned i=0; i<nrows1; ++i) {

      xin1 = xin[0] + row1[i].start;
      const unsigned * ai_col = row1[i].cols.get();
      unsigned ai_size = row1[i].size;

      for (unsigned ik=0; ik<ai_size; ++ik) {

        unsigned k = ai_col[ik];
        const UInt x1 = *xin1;
        xin2 = xin[1] + row2[k].start;
        const unsigned * bk_col = row2[k].cols.get();
        unsigned bk_size = row2[k].size;

        for (unsigned ij=0; ij<bk_size; ++ij) {
          unsigned j = bk_col[ij];
          const UInt x2 = *xin2;
          UInt & xo = xout[i*ncols2 + j];
          xo = add_mod(xo, mul_mod(x1, x2, mod), mod);
          ++xin2;
        }

        ++xin1;
      }

    }

    return SUCCESS;
  }

} // namespace fflow

