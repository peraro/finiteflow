// Basic algorithm for manipulating lists of arguments.  These do not
// need any AnlgorithmData, so we can just pass a nullptr for that

#ifndef FFLOW_ALG_LISTS_HH
#define FFLOW_ALG_LISTS_HH

#include <fflow/algorithm.hh>

namespace fflow {

  class Chain : public Algorithm {
  public:

    void init(const unsigned npars[], unsigned npars_size);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;
  };


  class Take : public Algorithm {
  public:

    struct InputEl {
      unsigned list;
      unsigned el;
    };

    Ret init(const unsigned npars[], unsigned npars_size,
             std::vector<InputEl> && elems);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

  private:
    std::vector<InputEl> elems_;
  };


  class Slice : public Algorithm {
  public:

    Ret init(unsigned npars, unsigned start, unsigned end);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

  private:
    unsigned start_=0, end_=0;
  };


  class MatrixMul : public Algorithm {
  public:

    void init(unsigned nrows1, unsigned ncols1, unsigned ncols2);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

  private:
    unsigned nr1_, nc1_, nc2_;
  };


  class Add : public Algorithm {
  public:

    void init(unsigned nlists, unsigned list_len);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;
  };


  class Mul : public Algorithm {
  public:

    void init(unsigned nlists, unsigned list_len);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;
  };


  class NonZeroes : public Algorithm {
  public:

    void init(unsigned input_len);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    virtual Ret learn(Context * ctxt, AlgInput * xin, Mod mod,
                      AlgorithmData * data) override;

    virtual UInt min_learn_times() override { return 2; }

    const unsigned * non_zeroes() const
    {
      return nonzero_.get();
    }

  public:
    std::unique_ptr<unsigned[]> nonzero_;
  };


  class TakeAndAdd : public Algorithm {
  public:

    struct InputEl {
      unsigned list;
      unsigned el;
    };

    Ret init(const unsigned npars[], unsigned npars_size,
             std::vector<std::vector<InputEl>> && elems);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

  private:
    std::vector<std::vector<InputEl>> elems_;
  };


  class TakeAndAddBL : public Algorithm {
  public:

    struct InputEl {
      unsigned list1;
      unsigned el1;
      unsigned list2;
      unsigned el2;
    };

    Ret init(const unsigned npars[], unsigned npars_size,
             std::vector<std::vector<InputEl>> && elems);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

  private:
    std::vector<std::vector<InputEl>> elems_;
  };


  class SparseMatrixMul : public Algorithm {
  public:

    // Call this after setting 'size' and 'cols' for 'row1' and 'row2'
    Ret init(unsigned nrows1, unsigned ncols1, unsigned ncols2);

    virtual Ret evaluate(Context * ctxt,
                         AlgInput xin[], Mod mod, AlgorithmData * data,
                         UInt xout[]) const override;

    struct RowInfo {
    private:
      std::size_t start = 0;
    public:
      std::size_t size = 0;
      std::unique_ptr<unsigned[]> cols;
    private:
      friend class SparseMatrixMul;
    };

    std::vector<RowInfo> row1, row2;

  private:
    unsigned nr1_, nc1_, nc2_;
  };

} // namespace fflow


#endif // FFLOW_ALG_LISTS_HH
