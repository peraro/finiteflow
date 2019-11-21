#ifndef FFLOW_ALG_RECONSTRUCTION_HH
#define FFLOW_ALG_RECONSTRUCTION_HH

#include <unordered_set>
#include <fflow/algorithm.hh>
#include <fflow/multivariate_reconstruction.hh>
#include <fflow/function_cache.hh>
#include <fflow/integer_math.hh>

namespace fflow {

  // assumes xptr is either nullptr or an array with the rigth size = nv+1
  void copy_input_params(const UInt x[], Mod mod, unsigned nv,
                         std::unique_ptr<UInt[]> & xptr);


  class AlgorithmRatFun : public RatFun {
  public:

    virtual UInt evaluate(const UInt x[], Mod mod) override;

    void set_algorithm(const Algorithm & alg)
    {
      unsigned npars_in = alg.nparsin[0];
      unsigned npars_out = alg.nparsout;
      xin_.reset(new UInt[npars_in+1]);
      xout_.reset(new UInt[npars_out]);
      cache_ = FunctionCache<0>(binomial(npars_in+5,5),
                                IntArrayXiHash<0>{npars_in+1},
                                IntArrayXiEqual<0>{npars_in+1});
      alg_ = &alg;
    }

    void set_algorithm_data(AlgorithmData * alg_data)
    {
      alg_data_ = alg_data;
    }

    void set_context(Context * ctxt)
    {
      ctxt_ = ctxt;
    }

    const Algorithm * get_alg() const
    {
      return alg_;
    }

  private:
    FunctionCache<0> cache_;
    std::unique_ptr<UInt[]> xin_;
    std::unique_ptr<UInt[]> xout_;
    const Algorithm * alg_ = nullptr;
    AlgorithmData * alg_data_ = nullptr;
    Context * ctxt_ = nullptr;

  public:
    unsigned idx = 0;
  };


  class SampledRatFun : public RatFun {
  public:

    SampledRatFun() : samples_(nullptr), missing_mods_(), idx(0) {}

    void set_function_cache(const UIntCache * samples)
    {
      samples_ = samples;
    }

    virtual UInt evaluate(const UInt x[], Mod mod) override;

    unsigned nvars() const
    {
      return samples_->nparsin() - 1;
    }

    void register_missing(UInt mod);

    const std::vector<UInt> & missing_mods()
    {
      return missing_mods_;
    }

  private:
    const UIntCache * samples_;
    std::vector<UInt> missing_mods_;
    std::unique_ptr<UInt[]> x_;

  public:
    unsigned idx;
  };


  class SparselySampledRatFun : public RatFun {
  public:

    SparselySampledRatFun()
      : samples_(nullptr), missing_mods_(),
        nparsin_(0), flags_size_(0), idx(0) {}

    void set_function_cache(const UIntCache * samples,
                            unsigned nparsin, unsigned nparsout)
    {
      samples_ = samples;
      nparsin_ = nparsin;
      flags_size_ = bit_array_u64size(nparsout);
    }

    virtual UInt evaluate(const UInt x[], Mod mod) override;

    unsigned nvars() const
    {
      return samples_->nparsin() - 1;
    }

    void register_missing(UInt mod);

    const std::vector<UInt> & missing_mods()
    {
      return missing_mods_;
    }

  private:
    const UIntCache * samples_;
    std::vector<UInt> missing_mods_;
    std::unique_ptr<UInt[]> x_;
    unsigned nparsin_, flags_size_;

  public:
    unsigned idx;
  };


  typedef std::vector<std::unique_ptr<UInt[]>> SamplePointsVector;
  typedef std::unique_ptr<UIntCache> CacheUPtr;


  class GenerateSamplePoints : public RatFun {
  public:
    GenerateSamplePoints(unsigned nvars)
      : set_(binomial(nvars+6,6),
             IntArrayXiHash<0>{nvars+1},
             IntArrayXiEqual<0>{nvars+1}),
        complement_(nullptr),
        x_(nullptr) {}

    virtual UInt evaluate(const UInt x[], Mod mod) override;

    unsigned nvars() const
    {
      return set_.hash_function().n_ - 1;
    }

    // If this is set, any sample point already in the cache is not
    // recorded
    void set_complement(const UIntCache * cache)
    {
      complement_ = cache;
    }

    void append_to_vector(SamplePointsVector & vec,
                          std::size_t extra_entries = 0) const;

  private:
    typedef std::unordered_set< IntArrayXiHash<0>::Key,
                                IntArrayXiHash<0>,
                                IntArrayXiEqual<0> > Set;

  private:
    Set set_;
    const UIntCache * complement_;
    std::unique_ptr<UInt[]> x_;
  };


  class VerifySamplePoints : public RatFun {
  public:
    VerifySamplePoints(unsigned nvars)
      : samples_(), positions_(), x_(nullptr), nvars_(nvars), idx(0) {}

    virtual UInt evaluate(const UInt x[], Mod mod) override;

    unsigned nvars() const
    {
      return nvars_;
    }

    void set_samples(std::unique_ptr<UInt[]> * samples,
                     const UIntCache * positions)
    {
      samples_ = samples;
      positions_ = positions;
    }

  private:
    typedef std::unordered_set< IntArrayXiHash<0>::Key,
                                IntArrayXiHash<0>,
                                IntArrayXiEqual<0> > Set;

  private:
    std::unique_ptr<UInt[]> * samples_;
    const UIntCache * positions_;
    std::unique_ptr<UInt[]> x_;
    unsigned nvars_;

  public:
    unsigned idx;
  };


  class SamplePointsGenerator {
  public:
    virtual Ret load_samples(unsigned nparsin, unsigned nparsout,
                             SamplePointsVector & samples) = 0;
    virtual ~SamplePointsGenerator() {}
  };


  class SamplePointsFromFile : public SamplePointsGenerator {
  public:

    explicit SamplePointsFromFile(const char * file,
                                  std::size_t samples_start = 0,
                                  std::size_t samples_size = 0)
      : file_(file),
        samples_start_(samples_start),
        samples_size_(samples_size) {}

    virtual Ret load_samples(unsigned nparsin, unsigned nparsout,
                             SamplePointsVector & samples);

  private:
    const char * file_ = nullptr;
    std::size_t samples_start_ = 0;
    std::size_t samples_size_ = 0;
  };


  class SamplePointsFromVector : public SamplePointsGenerator {
  public:
    void setVector(SamplePointsVector && samples)
    {
      vec_ = std::move(samples);
    }

    virtual Ret load_samples(unsigned nparsin, unsigned nparsout,
                             SamplePointsVector & samples);

  private:
    SamplePointsVector vec_;
  };


  void sort_by_mod(std::unique_ptr<UInt[]> * start,
                   std::unique_ptr<UInt[]> * end,
                   unsigned nvars);


  // for best performance, use on distinct sample points sorted by mod
  void evaluate_and_cache_samples(const Algorithm & alg,
                                  AlgorithmData * alg_data,
                                  Context * ctxt,
                                  const std::unique_ptr<UInt[]> * start,
                                  const std::unique_ptr<UInt[]> * end,
                                  UIntCache & cache);

  void evaluate_and_sparse_cache_samples(const Algorithm & alg,
                                         AlgorithmData * alg_data,
                                         Context * ctxt,
                                         const std::unique_ptr<UInt[]> * start,
                                         const std::unique_ptr<UInt[]> * end,
                                         UIntCache & cache);

  inline void merge_function_caches(CacheUPtr cachein[],
                                    unsigned ncaches,
                                    UIntCache & cacheout)
  {
    cacheout.merge_caches(cachein, ncaches);
  }

  Ret dump_samples(const char * filename,
                   unsigned npars_in,
                   const SamplePointsVector & vec);

  UInt samples_file_size(const char * filename,
                         unsigned npars_in = (~unsigned(0)));

  Ret load_samples(const char * filename,
                   unsigned npars_in,
                   SamplePointsVector & vec,
                   std::size_t start=0,
                   std::size_t n_samples=0);

  Ret dump_evaluations(const char * filename,
                       unsigned npars_in,
                       unsigned npars_out,
                       const UIntCache & cache);

  Ret load_evaluations(const char * filename,
                       unsigned npars_in,
                       unsigned npars_out,
                       UIntCache & cache);

  Ret load_evaluations_for_var(const char * filename,
                               unsigned npars_in,
                               unsigned npars_out,
                               UIntCache & cache,
                               unsigned var);

} // namespace fflow


#endif // FFLOW_ALG_RECONSTRUCTION_HH

