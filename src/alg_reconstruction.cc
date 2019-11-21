#include <algorithm>
#include <fstream>
#include <fflow/alg_reconstruction.hh>

namespace fflow {

  void copy_input_params(const UInt x[], unsigned nv, Mod mod,
                         std::unique_ptr<UInt[]> & xptr)
  {
    if (xptr.get() == nullptr)
      xptr.reset(new UInt[nv+1]);
    std::copy(x, x+nv, xptr.get());
    xptr[nv] = mod.n();
  }

  void copy_params(const UInt x[], unsigned nv,
                   std::unique_ptr<UInt[]> & xptr)
  {
    if (x == nullptr)
      xptr.reset(nullptr);
    if (xptr.get() == nullptr)
      xptr.reset(new UInt[nv]);
    std::copy(x, x+nv, xptr.get());
  }


  UInt AlgorithmRatFun::evaluate(const UInt x[], Mod mod)
  {
    copy_input_params(x, alg_->nparsin[0], mod, xin_);
    auto lookup = cache_.find(xin_);

    if (lookup == cache_.end()) {

      if (xout_.get() == nullptr)
        xout_.reset(new UInt[alg_->nparsout]);
      Ret ret = alg_->evaluate(ctxt_, &x, mod, alg_data_, xout_.get());
      if (ret == FAILED) {
        cache_[std::move(xin_)] = std::unique_ptr<UInt[]>(nullptr);
        return FAILED;
      }

      UInt res = xout_[idx];
      cache_[std::move(xin_)] = std::move(xout_);
      return res;

    }

    if ((*lookup).second.get() == nullptr)
      return FAILED;

    return (*lookup).second[idx];
  }


  UInt SampledRatFun::evaluate(const UInt x[], Mod mod)
  {
    const unsigned nv = nvars();

    copy_input_params(x, nv, mod, x_);
    UInt * xout;
    bool lookup = samples_->find(x_.get(), &xout);

    if (!lookup) {
      register_missing(mod.n());
      return FAILED;
    }

    if (xout[0] == FAILED)
      return FAILED;

    return xout[idx];
  }


  void SampledRatFun::register_missing(UInt mod)
  {
    auto lookup = std::find(missing_mods_.begin(), missing_mods_.end(), mod);
    if (lookup == missing_mods_.end())
      missing_mods_.push_back(mod);
  }


  UInt SparselySampledRatFun::evaluate(const UInt x[], Mod mod)
  {
    const unsigned nv = nvars();

    copy_input_params(x, nv, mod, x_);
    UInt * xout;
    bool lookup = samples_->find(x_.get(), &xout);

    if (!lookup) {
      register_missing(mod.n());
      return FAILED;
    }

    if (xout[flags_size_] == FAILED)
      return FAILED;

    unsigned out_idx=0;
    ConstBitArrayView av(xout);
    for (unsigned j=0; j<idx; ++j)
      if (av.get(j))
        ++out_idx;

    return xout[flags_size_+out_idx];
  }


  void SparselySampledRatFun::register_missing(UInt mod)
  {
    auto lookup = std::find(missing_mods_.begin(), missing_mods_.end(), mod);
    if (lookup == missing_mods_.end())
      missing_mods_.push_back(mod);
  }


  UInt GenerateSamplePoints::evaluate(const UInt x[], Mod mod)
  {
    const unsigned nv = nvars();

    copy_input_params(x, nv, mod, x_);

    if (complement_) {
      UInt * xout;
      bool clookup = complement_->find(x_.get(), &xout);
      if (clookup)
        return 0;
    }

    auto lookup = set_.find(x_);
    if (lookup == set_.end()) {
      set_.insert(std::move(x_));
    }

    return 0;
  }


  void GenerateSamplePoints::append_to_vector(SamplePointsVector & vec,
                                              std::size_t extra_entries) const
  {
    const unsigned nv = nvars();
    const std::size_t size = set_.size();

    vec.reserve(vec.size() + size);

    for (const auto & x : set_) {
      std::unique_ptr<UInt[]> xv(new UInt[nv + extra_entries + 1]());
      std::copy(x.get(), x.get()+nv+1, xv.get());
      vec.push_back(std::move(xv));
    }
  }


  UInt VerifySamplePoints::evaluate(const UInt x[], Mod mod)
  {
    const unsigned nv = nvars();

    copy_input_params(x, nv, mod, x_);

    UInt pos = FAILED;
    {
      UInt * pos_ptr;
      bool clookup = positions_->find(x_.get(), &pos_ptr);
      if (clookup)
        pos = pos_ptr[0];
    }

    if (pos != FAILED) {
      BitArrayView barr(samples_[pos].get() + nvars_ + 1);
      barr.set(idx);
    }

    return 0;
  }


  Ret SamplePointsFromFile::load_samples(unsigned nparsin,
                                         unsigned nparsout,
                                         SamplePointsVector & samples)
  {
    return ::fflow::load_samples(file_,
                                 nparsin + bit_array_u64size(nparsout),
                                 samples, samples_start_, samples_size_);
  }

  Ret SamplePointsFromVector::load_samples(unsigned, unsigned,
                                           SamplePointsVector & samples)
  {
    samples = std::move(vec_);
    return SUCCESS;
  }


  void sort_by_mod(std::unique_ptr<UInt[]> * start,
                   std::unique_ptr<UInt[]> * end,
                   unsigned nvars)
  {
    std::sort(start, end,
              [nvars] (const std::unique_ptr<UInt[]> & a,
                       const std::unique_ptr<UInt[]> & b)
              {
                return a[nvars] > b[nvars];
              });
  }


  void evaluate_and_cache_samples(const Algorithm & alg,
                                  AlgorithmData * alg_data,
                                  Context * ctxt,
                                  const std::unique_ptr<UInt[]> * start,
                                  const std::unique_ptr<UInt[]> * end,
                                  UIntCache & cache)
  {
    Mod mod;
    const unsigned nparsin = alg.nparsin[0];
    const unsigned nparsout = alg.nparsout;

    for (; start != end; ++start) {

      UInt * xin = cache.get_new_inptr();
      std::copy((*start).get(), (*start).get()+nparsin+1, xin);

      UInt * xout = cache.new_entry(xin, nparsout);

      if (mod.n() != xin[nparsin])
        mod = Mod(xin[nparsin]);
      Ret ret = alg.evaluate(ctxt, &xin, mod, alg_data, xout);

      if (ret == FAILED)
        xout[0] = FAILED;
    }
  }


  void evaluate_and_sparse_cache_samples(const Algorithm & alg,
                                         AlgorithmData * alg_data,
                                         Context * ctxt,
                                         const std::unique_ptr<UInt[]> * start,
                                         const std::unique_ptr<UInt[]> * end,
                                         UIntCache & cache)
  {
    Mod mod;
    const unsigned nparsin = alg.nparsin[0];
    const unsigned nparsout = alg.nparsout;
    unsigned flags_size = bit_array_u64size(nparsout);
    std::unique_ptr<UInt[]> allxout(new UInt[nparsout]);

    for (; start != end; ++start) {

      ConstBitArrayView needed((*start).get()+nparsin+1);

      UInt * xin = nullptr;
      UInt * xout = nullptr;
      unsigned this_nout = bit_array_nonzeroes(needed, nparsout);

      xin = cache.get_new_inptr();
      std::copy((*start).get(), (*start).get()+nparsin+1, xin);
      xout = cache.new_entry(xin, this_nout+flags_size);

      if (mod.n() != xin[nparsin])
        mod = Mod(xin[nparsin]);
      Ret ret = alg.evaluate(ctxt, &xin, mod, alg_data, allxout.get());

      if (ret == FAILED) {
        xout[flags_size] = FAILED;
      } else {
        std::copy((*start).get()+nparsin+1,
                  (*start).get()+nparsin+1+flags_size, xout);
        unsigned xidx = flags_size;
        for (unsigned j=0; j<nparsout; ++j)
          if (needed.get(j))
            xout[xidx++] = allxout[j];
      }
    }
  }


  Ret dump_samples(const char * filename,
                   unsigned npars_in,
                   const SamplePointsVector & vec)
  {
    std::ofstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto dump = [&file](UInt z)
      {
        file.write(reinterpret_cast<char*>(&z), sizeof(UInt));
      };

	dump(npars_in);
    dump(vec.size());

    for (const auto & xi : vec) {
      for (unsigned i=0; i<npars_in+1; ++i)
        dump(xi[i]);
    }

    return SUCCESS;
  }

  UInt samples_file_size(const char * filename,
                         unsigned npars_in)
  {
    std::ifstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto load = [&file](UInt & z)
      {
        file.read(reinterpret_cast<char*>(&z), sizeof(UInt));
      };

    UInt check_npars_in, size;

    load(check_npars_in);
    if (file.fail() ||
        (npars_in != (~unsigned(0)) && check_npars_in != npars_in))
      return FAILED;

    load(size);
    if (file.fail())
      return FAILED;

    return size;
  }

  Ret load_samples(const char * filename,
                   unsigned npars_in,
                   SamplePointsVector & vec,
                   std::size_t start,
                   std::size_t n_samples)
  {
    std::ifstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto load = [&file](UInt & z)
      {
        file.read(reinterpret_cast<char*>(&z), sizeof(UInt));
      };

    UInt check_npars_in, size;

    load(check_npars_in);
    if (file.fail() || check_npars_in != npars_in)
      return FAILED;

    load(size);
    if (file.fail())
      return FAILED;
    if (!n_samples)
      n_samples = size;
    else if (start + n_samples > size)
      return FAILED;
    vec.reserve(vec.size() + n_samples);

    if (start)
      file.ignore((npars_in+1)*sizeof(UInt)*start);

    for (unsigned i=0; i<n_samples; ++i) {
      std::unique_ptr<UInt[]> x(new UInt[npars_in + 1]);
      for (unsigned j=0; j<npars_in+1; ++j) {
        load(x[j]);
        if (file.fail())
          return FAILED;
      }
      vec.push_back(std::move(x));
    }

    return SUCCESS;
  }

  Ret dump_evaluations(const char * filename,
                       unsigned npars_in,
                       unsigned npars_out,
                       const UIntCache & cache)
  {
    std::ofstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto dump = [&file](UInt z)
      {
        file.write(reinterpret_cast<char*>(&z), sizeof(UInt));
      };

	dump(npars_in);
    dump(npars_out);
    dump(cache.size());
    dump(cache.allocated_uint64());

    for (const auto xi : cache) {
      for (unsigned i=0; i<npars_in+1; ++i)
        dump(xi.first[i]);
      if (xi.second[0] == FAILED) {
        for (unsigned i=0; i<npars_out; ++i)
          dump(FAILED);
      } else {
        for (unsigned i=0; i<npars_out; ++i)
          dump(xi.second[i]);
      }
    }

    return SUCCESS;
  }

  Ret load_evaluations(const char * filename,
                       unsigned npars_in,
                       unsigned npars_out,
                       UIntCache & cache)
  {
    std::ifstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto load = [&file](UInt & z)
      {
        file.read(reinterpret_cast<char*>(&z), sizeof(UInt));
      };

    UInt check_npars_in, check_npars_out, size, alloc_u64;

    load(check_npars_in);
    if (file.fail() || check_npars_in != npars_in)
      return FAILED;

    load(check_npars_out);
    if (file.fail() || check_npars_out != npars_out)
      return FAILED;

    load(size);
    if (file.fail())
      return FAILED;
    load(alloc_u64);
    if (file.fail())
      return FAILED;
    cache.reserve_more(size, alloc_u64);

    for (unsigned i=0; i<size; ++i) {
      UInt * xin = cache.get_new_inptr();
      for (unsigned j=0; j<npars_in+1; ++j) {
        load(xin[j]);
        if (file.fail())
          return FAILED;
      }
      UInt * xout = cache.new_entry(xin, npars_out);
      for (unsigned j=0; j<npars_out; ++j) {
        load(xout[j]);
        if (file.fail())
          return FAILED;
      }
    }

    return SUCCESS;
  }

  Ret load_evaluations_for_var(const char * filename,
                               unsigned npars_in,
                               unsigned npars_out,
                               UIntCache & cache,
                               unsigned var)
  {
    std::ifstream file;
    file.open(filename, std::ios::binary);
    if (file.fail())
      return FAILED;

    auto load = [&file](UInt & z)
      {
        file.read(reinterpret_cast<char*>(&z), sizeof(UInt));
      };

    UInt check_npars_in, check_npars_out, size, alloc_u64;

    load(check_npars_in);
    if (file.fail() || check_npars_in != npars_in)
      return FAILED;

    load(check_npars_out);
    if (file.fail() || check_npars_out != npars_out)
      return FAILED;

    load(size);
    if (file.fail())
      return FAILED;
    load(alloc_u64);
    if (file.fail())
      return FAILED;
    cache.reserve_more(size,size*(npars_in+1));

    for (unsigned i=0; i<size; ++i) {
      UInt * xin = cache.get_new_inptr();
      for (unsigned j=0; j<npars_in+1; ++j) {
        load(xin[j]);
        if (file.fail())
          return FAILED;
      }
      UInt * xout = cache.new_entry(xin, npars_out);
      if (var)
        file.ignore(var*sizeof(UInt));
      if (file.fail())
        return FAILED;
      load(xout[0]);
      if (file.fail())
        return FAILED;
      file.ignore((npars_out-var-1)*sizeof(UInt));
      if (file.fail())
        return FAILED;
    }

    return SUCCESS;
  }

} // namespace fflow
