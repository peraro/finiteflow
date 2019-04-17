#ifndef FFLOW_FUNCTION_CACHE_HH
#define FFLOW_FUNCTION_CACHE_HH

#include <unordered_map>
#include <memory>
#include <fflow/common.hh>

namespace fflow {

  template <std::size_t size> struct IntArrayXiEqual;
  template <std::size_t size> struct IntArrayXiHash;

  template <std::size_t size>
  struct IntArrayXiEqual {

    static const std::size_t n_ = size;

    typedef std::array<UInt,n_> Key;

    bool operator()(const Key & lhs, const Key & rhs) const
    {
      for (unsigned i=0; i<n_; ++i)
        if (lhs[i] != rhs[i])
          return false;
      return true;
    }
  };

  template <std::size_t size>
  struct IntArrayXiHash {

    static const std::size_t n_ = size;

    typedef std::array<UInt,n_> Key;

    std::size_t operator()(const Key & k) const
    {
      UInt h = k[0];
      for (unsigned i=1; i<n_; ++i)
        h = k[i] + 0x9e3779b9 + (h << 6) + (h >> 2);
      return std::size_t(h);
    }
  };


  // size=0 --> dynamic case

  template<>
  struct IntArrayXiEqual<0> {

    typedef std::unique_ptr<UInt[]> Key;

    bool operator()(const Key & lhs, const Key & rhs) const
    {
      for (unsigned i=0; i<n_; ++i)
        if (lhs[i] != rhs[i])
          return false;
      return true;
    }

    std::size_t n_;
  };

  template<>
  struct IntArrayXiHash<0> {

    typedef std::unique_ptr<UInt[]> Key;

    std::size_t operator()(const Key & k) const
    {
      UInt h = k[0];
      for (unsigned i=1; i<n_; ++i)
        h = k[i] + 0x9e3779b9 + (h << 6) + (h >> 2);
      return std::size_t(h);
    }

    std::size_t n_;
  };

  template<std::size_t size>
  using FunctionCache = std::unordered_map< typename IntArrayXiHash<size>::Key,
                                            std::unique_ptr<UInt[]>,
                                            IntArrayXiHash<size>,
                                            IntArrayXiEqual<size> >;


  // -- public domain code from: lookup3.c, by Bob Jenkins, May 2006 --

#define FFLOW_HASH_ROT(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

#define FFLOW_HASH_FINAL_MIX(a,b,c)             \
  {                                             \
    c ^= b; c -= FFLOW_HASH_ROT(b,14);          \
    a ^= c; a -= FFLOW_HASH_ROT(c,11);          \
    b ^= a; b -= FFLOW_HASH_ROT(a,25);          \
    c ^= b; c -= FFLOW_HASH_ROT(b,16);          \
    a ^= c; a -= FFLOW_HASH_ROT(c,4);           \
    b ^= a; b -= FFLOW_HASH_ROT(a,14);          \
    c ^= b; c -= FFLOW_HASH_ROT(b,24);          \
  }

  // -- end of public domain code --

  inline UInt hash_uint(UInt key)
  {
    std::uint32_t a = key & 0xffffffff, b = key >> 32, c = 0;
    FFLOW_HASH_FINAL_MIX(a, b, c);
    return (UInt(c) << 32) + b;
  }

  // some offsets for internal use:
  //
  //   OFFSET == 1, i == 0 : default t0
  //   OFFSET_1 : for sampling in z (the main variables)
  //   OFFSET_2 : for sampling in t (additional parameters)
  //   OFFSET_3 : for auxiliary variables
  //   OFFSET_4 : for generating shifts in variables
  //   OFFSET_5 : used as default seed for hash functions
  const UInt OFFSET_MAX_ = ~UInt(0) - UInt(7);
  const UInt OFFSET_0 = 3;
  const UInt OFFSET_1 = OFFSET_MAX_/UInt(10);
  const UInt OFFSET_2 = UInt(2)*(OFFSET_MAX_/UInt(10));
  const UInt OFFSET_3 = UInt(3)*(OFFSET_MAX_/UInt(10));
  const UInt OFFSET_4 = UInt(4)*(OFFSET_MAX_/UInt(10));
  const UInt OFFSET_5 = UInt(5)*(OFFSET_MAX_/UInt(10));

  // sample_uint(1, 0, Mod(smallest_prime))
  const UInt OFFSET_0_HASH = 4680351337364691545ULL;

  // We use this to generate pseudo-random integers with
  //
  //    sample_uint(offset, i, mod)
  //
  // with offset > 0.  Note that we don't have any statistical
  // requirement, and this is just used to generated reasonably
  // independent sample points in a deterministic way.
  inline UInt sample_uint(UInt offset, std::uint32_t i, Mod mod)
  {
    return red_mod(hash_uint(offset+UInt(i)), mod);
  }



  class UIntAlloc {
  public:

    UIntAlloc() = default;

    UInt * alloc(std::size_t n)
    {
      if (free_mem_() < n)
        reserve_more(DEFAULT_CHUNK_SIZE*(n/DEFAULT_CHUNK_SIZE+1));
      UInt * ret = free_mem_start_;
      free_mem_start_ += n;
      return ret;
    }

    void reserve_more(std::size_t n);

    void reserve_more_chunks(std::size_t n)
    {
      chunks_.reserve(chunks_.size() + n);
    }

    std::size_t n_chuncks() const
    {
      return chunks_.size();
    }

    void clear()
    {
      chunks_.clear();
      free_mem_start_ = free_mem_end_ = nullptr;
    }

    void merge(UIntAlloc && oth);

  private:

    std::size_t free_mem_() const
    {
      return free_mem_end_ - free_mem_start_;
    }

    const std::size_t DEFAULT_CHUNK_SIZE = 4096/sizeof(UInt);
    // Data for memory allocation
    std::vector<std::unique_ptr<UInt[]>> chunks_;
    UInt *free_mem_start_ =  nullptr, *free_mem_end_ = nullptr;
  };


  class UIntCacheIter;

  class UIntCache {
  public:

    UIntCache() = default;

    unsigned nparsin() const
    {
      return nparsin_;
    }

    void init(unsigned nparsin)
    {
      nparsin_ = nparsin;
    }

    UInt * get_new_inptr()
    {
      alloc_uint64_  += nparsin_;
      return alloc_.alloc(nparsin_);
    }

    UInt * get_new_outptr(unsigned nparsout)
    {
      alloc_uint64_ += nparsout;
      return alloc_.alloc(nparsout);
    }

    std::size_t size() const
    {
      return size_;
    }

    std::size_t allocated_uint64() const
    {
      return alloc_uint64_;
    }

    // creates a new entry and returns the val pointer, assuming it is
    // not already there
    UInt * new_entry(UInt * xin, unsigned nparsout);

    // if found *xout will point to val
    bool find(const UInt * xin, UInt ** xout) const;

    // looks for an entry and returns the registered regxin and a new
    // xout to fill in
    bool set_entry(const UInt * xin, unsigned nparsout,
                   UInt ** reg_xin, UInt ** new_xout);

    void reserve_more(std::size_t n, std::size_t allocate);

    // this will clear the other caches
    void merge_caches(std::unique_ptr<UIntCache> caches[], std::size_t n);

    void clear()
    {
      data_.clear();
      mask_ = 0;
      size_ = 0;
      log2cap_ = EMPTY_;
      alloc_.clear();
    }

    UIntCacheIter begin() const;
    UIntCacheIter end() const;

    // sets the first value of all entries to "val"
    void set_first_value_entries(UInt val);

  private:

    struct Entry {
      std::size_t hashv = 0;
      UInt * key = nullptr;
      UInt * val = nullptr;

      bool empty() const
      {
        return key == nullptr;
      }
    };

  private:

    friend class UIntCacheIter;

    void init_first_(unsigned log2cap = DEFAULT_LOG2_CAP_);
    UInt get_hash_(const UInt * x) const;
    void rehash_(unsigned log2cap);
    bool is_in_pos_(const UInt * x, UInt hash, unsigned pos) const;

    static std::size_t get_cap_from_(unsigned log2cap)
    {
      return std::size_t(1) << log2cap;
    }

    std::size_t get_cap_() const
    {
      return get_cap_from_(log2cap_);
    }

  private:

    static constexpr double load_factor_ = 0.5;
    static const unsigned EMPTY_ = ~unsigned(0);
    static const unsigned DEFAULT_LOG2_CAP_ = 7;

    std::vector<Entry> data_;
    std::size_t mask_ = 0;
    std::size_t size_ = 0;
    std::size_t alloc_uint64_ = 0;
    unsigned nparsin_ = 0;
    unsigned log2cap_ = EMPTY_;
    UIntAlloc alloc_;
  };


  class UIntCacheIter {
  public:

    UIntCacheIter() = default;
    UIntCacheIter(const UIntCacheIter & oth) = default;

    UIntCacheIter & operator ++()
    {
      const auto & data = cache_->data_;
      do {
        ++i_;
      } while (i_ < data.size() && data[i_].empty());
      return *this;
    }

    bool operator ==(const UIntCacheIter & oth) const
    {
      return i_ == oth.i_;
    }

    bool operator !=(const UIntCacheIter & oth) const
    {
      return i_ != oth.i_;
    }

    std::pair<const UInt *, const UInt *> operator* () const
    {
      const auto & data = cache_->data_;
      return std::make_pair(data[i_].key, data[i_].val);
    }

  private:
    friend class UIntCache;

    UIntCacheIter(const UIntCache * cache, std::size_t i)
      : cache_(cache), i_(i) {}

  private:
    const UIntCache * cache_ = nullptr;
    std::size_t i_ = 0;
  };

  inline UIntCacheIter UIntCache::begin() const
  {
    std::size_t first=0;
    while (first<data_.size() && data_[first].empty())
      ++first;
    return UIntCacheIter(this, first);
  }

  inline UIntCacheIter UIntCache::end() const
  {
    return UIntCacheIter(this, data_.size());
  }


} // namespace fflow

#endif // FFLOW_FUNCTION_CACHE_HH
