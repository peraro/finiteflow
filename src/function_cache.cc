#include <fflow/function_cache.hh>
#include "spooky_v2.hh"

namespace fflow {

  namespace  {

    unsigned next_pow_of_2(std::size_t x)
    {
      unsigned i=1;
      while (unsigned(1)<<i < x)
        ++i;
      return i;
    }

  } // namespace


  void UIntAlloc::reserve_more(std::size_t n)
  {
    free_mem_start_ = new UInt[n];
    free_mem_end_ = free_mem_start_ + n;
    chunks_.push_back(std::unique_ptr<UInt[]>(free_mem_start_));
  }

  void UIntAlloc::merge(UIntAlloc && oth)
  {
    chunks_.reserve(chunks_.size() + oth.chunks_.size());
    for (auto & c : oth.chunks_)
      chunks_.push_back(std::move(c));
    free_mem_start_ = oth.free_mem_start_;
    free_mem_end_ = oth.free_mem_end_;
    oth.clear();
  }

  constexpr double UIntCache::load_factor_;

  void UIntCache::reserve_more(std::size_t n, std::size_t allocate)
  {
    std::size_t newsize = size_ + n;
    std::size_t cap = get_cap_();

    if (log2cap_ == EMPTY_)
      init_first_(next_pow_of_2(2*(newsize+1)));
    else if (double(newsize)/double(cap) > load_factor_)
      rehash_(next_pow_of_2(2*(newsize+1)));

    if (allocate)
      alloc_.reserve_more(allocate);
  }

  UInt UIntCache::get_hash_(const UInt * x) const
  {
    return SpookyHash::Hash64((void*)(x), nparsin_*sizeof(UInt), OFFSET_5);
  }

  void UIntCache::init_first_(unsigned log2cap)
  {
    log2cap_ = log2cap;
    std::size_t cap = get_cap_();
    mask_ = cap - 1;
    data_.resize(cap);
  }

  bool UIntCache::is_in_pos_(const UInt * x, UInt hash, unsigned pos) const
  {
    const Entry & entry = data_[pos];
    if (entry.hashv == hash && std::equal(x, x+nparsin_, entry.key))
      return true;
    return false;
  }

  bool UIntCache::find(const UInt * xin, UInt ** xout) const
  {
    if (log2cap_ == EMPTY_)
      return false;

    UInt hash = get_hash_(xin);
    std::size_t pos = hash & mask_;

    bool found = false;
    while (1) {
      if (data_[pos].empty()) {
        break;
      } else if (is_in_pos_(xin, hash, pos)) {
        found = true;
        *xout = data_[pos].val;
        break;
      }
      pos = (pos + 1) & mask_;
    };

    return found;
  }


  // looks for an entry and returns the registered regxin and a new
  // xout to fill in
  bool UIntCache::set_entry(const UInt * xin, unsigned nparsout,
                            UInt ** reg_xin, UInt ** new_xout)
  {
    if (log2cap_ == EMPTY_)
      return false;

    UInt hash = get_hash_(xin);
    std::size_t pos = hash & mask_;

    bool found = false;
    while (1) {
      if (data_[pos].empty()) {
        break;
      } else if (is_in_pos_(xin, hash, pos)) {
        found = true;
        *reg_xin = data_[pos].key;
        *new_xout = data_[pos].val = get_new_outptr(nparsout);
        break;
      }
      pos = (pos + 1) & mask_;
    };

    if (!found) {

      UInt * nxin = get_new_inptr();
      std::copy(xin, xin+nparsin_, nxin);
      *reg_xin = nxin;
      *new_xout = new_entry(nxin, nparsout);

    }

    return found;
  }


  UInt * UIntCache::new_entry(UInt * xin, unsigned nparsout)
  {
    if (log2cap_ == EMPTY_)
      init_first_();

    UInt hash = get_hash_(xin);
    std::size_t pos = hash & mask_;
    std::size_t cap = get_cap_();

    if (double(size_+1)/double(cap) > load_factor_) {
      rehash_(log2cap_+1);
      pos = hash & mask_;
      cap = get_cap_();
    }

    while (!data_[pos].empty())
      pos = (pos + 1) & mask_;

    data_[pos].key = xin;
    data_[pos].val = get_new_outptr(nparsout);
    data_[pos].hashv = hash;
    ++size_;

    return data_[pos].val;
  }


  void UIntCache::rehash_(unsigned log2cap)
  {
    std::size_t newcap = get_cap_from_(log2cap);
    std::vector<Entry> newdata(newcap);
    std::size_t newmask = newcap - 1;

    for (auto & entry : data_) {
      std::size_t newpos = entry.hashv & newmask;
      while (1) {
        if (newdata[newpos].empty()) {
          newdata[newpos] = entry;
          break;
        }
        newpos = (newpos + 1) & newmask;
      };
    }

    data_.swap(newdata);
    log2cap_ = log2cap;
    mask_ = newmask;
  }


  void UIntCache::merge_caches(std::unique_ptr<UIntCache> caches[],
                               std::size_t n)
  {
    std::size_t additional_size = 0;
    std::size_t additional_chunks = 0;
    for (unsigned i=0; i<n; ++i) {
      additional_size += (*caches[i]).size_;
      additional_chunks += (*caches[i]).alloc_.n_chuncks();
    }
    reserve_more(additional_size, 0);

    alloc_.reserve_more_chunks(additional_chunks);

    for (unsigned i=0; i<n; ++i) {

      UIntCache & c = *caches[i];
      unsigned cap = c.get_cap_();

      for (unsigned i=0; i<cap; ++i) {
        Entry & entry = c.data_[i];
        if (!entry.empty()) {
          std::size_t pos = entry.hashv & mask_;
          while (1) {
            if (data_[pos].empty()) {
              data_[pos] = entry;
              break;
            }
            pos = (pos + 1) & mask_;
          }
          ++size_;
        }
      }

      alloc_.merge(std::move(c.alloc_));
      c.clear();
    }

  }


  void UIntCache::set_first_value_entries(UInt val)
  {
    for (auto entry : data_)
      if (!entry.empty())
        entry.val[0] = val;
  }

} // namespace fflow
