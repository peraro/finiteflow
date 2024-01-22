#ifndef FFLOW_COMMON_HH
#define FFLOW_COMMON_HH

#include <cstdint>
#include <ostream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <flint/flint.h>
#if defined(__FLINT_VERSION) && (__FLINT_VERSION >= 3)
# include <flint/nmod.h>
#endif
#include <flint/nmod_vec.h>
#include <fflow/config.hh>
#include <fflow/format.h>
#include <fflow/ostream.h>
#ifndef FFLOW_NO_DBG
# include <fflow/debug.hh>
#endif

#ifdef FFLOW_USE_EXCEPTIONS
# define FFLOW_THROW(exc) throw exc
#else
# define FFLOW_THROW(exc) std::terminate()
#endif

#define FFLOW_LOOP_UNROLL 1

namespace fflow {

  typedef std::int64_t Int;
  typedef std::uint64_t UInt;

  static_assert(sizeof(mp_limb_t) == sizeof(UInt),
                "Flint integer type is not a 64-bit integer");

  // Shoup’s modular multiplication algorithm:
  // http://web.maths.unsw.edu.au/~davidharvey/talks/fastntt-2-talk.pdf
  //
  // Given w in Z_p, with p<2^63, if we need to compute
  //
  //   t * w mod p
  //
  // for several values of t in Z_p, we can precompute
  //
  //   wp = quotient(2^64 * w, p)
  //
  // with precomp_mul_shoup, and then use the fast multiplication
  // algoritm mul_mod_shoup
  inline UInt precomp_mul_shoup_impl(UInt w, UInt p)
  {
    UInt q, r;
    udiv_qrnnd(q, r, w, 0, p);
    return q;
  }
  inline UInt mul_mod_shoup_impl(UInt t, UInt w, UInt wp, UInt p)
  {
    UInt q, r;

    // q = high bits of wp*t (r is ignored)
    umul_ppmm(q, r, wp, t);

    // note: the products below may overflow, but it's okay since we are
    // supposed to compute this modulo 2^64
    r = w*t - q*p;

    if (r >= p)
      r -= p;

    return r;
  }

  // returns 2^64 mod p, assuming 2 p < 2^64 < 3 p
  inline UInt beta_mod(UInt p)
  {
    return -(p<<1);
  }

  // This must be initialized with a prime "p" satisfying
  //
  //    2 p < 2^64 < 3 p
  //
  // These inequalities are assumed to be valid in the whole fflow
  // library.
  class Mod {
  public:
    Mod() : mod_{0,0,0} {}

    explicit Mod(UInt mod)
    {
      nmod_init(&mod_, mod);
    }

    UInt n() const
    {
      return mod_.n;
    }

    const nmod_t & nmod() const
    {
      return mod_;
    }

    void norm_to_one_()
    {
      mod_.norm = 1;
    }

  private:
    nmod_t mod_;
  };

  // this is used when the return value is just a flag
  typedef UInt Ret;

  // special constants
  const UInt SUCCESS = 0;
  const UInt FAILED = ~UInt(0);
  const UInt BAD_SHIFT = FAILED-1;
  const UInt MISSING_SAMPLES = FAILED-2;
  const UInt MISSING_PRIMES = FAILED-3;
  const UInt MUTABILITY_ERROR = FAILED-4;
  const UInt INVALID_ID = FAILED-5;

  template <typename T>
  inline void cswap(T & a, T & b)
  {
    T tmp = a;
    a = b;
    b = tmp;
  }

  inline UInt iabs(Int n)
  {
    return n < 0 ? -n : n;
  }

  inline Int sign(Int n)
  {
    return n < 0 ? -1 : 1;
  }

  // assumes a < 2 p
  inline UInt red_mod_(UInt a, Mod mod)
  {
    return a < mod.n() ? a : a - mod.n();
  }
  // faster a % mod
  inline UInt red_mod(UInt a, Mod mod)
  {
    UInt two_p = mod.n() << 1;
    return a >= two_p ? (a - two_p) : red_mod_(a, mod);
  }
  inline UInt add_mod(UInt a, UInt b, Mod mod)
  {
    UInt ret = a+b-mod.n();
    return Int(ret) < 0 ? ret + mod.n() : ret;
  }
  inline UInt sub_mod(UInt a, UInt b, Mod mod)
  {
    return _nmod_sub(a, b, nmod_t{mod.n(),0,0});
  }
  inline UInt neg_mod(UInt a, Mod mod)
  {
    return _nmod_sub(0, a, nmod_t{mod.n(),0,0});
  }
  // (ah*2^64 + al) % mod, assuming ah < mod
  inline UInt red_mod2(UInt ah, UInt al, Mod mod)
  {
    mod.norm_to_one_();
    const nmod_t & nmod = mod.nmod();
    NMOD_RED2(ah, ah, al, nmod);
    return ah;
  }
  // (ah*2^128 + am*2^64 + al) % mod, assuming ah < mod
  inline UInt red_mod3(UInt ah, UInt am, UInt al, Mod mod)
  {
    ah = red_mod2(ah, am, mod);
    return red_mod2(ah, al, mod);
  }
  // (ah*2^64 + al) % mod, no assumption on ah, al
  inline UInt red2_mod2(UInt ah, UInt al, Mod mod)
  {
    return red_mod2(red_mod(ah, mod), al, mod);
  }
  inline UInt mul_mod(UInt a, UInt b, Mod mod)
  {
    UInt resh, resl;
    umul_ppmm(resh, resl, a, b);
    return red_mod2(resh, resl, mod);
  }
  // a + b*c
  inline UInt apbc_mod(UInt a, UInt b, UInt c, Mod mod)
  {
    UInt prodh, prodl;
    umul_ppmm(prodh, prodl, b, c);
    add_ssaaaa(prodh, prodl, prodh, prodl, UInt(0), a);
    return red_mod2(prodh, prodl, mod);
  }
  // a - b*c
  inline UInt ambc_mod(UInt a, UInt b, UInt c, Mod mod)
  {
    b = neg_mod(b, mod);
    return apbc_mod(a, b, c, mod);
  }
  // a = a + b*c
  inline void addmul_mod(UInt & a, UInt b, UInt c, Mod mod)
  {
    a = apbc_mod(a, b, c, mod);
  }
  // a = a - b*c
  inline void submul_mod(UInt & a, UInt b, UInt c, Mod mod)
  {
    a = ambc_mod(a, b, c, mod);
  }
  // a*b - c*d
  inline UInt abmcd_mod_(UInt a, UInt b, UInt c, UInt d, Mod mod)
  {
    return sub_mod(mul_mod(a,b,mod), mul_mod(c,d,mod), mod);
  }
  inline UInt abmcd_mod(UInt a, UInt b, UInt c, UInt d, Mod mod)
  {
    UInt prodh, prodl, sumh=0, suml=0;
    UInt negc = neg_mod(c, mod);
    umul_ppmm(sumh, suml, a, b);
    umul_ppmm(prodh, prodl, negc, d);
    add_ssaaaa(sumh, suml, sumh, suml, prodh, prodl);
    return red_mod2(sumh, suml, mod);
  }
  inline UInt add3_mod(UInt a, UInt b, UInt c, Mod mod)
  {
    return add_mod(add_mod(a,b,mod),c,mod);
  }
  inline UInt mul3_mod(UInt a, UInt b, UInt c, Mod mod)
  {
    return mul_mod(mul_mod(a,b,mod),c,mod);
  }

  // fast calculation of 1/2 and -1/2 (note: mod is an odd prime)
  inline UInt half_mod(Mod mod)
  {
    return (mod.n()+1)/2;
  }
  inline UInt minus_half_mod(Mod mod)
  {
    return (mod.n()-1)/2;
  }

  // twice
  inline UInt twice_mod(UInt a, Mod mod)
  {
    return red_mod_(a<<1, mod);
  }


  inline UInt precomp_mul_shoup(UInt w, Mod mod)
  {
    return precomp_mul_shoup_impl(w, mod.n());
  }
  inline UInt mul_mod_shoup(UInt t, UInt w, UInt wp, Mod mod)
  {
    return mul_mod_shoup_impl(t, w, wp, mod.n());
  }
  inline void precomp_array_mul_shoup(const UInt w[], std::size_t n, Mod mod,
                                      UInt wp[])
  {
    for (unsigned i=0; i<n; ++i)
      wp[i] = precomp_mul_shoup(w[i], mod);
  }


  // assuming |n| < mod
  inline UInt iabs_mod(Int n, Mod mod)
  {
    return n < 0 ? neg_mod(-n,mod) : n;
  }

  // a minimalistic structure for storing rational numbers
  struct Rational {
    Int num, den;
  };

  inline bool operator == (const Rational & a, const Rational & b)
  {
    return (a.num == b.num) && (a.den == b.den);
  }

  inline std::ostream & operator << (std::ostream & os, const Rational & r)
  {
    return os << "(" << r.num << ")/(" << r.den << ")";
  }


  // basic compare
  template <typename T, typename F>
  int compare(const T & x, const T & y, F & isless)
  {
    return isless(y,x) - isless(x,y);
  }
  template <typename T>
  int compare(const T & x, const T & y)
  {
    return (y<x) - (x<y);
  }


  // Simple utility class for changing the value of a variable until
  // the end of the scope, and then restoring the old value (note:
  // must be copiable and assignable)
  template <typename T>
  class ScopeValue {
  public:
    ScopeValue(T & var, T value) : var_(var), old_val_(var)
    {
      var_ = value;
    }
    ~ScopeValue()
    {
      var_ = old_val_;
    }
  private:
    T & var_;
    T old_val_;
  };


  // bit masks
  template <std::size_t n> struct BitMask {
    static const UInt val = (UInt(1) << n) -1;
  };
  template <> struct BitMask<sizeof(UInt)> {
    static const UInt val = ~UInt(0);
  };


  template <typename T, std::size_t N_R, std::size_t N_C = N_R>
  struct StaticMatrixT {

    std::size_t rows() const
    {
      return N_R;
    }

    std::size_t columns() const
    {
      return N_C;
    }

    T & operator() (std::size_t i, std::size_t j)
    {
      return data_[i*N_C + j];
    }

    const T & operator() (std::size_t i, std::size_t j) const
    {
      return data_[i*N_C + j];
    }

    T data_[N_R*N_C];
  };

  template <typename T>
  struct DynamicMatrixT {

    void resize(std::size_t n)
    {
      resize(n,n);
    }

    void resize(std::size_t n, std::size_t m)
    {
      m_ = m;
      data_.resize(n*m);
    }

    void shrink_to_fit()
    {
      data_.shrink_to_fit();
    }

    void clear()
    {
      m_ = 0;
      data_.clear();
    }

    std::size_t rows() const
    {
      return m_ ? data_.size()/m_ : 0;
    }

    const T * row(std::size_t i) const
    {
      return data_.data() + i*m_;
    }

    T * row(std::size_t i)
    {
      return data_.data() + i*m_;
    }

    std::size_t columns() const
    {
      return m_;
    }

    T & operator() (std::size_t i, std::size_t j)
    {
      return data_[i*m_ + j];
    }

    const T & operator() (std::size_t i, std::size_t j) const
    {
      return data_[i*m_ + j];
    }

    std::vector<T> data_;
    std::size_t m_;
  };


  template <typename Iter>
  std::string join(Iter begin, Iter end, const std::string & separator)
  {
    MemoryWriter result;
    if (begin != end)
      result << *begin++;
    while (begin != end)
      result << separator << *begin++;
    return result.str();
  }


  // checked wrapper of malloc and realloc
  inline void * crealloc(void * ptr, std::size_t size)
  {
    ptr = std::realloc(ptr, size);
    if (!ptr)
      std::terminate();
    return ptr;
  }


  template <typename T>
  class MallocArray {
  public:
    MallocArray() : ptr_(nullptr) {}

    explicit MallocArray(std::size_t n)
      : ptr_((T*)crealloc(nullptr,n*sizeof(T))) {}

    MallocArray(const MallocArray & oth) = delete;

    MallocArray(MallocArray && oth) : ptr_(nullptr)
    {
      std::swap(ptr_, oth.ptr_);
    }

    ~MallocArray()
    {
      if (ptr_)
        std::free(ptr_);
    }

    MallocArray & operator=(const MallocArray & oth) = delete;

    MallocArray & operator=(MallocArray && oth)
    {
      MallocArray arr(std::move(oth));
      std::swap(*this,arr);
      return *this;
    }

    const T * get() const
    {
      return ptr_;
    }

    T * get()
    {
      return ptr_;
    }

    void clear()
    {
      if (ptr_) {
        std::free(ptr_);
        ptr_ = nullptr;
      }
    }

    void swap(MallocArray & oth)
    {
      std::swap(ptr_,oth.ptr_);
    }

    void resize(std::size_t n)
    {
      ptr_ = static_cast<T*>(crealloc(ptr_, sizeof(T)*n));
    }

    T operator[](std::size_t i) const
    {
      return ptr_[i];
    }

    T & operator[](std::size_t i)
    {
      return ptr_[i];
    }

  private:
    T * ptr_;
  };


  class BitArrayView {
  public:
    explicit BitArrayView(UInt * data) : data_(data) {}
    BitArrayView(BitArrayView & oth) : data_(oth.data_) {}

    void set(unsigned i)
    {
      data_[(i >> 6)] = data_[(i >> 6)] | (UInt(1) << (i & BitMask<6>::val));
    }

    void unset(unsigned i)
    {
      data_[(i >> 6)] = data_[(i >> 6)] & ~(UInt(1) << (i & BitMask<6>::val));
    }

    bool get(unsigned i) const
    {
      return data_[(i >> 6)] & (UInt(1) << (i & BitMask<6>::val));
    }

  private:
    friend class ConstBitArrayView;

  private:
    UInt * data_;
  };

  class ConstBitArrayView {
  public:
    explicit ConstBitArrayView(const UInt * data) : data_(data) {}
    ConstBitArrayView(ConstBitArrayView & oth) : data_(oth.data_) {}
    ConstBitArrayView(BitArrayView & oth) : data_(oth.data_) {}

    bool get(unsigned i) const
    {
      return data_[(i >> 6)] & (UInt(1) << (i & BitMask<6>::val));
    }

    const UInt * data() const
    {
      return data_;
    }

  private:
    const UInt * data_;
  };

  inline unsigned bit_array_u64size(unsigned size)
  {
    return (size >> 6) + !!(size & BitMask<6>::val);
  }

  inline unsigned bit_array_nonzeroes(ConstBitArrayView v, unsigned size)
  {
    unsigned ret = 0;
    for (unsigned i=0; i<size; ++i)
      if (v.get(i))
        ++ret;
    return ret;
  }

} // namespace fflow


#endif // FFLOW_COMMON_HH

