#ifndef FFLOW_MP_TYPES_HH
#define FFLOW_MP_TYPES_HH

#include <cstring>
#include <cstddef>
#include <gmp.h>
#include <fflow/common.hh>

namespace fflow {

  class MPInt {
  public:

    MPInt()
    {
      mpz_init(z_);
    }

    MPInt(const MPInt & oth)
    {
      mpz_init_set(z_, oth.z_);
    }

    MPInt(MPInt && oth)
    {
      mpz_init(z_);
      swap(oth);
    }

    explicit MPInt(UInt n)
    {
      mpz_init_set_ui(z_, n);
    }

    explicit MPInt(const char * str, unsigned base = 10)
    {
      mpz_init_set_str(z_, str, base);
    }

    ~MPInt()
    {
      mpz_clear(z_);
    }

    MPInt & operator= (const MPInt & oth)
    {
      mpz_set(z_, oth.z_);
      return *this;
    }

    MPInt & operator= (MPInt && oth)
    {
      swap(oth);
      return *this;
    }

    MPInt & operator= (UInt n)
    {
      mpz_set_ui(z_, n);
      return *this;
    }

    MPInt & operator= (const char * str)
    {
      set(str);
      return *this;
    }

    int set(const char * str, unsigned base = 10)
    {
      return mpz_set_str(z_, str, base);
    }

    Int to_int() const
    {
      return mpz_get_si(z_);
    }

    UInt to_uint() const
    {
      return mpz_get_ui(z_);
    }

    mpz_t & get()
    {
      return z_;
    }

    const mpz_t & get() const
    {
      return z_;
    }

    int cmp(const MPInt & oth) const
    {
      return mpz_cmp(z_, oth.z_);
    }

    int cmp(UInt oth) const
    {
      return mpz_cmp_ui(z_, oth);
    }

    int sign() const
    {
      return mpz_sgn(z_);
    }

    void swap(MPInt & oth)
    {
      mpz_swap(z_, oth.z_);
    }

    template <typename WriteStream>
    void print(WriteStream & os, unsigned base = 10) const;

  private:
    mpz_t z_;
  };


  class MPRational {
  public:

    MPRational()
    {
      mpq_init(q_);
    }

    MPRational(const MPRational & oth)
    {
      mpq_init(q_);
      mpq_set(q_, oth.q_);
    }

    MPRational(MPRational && oth)
    {
      mpq_init(q_);
      swap(oth);
    }

    MPRational(const MPInt & z)
    {
      mpq_init(q_);
      mpq_set_z(q_, z.get());
    }

    explicit MPRational(Int num, Int den=1)
    {
      mpq_init(q_);
      mpq_set_si(q_, num, den);
      canonicalize();
    }

    explicit MPRational(const Rational & q)
    {
      mpq_init(q_);
      mpq_set_si(q_, q.num, q.den);
      canonicalize();
    }

    explicit MPRational(const char * str, unsigned base = 10)
    {
      mpq_init(q_);
      set(str, base);
    }

    MPRational & operator= (const MPRational & oth)
    {
      mpq_set(q_, oth.q_);
      return *this;
    }

    MPRational & operator= (MPRational && oth)
    {
      swap(oth);
      return *this;
    }

    MPRational & operator= (const MPInt & z)
    {
      mpq_set_z(q_, z.get());
      return *this;
    }

    MPRational & operator =(const Rational & q)
    {
      mpq_set_si(q_, q.num, q.den);
      canonicalize();
      return *this;
    }

    MPRational & operator =(Int z)
    {
      mpq_set_si(q_, z, 1);
      canonicalize();
      return *this;
    }

    MPRational & operator= (const char * str)
    {
      set(str);
      return *this;
    }

    ~MPRational()
    {
      mpq_clear(q_);
    }

    int set(const char * str, unsigned base = 10)
    {
      int ret = mpq_set_str(q_, str, base);
      canonicalize();
      return ret;
    }

    mpq_t & get()
    {
      return q_;
    }

    const mpq_t & get() const
    {
      return q_;
    }

    void canonicalize()
    {
      mpq_canonicalize(q_);
    }

    void swap(MPRational & oth)
    {
      mpq_swap(q_, oth.q_);
    }

    int cmp(const MPRational & oth) const
    {
      return mpq_cmp(q_, oth.q_);
    }

#if 0
    int cmp(const MPInt & oth) const
    {
      return mpq_cmp_z(q_, oth.get());
    }
#endif

    int cmp(Int num, Int den=1) const
    {
      return mpq_cmp_si(q_, num, den);
    }

    int cmp(const Rational & q) const
    {
      Int num = q.num, den = q.den;
      return mpq_cmp_si(q_, num, den);
    }

    int sign() const
    {
      return mpq_sgn(q_);
    }

    void to_int(Int & num, Int & den) const
    {
      num = mpz_get_si(mpq_numref(q_));
      den = mpz_get_si(mpq_denref(q_));
    }

    bool operator== (const MPRational & oth) const
    {
      return mpq_equal(q_, oth.q_);
    }

    template <typename WriteStream>
    void print(WriteStream & os, unsigned base = 10) const;

  private:
    mpq_t q_;
  };


  template <typename WriteStream>
  inline void MPInt::print(WriteStream & os, unsigned base) const
  {
    void (*freefunc) (void *, size_t);
    mp_get_memory_functions (0, 0, &freefunc);

    char * str = mpz_get_str(0, base, z_);

    os << str;

    (*freefunc)(str, std::strlen(str)+1);
  }

  template <typename WriteStream>
  inline void MPRational::print(WriteStream & os, unsigned base) const
  {
    void (*freefunc) (void *, size_t);
    mp_get_memory_functions (0, 0, &freefunc);

    char * str = mpq_get_str(0, base, q_);

    os << str;

    (*freefunc)(str, std::strlen(str)+1);
  }


  inline std::ostream & operator<< (std::ostream & os, const MPInt & z)
  {
    z.print(os);
    return os;
  }

  inline std::ostream & operator<< (std::ostream & os, const MPRational & q)
  {
    q.print(os);
    return os;
  }

} // namespace fflow

#endif // FFLOW_MP_TYPES_HH
