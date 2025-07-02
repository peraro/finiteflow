#ifndef FFLOW_PRIMES_HH
#define FFLOW_PRIMES_HH

#include <vector>
#include <cmath>
#include <fflow/common.hh>

namespace fflow {


  // The 2048 largest primes which fit in a 63-bit integer.  They can
  // be generated with the following mathematica command
  //
  //    Select[2^63 - Range[89121], PrimeQ]
  //
  // The Mod argument in the routines of these library is assumed to
  // be initialized with one of these integers.
  const unsigned BIG_UINT_PRIMES_SIZE = 2048;
  const unsigned DEFAULT_MAX_REC_PRIMES = BIG_UINT_PRIMES_SIZE - 1;
  extern const UInt BIG_UINT_PRIMES[BIG_UINT_PRIMES_SIZE];
  //const UInt LARGEST_PRIME = BIG_UINT_PRIMES[0];
  //const UInt SMALLEST_PRIME = BIG_UINT_PRIMES[BIG_UINT_PRIMES_SIZE-1];


  inline UInt prime_no(unsigned i)
  {
    return BIG_UINT_PRIMES[i % BIG_UINT_PRIMES_SIZE];
  }


  // A vector of prime numbers computed using the segmented sieve
  // method
  class PrimeVector {
  public:

    static const std::size_t DEFAULT_SEGMENT_SIZE = 32768;

    PrimeVector() : primes_() {}

    PrimeVector(UInt limit, std::size_t segmentsize = DEFAULT_SEGMENT_SIZE)
      : primes_()
    {
      segmented_sieve_(limit, segmentsize);
    }

    void generate_primes(UInt limit,
                         std::size_t segmentsize = DEFAULT_SEGMENT_SIZE)
    {
      clear();
      segmented_sieve_(limit, segmentsize);
    }

    void clear()
    {
      primes_.clear();
    }

    std::size_t size()
    {
      return primes_.size();
    }

    UInt operator[] (std::size_t i) const
    {
      return primes_[i];
    }

    UInt * data()
    {
      return primes_.data();
    }

    const UInt * data() const
    {
      return primes_.data();
    }

    operator bool()
    {
      return !primes_.empty();
    }

  private:
    void segmented_sieve_(UInt limit,
                          std::size_t segmentsize = DEFAULT_SEGMENT_SIZE);

  private:
    std::vector<UInt> primes_;
  };


  // Utility function.  Satisfies prime_pi(invprime_bound(n)) > n,
  // where prime_pi is the prime counting function, by some
  // significant but not too large margin.  Useful if you need n
  // primes, in which case you use safely use invprime_bound(n) as
  // upper limit in the generator.
  inline UInt invprime_bound(UInt n_primes)
  {
    return UInt(1.25 * n_primes * std::log(n_primes) + 20);
  }



} // namespace fflow


#endif // FFLOW_PRIMES_HH
