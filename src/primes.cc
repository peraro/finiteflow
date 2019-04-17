#include <fflow/primes.hh>
#include <cmath>

namespace fflow {


  // Simple implementation of the segmented sieve of Eratosthenes,
  // based on the one of Kim Walisch which is in the public domain and
  // can be obtained at http://primesieve.org/segmented_sieve.html.

  void PrimeVector::segmented_sieve_(UInt limit, std::size_t segment_size)
  {
    UInt sqrt = static_cast<UInt>(std::sqrt(limit));
    UInt s = 2;
    UInt n = 3;

    // vector used for sieving
    std::vector<char> sieve(segment_size);

    // generate small primes <= sqrt
    std::vector<char> is_prime(sqrt + 1, 1);
    for (UInt i = 2; i * i <= sqrt; i++)
      if (is_prime[i])
        for (UInt j = i * i; j <= sqrt; j += i)
          is_prime[j] = 0;

    std::vector<UInt> primes;
    std::vector<UInt> next;

    if (2u <= limit)
      primes_.push_back(2);

    for (UInt low = 0; low <= limit; low += segment_size) {

        std::fill(sieve.begin(), sieve.end(), 1);

        // current segment = interval [low, high]
        UInt high = std::min(low + segment_size - 1, limit);

        // store small primes needed to cross off multiples
        for (; s * s <= high; s++) {
            if (is_prime[s]) {
                primes.push_back(s);
                next.push_back(s * s - low);
            }
        }
        // sieve the current segment
        for (std::size_t i = 1; i < primes.size(); i++) {
            UInt j = next[i];
            for (std::size_t k = primes[i] * 2; j < segment_size; j += k)
              sieve[j] = 0;
            next[i] = j - segment_size;
        }

        for (; n <= high; n += 2)
          if (sieve[n - low]) // n is a prime
            primes_.push_back(n);
      }

  }


} // namespace fflow
