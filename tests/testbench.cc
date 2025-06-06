#include <chrono>
#include <vector>
#include <random>
#include <fflow/primes.hh>
#include <fflow/common.hh>
#include <fflow/matrix.hh>
using namespace fflow;

int main(int argc, char** argv)
{
  std::cout << sizeof(LSVar) << std::endl;

  const unsigned w_size = 10000000;
  unsigned t_size = 0;

  if (argc > 3 || argc < 2) {
    std::cout << format("Usage: {} t_size [--shoup]", argv[0]) << std::endl;
    return -1;
  }
  t_size = std::atoi(argv[1]);

  bool shoup = false;
  if (argc == 3 && std::string(argv[2]) == std::string("--shoup"))
    shoup = true;

  std::vector<UInt> w(w_size);
  std::vector<UInt> t(t_size);

  Mod mod;
  {
    std::mt19937_64 gen;
    std::uniform_int_distribution<UInt> primei(0, 200);
    mod = Mod(BIG_UINT_PRIMES[primei(gen)]);
    std::uniform_int_distribution<UInt> dis(0, mod.n()-1);
    for (unsigned j=0; j<w_size; ++j)
      w[j] = dis(gen);
    for (unsigned j=0; j<t_size; ++j)
      t[j] = dis(gen);
  }

  if (shoup) {

    UInt res = 0;
    auto start = std::chrono::system_clock::now();
    for (unsigned j=0; j<w_size; ++j) {
      UInt wp = precomp_mul_shoup(w[j], mod);
      for (unsigned k=0; k<t_size; ++k)
        res = add_mod(res, mul_mod_shoup(t[k], w[j], wp, mod), mod);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    dbgprint(format("Shoup eval: {} sec. (res={})",
                    elapsed_seconds.count(), res));

  } else {

    UInt res = 0;
    auto start = std::chrono::system_clock::now();
    for (unsigned j=0; j<w_size; ++j) {
      for (unsigned k=0; k<t_size; ++k)
        res = apbc_mod(res, t[k], w[j], mod);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    dbgprint(format("Normal eval: {} sec. (res={})",
                    elapsed_seconds.count(), res));

  }

  return 0;
}
