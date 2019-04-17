#ifndef FFLOW_RATFUN_PARSER_HH
#define FFLOW_RATFUN_PARSER_HH

#include <string>
#include <fflow/common.hh>

namespace fflow {

  class Monomial;
  class MPRational;

  Ret parse_rat_fun(const std::string vars[], unsigned nvars,
                    const char * start, const char * end,
                    std::vector<Monomial> & num,
                    std::vector<MPRational> & cnum,
                    std::vector<Monomial> & den,
                    std::vector<MPRational> & cden);

  Ret parse_rat_fun_mod(const std::string vars[], unsigned nvars,
                        const char * start, const char * end,
                        UInt mod,
                        std::vector<Monomial> & num,
                        std::vector<Monomial> & den);

} // namespace fflow

#endif // FFLOW_RATFUN_PARSER_HH
