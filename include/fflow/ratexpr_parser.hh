#include <vector>
#include <fflow/alg_functions.hh>

namespace fflow {

  typedef const char * ExprCStr;

  Ret parse_ratexpr_list(unsigned npars,
                         const std::string & vars_prefix,
                         const std::string vars[],
                         const ExprCStr inputs[], const unsigned len[],
                         unsigned nfunctions,
                         std::vector<std::vector<Instruction>> & bytecode,
                         std::vector<MPRational> & bignums,
                         std::vector<AnalyticExpression::VarPow> & varpows);

}
