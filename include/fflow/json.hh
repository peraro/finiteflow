#ifndef FFLOW_JSON_HH
#define FFLOW_JSON_HH

#include <string>
#include <fflow/common.hh>

namespace fflow {

  class AnalyticSparseSolver;
  struct AnalyticSparseSolverData;
  class AnalyticFunction;
  struct AnalyticFunctionData;
  class MPReconstructedRatFun;

  std::string read_file(const char * filename);

  Ret json_sparse_system(const char * json_info_file,
                         AnalyticSparseSolver & sys,
                         AnalyticSparseSolverData & data,
                         unsigned & required_workspace);

  Ret json_sparse_system_with_info(const char * json_info_file,
                                   const char * binary_info_file,
                                   AnalyticSparseSolver & sys,
                                   AnalyticSparseSolverData & data,
                                   unsigned & required_workspace);

  Ret json_sparse_ratfun(const char * file,
                         AnalyticFunction & fun,
                         AnalyticFunctionData & data,
                         unsigned & required_workspace);

  Ret json_integer_list(const char * filename,
                        std::vector<unsigned> & list);

  Ret json_write_ratfun(const char * file,
                        const MPReconstructedRatFun * fun,
                        unsigned len,
                        unsigned nvars);

  Ret json_write_integer_list(const char * filename,
                              const unsigned * list,
                              unsigned len);

} // namespace fflow

#endif // FFLOW_JSON_HH
