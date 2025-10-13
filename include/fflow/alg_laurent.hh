#ifndef FFLOW_ALG_LAURENT_HH
#define FFLOW_ALG_LAURENT_HH

#include <fflow/subgraph_reconstruct.hh>

namespace fflow {

  typedef BaseSubgraphUniRecData LaurentExpansionData;

  class LaurentExpansion : public BaseSubgraphUniRec {
  public:
    LaurentExpansion()
    {
      laurent_expand_ = true;
    }
  };

} // namespace fflow

#endif // FFLOW_ALG_LAURENT_HH
