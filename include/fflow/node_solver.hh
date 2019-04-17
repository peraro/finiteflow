#ifndef FFLOW_ALG_NODE_SOLVER_HH
#define FFLOW_ALG_NODE_SOLVER_HH

#include<fflow/alg_linear_solver.hh>

namespace fflow {

  class NodeDenseSolver;
  class NodeSparseSolver;


  struct NodeDenseSolverData : public DenseLinearSolverData {};

  class NodeDenseSolver : public DenseLinearSolver {
  public:

    void init_node_solver(std::size_t neqs, std::size_t nvars,
                          unsigned * needed_vars, unsigned needed_size,
                          NodeDenseSolverData & data)
    {
      nparsin.resize(1);
      nparsin[0] = neqs*(nvars+1);
      init(neqs, nvars, needed_vars, needed_size, data);
    }

    virtual Ret fill_matrix(Context * ctxt,
                            std::size_t n_rows, const std::size_t rows[],
                            AlgInput xi[], Mod mod,
                            AlgorithmData * data,
                            MatrixView & m) const override;

    virtual AlgorithmData::Ptr
    clone_data(const AlgorithmData * data) const override;
  };

} // namespace fflow


#endif // FFLOW_ALG_NODE_SOLVER_HH
