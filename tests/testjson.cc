#include <iostream>
#include <chrono>
#include <fflow/graph.hh>
#include <fflow/analytic_solver.hh>
#include <fflow/alg_functions.hh>
#include <fflow/alg_lists.hh>
#include <fflow/json.hh>
using namespace fflow;

class AbsoluteTiming {
public:
  AbsoluteTiming(const std::string & str) : start()
  {
    std::cout << str << std::endl;
    start = std::chrono::system_clock::now();
  }

  void end()
  {
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "- completed in " << elapsed_seconds.count()
              << " seconds" << std::endl;
  }

  ~AbsoluteTiming()
  {
    end();
  }

private:
  std::chrono::time_point<std::chrono::system_clock> start;
};

int main(int argc, char** argv)
{
  if (argc != 3) {
    std::cout << format("Usage: {} input output", argv[0]) << std::endl;
    return -1;
  }

  Ret ret = SUCCESS;

  Session session;

  unsigned graphid = session.new_graph();
  Graph * graph = session.graph(graphid);
  unsigned needed_workspace = 0;

  // derivatives
  std::unique_ptr<AnalyticFunction> alg_func(new AnalyticFunction());
  std::unique_ptr<AnalyticFunctionData> data_func(new AnalyticFunctionData());
  {
    AbsoluteTiming time("Parse JSON derivdata");
    ret = json_sparse_ratfun(argv[1], *alg_func, *data_func,
                             needed_workspace);
  }
  if (ret != SUCCESS) {
    std::cout << format("Alg. func. parsing failed with status {}", ret)
              << std::endl;
    return -1;
  }
  if (alg_func->nparsout == 0) {
    std::cout << format("No function specified") << std::endl;
    return 0;
  }

  unsigned input = graph->set_input_vars(alg_func->nparsin[0]);

  session.main_context()->ww.ensure_size(needed_workspace);
  unsigned functions = graph->new_node(std::move(alg_func),
                                       std::move(data_func), &input);

  session.set_output_node(graphid, functions);

  std::vector<MPReconstructedRatFun> res(graph->nparsout);
  ReconstructionOptions opt;
  opt.max_primes = 5;
  unsigned nthreads = 0; // Use 0 to select autamatically
  {
    AbsoluteTiming time("Full reconstruction");
    ret = session.full_reconstruction(graphid, res.data(),
                                      nthreads,
                                      opt);
  }
  if (ret != SUCCESS) {
    std::cout << format("Reconstruction failed with status {}", ret)
              << std::endl;
    return -1;
  }

  ret =json_write_ratfun(argv[2], res.data(), res.size(), graph->nparsin[0]);
  if (ret != SUCCESS) {
    std::cout << format("Couldn't write output to file {}", argv[2])
              << std::endl;
    return -1;
  }

  return 0;
}
