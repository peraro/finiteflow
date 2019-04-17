#ifndef FFLOW_ALG_MP_RECONSTRUCTION_HH
#define FFLOW_ALG_MP_RECONSTRUCTION_HH

#include <fflow/alg_reconstruction.hh>
#include <fflow/mp_multivariate_reconstruction.hh>

namespace fflow {


  struct LearningOptions {
    unsigned n_singular = 0;
    unsigned prime_no = 0;
  };


  struct ReconstructionOptions {
    unsigned n_checks = RatFunReconstruction::DEFAULT_N_CHECKS;
    unsigned n_uchecks = RatFunReconstruction::DEFAULT_N_UCHECKS;
    unsigned n_singular = 0;
    unsigned start_mod = 0;
    unsigned max_primes = 1;
    unsigned max_deg = RatFunReconstruction::DEFAULT_MAX_DEG;
    unsigned dbginfo = 0;
  };


  inline void set_rec_options(const ReconstructionOptions & opt,
                              const UInt shift[],
                              RatFunReconstruction & rec)
  {
    rec.setShift(shift);
    rec.n_checks = opt.n_checks;
    rec.n_uchecks = opt.n_uchecks;
    rec.n_singular = opt.n_singular;
    rec.setMaxDegree(opt.max_deg);
  }


  inline void set_mp_rec_options(const ReconstructionOptions & opt,
                                 const UInt shift[],
                                 MPRatFunReconstruction & rec)
  {
    rec.setShift(shift);
    rec.n_checks = opt.n_checks;
    rec.n_uchecks = opt.n_uchecks;
    rec.n_singular = opt.n_singular;
    rec.start_mod = opt.start_mod;
    rec.max_primes = opt.max_primes;
    rec.setMaxDegree(opt.max_deg);
  }


  Ret algorithm_get_degrees(Algorithm & alg,
                            AlgorithmData * data,
                            Context * ctxt,
                            const UInt shift[],
                            const ReconstructionOptions & opt,
                            unsigned numdeg[], unsigned dendeg[]);

  Ret algorithm_get_var_degrees(Algorithm & alg,
                                AlgorithmData * data,
                                Context * ctxt,
                                unsigned var,
                                const UInt shift[],
                                const ReconstructionOptions & opt,
                                RatFunVarDegrees degs[]);

  Ret algorithm_dump_degree_info(const char * filename,
                                 unsigned nparsin,
                                 unsigned nparsout,
                                 const unsigned numdeg[],
                                 const unsigned dendeg[],
                                 const RatFunVarDegrees degs[]);

  Ret algorithm_load_degree_info(const char * filename,
                                 unsigned nparsin,
                                 unsigned nparsout,
                                 unsigned numdeg[],
                                 unsigned dendeg[],
                                 RatFunVarDegrees degs[]);

  Ret algorithm_npars_from_degree_info(const char * filename,
                                       unsigned & npars_in,
                                       unsigned & npars_out);

  Ret algorithm_load_degree_info_for_var(const char * filename,
                                         unsigned nparsin,
                                         unsigned nparsout,
                                         unsigned var,
                                         unsigned & numdeg,
                                         unsigned & dendeg,
                                         RatFunVarDegrees & degs);

  void algorithm_generate_sample_points(GenerateSamplePoints & samples,
                                        unsigned nparsin, unsigned nparsout,
                                        const UInt shift[],
                                        const ReconstructionOptions & opt,
                                        const RatFunVarDegrees degs[],
                                        const unsigned numdegs[],
                                        const unsigned dendegs[]);

  void algorithm_verify_sample_points(std::unique_ptr<UInt[]> * samples,
                                      unsigned nsamples,
                                      unsigned nparsin, unsigned nparsout,
                                      const std::vector<unsigned> & needed,
                                      const UInt shift[],
                                      const ReconstructionOptions & opt,
                                      const RatFunVarDegrees degs[],
                                      const unsigned numdegs[],
                                      const unsigned dendegs[]);

  Ret algorithm_reconstruct(const UIntCache & cache,
                            unsigned nparsin, unsigned,
                            unsigned idx,
                            const UInt shift[],
                            const ReconstructionOptions & opt,
                            const RatFunVarDegrees degs[],
                            const unsigned numdeg[], const unsigned dendeg[],
                            MPReconstructedRatFun & recf);

  Ret algorithm_sparse_reconstruct(const UIntCache & cache,
                                   unsigned nparsin, unsigned,
                                   unsigned idx,
                                   const UInt shift[],
                                   const ReconstructionOptions & opt,
                                   const RatFunVarDegrees degs[],
                                   const unsigned numdeg[],
                                   const unsigned dendeg[],
                                   MPReconstructedRatFun & recf);

  Ret algorithm_reconstruct_univariate(const Algorithm & alg,
                                       AlgorithmData * data,
                                       Context * ctxt,
                                       const UInt shift[],
                                       const ReconstructionOptions & opt,
                                       MPReconstructedRatFun res[]);

  Ret algorithm_reconstruct_numeric(const Algorithm & alg,
                                    AlgorithmData * data,
                                    Context * ctxt,
                                    const ReconstructionOptions & opt,
                                    MPRational res[]);

  Ret algorithm_reconstruct_from_evals(const std::string ev_files[],
                                       unsigned nfiles,
                                       const char * degs_file,
                                       unsigned nparsin, unsigned nparsout,
                                       unsigned idx,
                                       const UInt shift[],
                                       const ReconstructionOptions & opt,
                                       MPReconstructedRatFun & recf);

} // namespace fflow

#endif // FFLOW_ALG_MP_RECONSTRUCTION_HH
