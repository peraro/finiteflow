#include <cstdio>
#include <mutex>
#include <fflow/debug.hh>

namespace {

  fflow::DBGPrint * dbgprint_ = &fflow::stdout_dbgprint;
  fflow::DBGPrint * logerr_ = &fflow::stderr_dbgprint;
  std::mutex logerr_mutex_;

} // namespace


namespace fflow {

  StdOutDBGPrint stdout_dbgprint;
  StdErrDBGPrint stderr_dbgprint;
  NullDBGPrint null_dbgprint;

  void StdOutDBGPrint::print(const std::string & msg)
  {
    std::puts(msg.c_str());
  }

  void StdErrDBGPrint::print(const std::string & msg)
  {
    std::fprintf(stderr, "FF_LOG: %s\n", msg.c_str());
  }

  void NullDBGPrint::print(const std::string & msg)
  {
    (void)(msg);
  }

  DBGPrintSet::DBGPrintSet(DBGPrint & p) : previous(dbgprint_)
  {
    dbgprint_ = &p;
  }

  DBGPrintSet::~DBGPrintSet()
  {
    dbgprint_ = previous;
  }

  void dbgprint(const std::string & msg)
  {
    dbgprint_->print(msg);
  }

  LogErrorSet::LogErrorSet(DBGPrint & p) : previous(logerr_)
  {
    logerr_ = &p;
  }

  LogErrorSet::~LogErrorSet()
  {
    logerr_ = previous;
  }

  void logerr(const std::string & msg)
  {
    std::lock_guard<std::mutex> lock(logerr_mutex_);
    logerr_->print(msg);
  }

} // namespace fflow
