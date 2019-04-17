#ifndef FFLOW_NO_DBG

#include <iostream>
#include <fflow/debug.hh>

namespace {

  fflow::StdOutDBGPrint stdout_dbgprint_;

  fflow::DBGPrint * dbgprint_ = &stdout_dbgprint_;

} // namespace


namespace fflow {

  void StdOutDBGPrint::print(const std::string & msg)
  {
    std::cout << msg << std::endl;
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

} // namespace fflow

#endif // !FFLOW_NO_DBG
