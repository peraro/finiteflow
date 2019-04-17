#ifndef FFLOW_DEBUG_HH
#define FFLOW_DEBUG_HH

#include <string>

namespace fflow {

  // inherit from this in order to create a custom class for printing
  // debug info
  class DBGPrint {
  public:
    virtual void print(const std::string & msg) = 0;
  };

  // default DBGPrint just prints to stdout
  class StdOutDBGPrint : public DBGPrint {
  public:
    StdOutDBGPrint() {}
    virtual void print(const std::string & msg);
  };


  // Create an instance of this, passing a DBGPrint object, in order
  // to set it for the current scope.  Call the static method
  // DBGPrintSet::global in order to set it on the global scope.
  class DBGPrintSet {
  public:

    explicit DBGPrintSet(DBGPrint & p);

    ~DBGPrintSet();

    static void global(DBGPrint & p);

  private:
    DBGPrint * previous;
  };


  // Calls DBGPrint::print on the currently active DBGPrint object
  void dbgprint(const std::string & msg);

} // namespace fflow

#endif // FFLOW_DEBUG_HH
