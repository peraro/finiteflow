#ifndef FFLOW_GLOBAL_HH
#define FFLOW_GLOBAL_HH

#include <new>

namespace fflow {

// Defines a variable whose destructor is never called.  The type T
// must be default-constructible.  It is used for global variables, as
// calling their destructors can cause issues, especially when
// unloading dynamically loaded libraries.
template<typename T>
class Global {

public:

  Global()
  {
    new(&bytes_) T();
  }

  Global(const Global & oth) = delete;
  Global & operator=(const Global & oth) = delete;

  const T & operator*() const
  {
    return *(reinterpret_cast<const T *>(bytes_));
  }

  T & operator*()
  {
    return *((reinterpret_cast<T *>(bytes_)));
  }

private:
  alignas(T) char bytes_[sizeof(T)];
};

}

#endif  // FFLOW_GLOBAL_HH
