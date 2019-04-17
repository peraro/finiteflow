#ifndef FFLOW_SHARED_ARRAY_HH
#define FFLOW_SHARED_ARRAY_HH

#include <memory>

namespace fflow {

  template <typename T>
  class SharedArray {
  public:
    SharedArray() : ptr_(nullptr) {}
    explicit SharedArray(std::nullptr_t) : ptr_(nullptr) {}
    explicit SharedArray(T * ptr) : ptr_(ptr, std::default_delete<T[]>()) {}

    SharedArray(const SharedArray & oth) : ptr_(oth.ptr_) {}
    SharedArray(SharedArray && oth) : ptr_(std::move(oth.ptr_)) {}

    SharedArray & operator=(const SharedArray & oth)
    {
      ptr_ = oth.ptr_;
      return *this;
    }

    SharedArray & operator=(SharedArray && oth)
    {
      ptr_ = std::move(oth.ptr_);
      return *this;
    }

    void reset()
    {
      ptr_.reset();
    }

    void reset(T * ptr)
    {
      ptr_.reset(ptr, std::default_delete<T[]>());
    }

    void swap(SharedArray & oth)
    {
      ptr_.swap(oth);
    }

    const T * get() const
    {
      return ptr_.get();
    }

    T * get()
    {
      return ptr_.get();
    }

    const T & operator[] (std::size_t i) const
    {
      return ptr_.get()[i];
    }

    T & operator[] (std::size_t i)
    {
      return ptr_.get()[i];
    }

  private:
    std::shared_ptr<T> ptr_;
  };

} // namespace fflow


#endif // FFLOW_SHARED_ARRAY_HH
