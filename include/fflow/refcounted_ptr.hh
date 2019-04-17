#ifndef FFLOW_REFCOUNTED_PTR_HH
#define FFLOW_REFCOUNTED_PTR_HH

#include <cstddef>
#include <new>

namespace fflow {

  // T must be a type with inc_ref_count() and dec_ref_count()
  // methods, which must update and return the ref. count after
  // increasing/decreasing it
  template <typename T>
  class RefCountedPtr {
  public:

    RefCountedPtr() : ptr_(nullptr) {}

    explicit RefCountedPtr(std::nullptr_t) : ptr_(nullptr) {}

    explicit RefCountedPtr(T * ptr) : ptr_(ptr)
    {
      if (ptr)
        ptr->inc_ref_count();
    }

    RefCountedPtr(RefCountedPtr & oth) : ptr_(oth.ptr_)
    {
      if (ptr_)
        ptr_->inc_ref_count();
    }

    RefCountedPtr(RefCountedPtr && oth) : ptr_(oth.ptr_)
    {
      oth.ptr_ = nullptr;
    }

    ~RefCountedPtr()
    {
      remove_ref_();
    }

    void swap(RefCountedPtr & oth)
    {
      std::swap(ptr_, oth.ptr_);
    }

    void reset()
    {
      remove_ref_();
      ptr_ = nullptr;
    }

    void reset(T * ptr)
    {
      RefCountedPtr newptr(ptr);
      newptr.swap(*this);
    }

    RefCountedPtr & operator=(RefCountedPtr & oth)
    {
      RefCountedPtr newptr(oth);
      newptr.swap(*this);
      return *this;
    }

    RefCountedPtr & operator=(RefCountedPtr && oth)
    {
      RefCountedPtr newptr(std::move(oth));
      newptr.swap(*this);
      return *this;
    }

    T * get()
    {
      return ptr_;
    }

    const T * get() const
    {
      return ptr_;
    }

    T & operator* ()
    {
      return *ptr_;
    }

    const T & operator* () const
    {
      return *ptr_;
    }

    T * operator-> ()
    {
      return ptr_;
    }

    const T * operator-> () const
    {
      return ptr_;
    }

  private:

    void remove_ref_()
    {
      if (ptr_ && !ptr_->dec_ref_count())
        delete ptr_;
    }

  private:
    T * ptr_;
  };

#define FFLOW_IMPLEM_REF_COUNT(refcount_)       \
  unsigned inc_ref_count()                      \
  {                                             \
    return ++refcount_;                         \
  }                                             \
  unsigned dec_ref_count()                      \
  {                                             \
    return --refcount_;                         \
  }                                             \
  unsigned ref_count() const                    \
  {                                             \
    return refcount_;                           \
  }

} // namespace fflow


#endif // FFLOW_REFCOUNTED_PTR_HH
