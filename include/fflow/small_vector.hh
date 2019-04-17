#ifndef FFLOW_SMALL_VECTOR_HH
#define FFLOW_SMALL_VECTOR_HH

#include <type_traits>
#include <algorithm>

namespace fflow {

  template <typename T, std::size_t N>
  class SmallVector {
  public:

    static_assert(std::is_pod<T>::value,
                  "fflow::SmallVector can only be used with POD types");

    SmallVector() : data_(), size_(0), cap_(N)
    {
      data_ = buff_;
    }

    explicit SmallVector(std::size_t n) : data_(), size_(0), cap_(N)
    {
      data_ = buff_;
      resize(n);
    }

    SmallVector(const SmallVector & oth) : data_(), size_(0), cap_(N)
    {
      data_ = buff_;
      *this = oth;
    }

    ~SmallVector()
    {
      if (data_ != buff_)
        delete [] data_;
    }

    SmallVector & operator=(const SmallVector & oth)
    {
      resize(oth.size());
      std::copy(oth.data_, oth.data_+size_, data_);
      return *this;
    }

    std::size_t size() const
    {
      return size_;
    }

    T & operator[] (std::size_t i)
    {
      return data_[i];
    }

    const T & operator[] (std::size_t i) const
    {
      return data_[i];
    }

    void reserve(std::size_t n);

    void resize(std::size_t n)
    {
      if (n > cap_)
        reserve(n);
      size_ = n;
    }

    void push_back(const T & el)
    {
      resize(size_+1);
      data_[size_-1] = el;
    }

    void clear()
    {
      size_ = 0;
    }

    T * begin()
    {
      return data_;
    }

    T * end()
    {
      return data_ + size_;
    }

    const T * begin() const
    {
      return data_;
    }

    const T * end() const
    {
      return data_ + size_;
    }

    T * data()
    {
      return data_;
    }

    const T * data() const
    {
      return data_;
    }

  private:
    T * data_;
    unsigned size_, cap_;
    T buff_[N];
  };


  template <typename T, std::size_t N>
  void SmallVector<T,N>::reserve(std::size_t n)
  {
    if (n <= cap_)
      return;

    unsigned new_cap = std::max(cap_*2, unsigned(n));
    T * new_data = new T[new_cap];
    std::copy(data_, data_+size_, new_data);

    if (data_ != buff_)
      delete [] data_;

    data_ = new_data;
    cap_ = new_cap;
  }


} // namespace fflow

#endif // FFLOW_SMALL_VECTOR_HH
