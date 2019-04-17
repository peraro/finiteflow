#include <vector>
#include <fflow/common.hh>

namespace fflow {

  class VarValues {
  public:
    enum {END=-1};

    VarValues(std::size_t n_vars, unsigned rank)
      : v_(n_vars), r_(rank), sum_(0) {}

    unsigned operator[] (std::size_t i) const
    {
      return v_[i];
    }

#if 0
    unsigned & operator[] (std::size_t i)
    {
      return v_[i];
    }
#endif

    void reset()
    {
      for (auto & i : v_)
        i = 0;
      sum_ = 0;
    }

    void reset(std::size_t n_vars, unsigned rank)
    {
      if (n_vars != v_.size())
        v_.resize(n_vars);
      r_ = rank;
      reset();
    }

    unsigned tdeg() const
    {
      return r_+1-sum_;
    }

    // Returns the variable i which got increased, or END
    int next();

    std::size_t nvars()
    {
      return v_.size();
    }

  private:
    std::vector<unsigned> v_;
    unsigned r_, sum_;
  };

  inline
  int VarValues::next()
  {
    unsigned i = nvars()-1;

    while (i>0) {

      if (sum_ < r_) {
        ++v_[i];
        ++sum_;
        return i;
      } else {
        sum_ -= v_[i];
        v_[i]=0;
      }

      --i;
    }
    return END;
  }


  class RatFunVarDegrees {
  public:

    RatFunVarDegrees() : x_(nullptr) {}

    RatFunVarDegrees(RatFunVarDegrees && oth) : x_(std::move(oth.x_)) {}

    void resize(std::size_t nvars)
    {
      x_.reset(new std::size_t[4*nvars]());
    }

    void copy(const RatFunVarDegrees & oth, unsigned nvars)
    {
      if (oth.x_.get() != nullptr) {
        resize(nvars);
        std::copy(oth.x_.get(), oth.x_.get() + 4*nvars, x_.get());
      }
    }

    operator bool() const
    {
      return x_.get();
    }

    RatFunVarDegrees & operator=(RatFunVarDegrees && oth)
    {
      x_ = std::move(oth.x_);
      return *this;
    }

    const std::size_t * num_mindegs(std::size_t) const
    {
      return x_.get();
    }

    std::size_t * num_mindegs(std::size_t)
    {
      return x_.get();
    }

    const std::size_t * num_maxdegs(std::size_t nvars) const
    {
      return x_.get() + nvars;
    }

    std::size_t * num_maxdegs(std::size_t nvars)
    {
      return x_.get() + nvars;
    }

    const std::size_t * den_mindegs(std::size_t nvars) const
    {
      return x_.get() + 2*nvars;
    }

    std::size_t * den_mindegs(std::size_t nvars)
    {
      return x_.get() + 2*nvars;
    }

    const std::size_t * den_maxdegs(std::size_t nvars) const
    {
      return x_.get() + 3*nvars;
    }

    std::size_t * den_maxdegs(std::size_t nvars)
    {
      return x_.get() + 3*nvars;
    }

  private:
    std::unique_ptr<std::size_t[]> x_;
  };

} // namespace fflow
