#include <cctype>
#include <string>
#include <cstdlib>
#include <fflow/ratfun_parser.hh>
#include <fflow/rational_function.hh>
#include <fflow/mp_common.hh>
#include <fflow/mp_gcd.hh>

namespace fflow {
  namespace {

    struct RatFunToken {

      enum {
        COEFF,
        VAR,
        EXPONENT,
        TIMES,
        DIV,
        OPEN_PAR,
        CLOSE_PAR,
        END,
        INVALID
      };

      bool is_end_of_poly() const
      {
        return type == DIV || type == CLOSE_PAR || type == END;
      }

      unsigned type;
      const char * start;
      unsigned size;
      unsigned val;
    };

    struct RatFunTokenizer {

      static RatFunToken invalid()
      {
        return RatFunToken{RatFunToken::INVALID, nullptr, 0, 0};
      }

      static RatFunToken end_token()
      {
        return RatFunToken{RatFunToken::END, nullptr, 0, 0};
      }

      RatFunToken parse_var();
      RatFunToken parse_coeff();
      RatFunToken next_token();

      const char * cur;
      const char * end;

      const std::string * vars;
      unsigned nvars;
    };

    RatFunToken RatFunTokenizer::next_token()
    {
      if (cur == end)
        return end_token();

      while (std::isspace(*cur))
        ++cur;

      if (cur == end || *cur == '\0')
        return end_token();

      if (*cur == '+' || *cur == '-' || std::isdigit(*cur))
        return parse_coeff();

      if (*cur == '*') {
        RatFunToken tok{RatFunToken::TIMES, cur, 1, 0};
        ++cur;
        return tok;
      }

      if (*cur == '/') {
        RatFunToken tok{RatFunToken::DIV, cur, 1, 0};
        ++cur;
        return tok;
      }

      if (*cur == '(') {
        RatFunToken tok{RatFunToken::OPEN_PAR, cur, 1, 0};
        ++cur;
        return tok;
      }

      if (*cur == ')') {
        RatFunToken tok{RatFunToken::CLOSE_PAR, cur, 1, 0};
        ++cur;
        return tok;
      }

      if (std::isalpha(*cur) || *cur == '_')
        return parse_var();

      if (*cur == '^') {
        unsigned val;
        ++cur;
        const char * start = cur;
        val = std::strtol(cur, const_cast<char**>(&cur), 10);
        return RatFunToken{RatFunToken::EXPONENT,
            start, unsigned(cur-start), val};
      }

      return invalid();
    }


    RatFunToken RatFunTokenizer::parse_coeff()
    {
      const char * start = cur;
      char sign = 0;

      // optional sign
      if (*cur == '+' || *cur == '-') {
        if (*cur == '+')
          ++start;
        sign = *cur;
        ++cur;
      }

      while (std::isspace(*cur))
        ++cur;

      // numerator
      const char * startdigits = cur;
      while (std::isdigit(*cur))
        ++cur;

      // case where only a sign is present
      if (cur == startdigits) {
        if (sign == '+')
          return RatFunToken{RatFunToken::COEFF, start, unsigned(cur-start), 1};
        else if (sign == '-')
          return RatFunToken{RatFunToken::COEFF, start, unsigned(cur-start), 2};
        else
          return invalid();
      }

      while (std::isspace(*cur))
        ++cur;

      // optional denominator
      const char * maybe_denstart = cur;
      if (*cur == '/') {
        ++cur;
        while (std::isspace(*cur))
          ++cur;
        if (!std::isdigit(*cur))
          cur = maybe_denstart;
        while (std::isdigit(*cur))
          ++cur;
      }

      return RatFunToken{RatFunToken::COEFF, start, unsigned(cur-start), 0};
    }


    RatFunToken RatFunTokenizer::parse_var()
    {
      std::size_t max_length = end-cur;
      for (unsigned i=0; i<nvars; ++i) {
        if (vars[i].size() > max_length)
          continue;
        if (vars[i].compare(0, vars[i].size(), cur, vars[i].size()) == 0) {
          const char * start = cur;
          cur += vars[i].size();
          return RatFunToken{RatFunToken::VAR, start, unsigned(cur-start), i};
        }
      }
      return invalid();
    }


    struct RatFunParser {

      void next_token()
      {
        curr_token = tokenizer.next_token();
      }

      Ret parse_poly_term(std::vector<Monomial> & num,
                          std::vector<MPRational> & cnum);

      Ret parse_poly(std::vector<Monomial> & num,
                     std::vector<MPRational> & cnum);

      Ret parse_num_or_den(std::vector<Monomial> & num,
                           std::vector<MPRational> & cnum);

      void set_denom_to_one(std::vector<Monomial> & den,
                            std::vector<MPRational> & cden);

      Ret parse(const char * start, const char * end,
                std::vector<Monomial> & num, std::vector<MPRational> & cnum,
                std::vector<Monomial> & den, std::vector<MPRational> & cden);

#define DO_NEXT_TOKEN()                                 \
      do {                                              \
        next_token();                                   \
        if (curr_token.type == RatFunToken::INVALID)    \
          return FAILED;                                \
      } while (0)

      RatFunTokenizer tokenizer;
      RatFunToken curr_token;
      const std::string * vars;
      unsigned nvars;
    };

    Ret RatFunParser::parse(const char * start, const char * end,
                                std::vector<Monomial> & num,
                                std::vector<MPRational> & cnum,
                                std::vector<Monomial> & den,
                                std::vector<MPRational> & cden)
    {
      num.clear();
      cnum.clear();
      den.clear();
      cden.clear();

      tokenizer = RatFunTokenizer{start, end, vars, nvars};
      DO_NEXT_TOKEN();

      // parse numerator
      Ret ret = parse_num_or_den(num, cnum);
      if (ret == FAILED)
        return FAILED;

      // optionally parse denominator
      if (curr_token.type == RatFunToken::DIV) {
        DO_NEXT_TOKEN();
        ret = parse_num_or_den(den, cden);
        if (ret == FAILED)
          return FAILED;
      } else {
        set_denom_to_one(den, cden);
      }

      if (curr_token.type == RatFunToken::END)
        return SUCCESS;

      return FAILED;
    }

    Ret RatFunParser::parse_num_or_den(std::vector<Monomial> & num,
                                           std::vector<MPRational> & cnum)
    {
      if (curr_token.type == RatFunToken::OPEN_PAR) {
        DO_NEXT_TOKEN();
        parse_poly(num, cnum);
        if (curr_token.type == RatFunToken::CLOSE_PAR)
          DO_NEXT_TOKEN();
        else
          return FAILED;
      } else {
        return parse_poly(num, cnum);
      }
      return SUCCESS;
    }

    Ret RatFunParser::parse_poly(std::vector<Monomial> & num,
                                     std::vector<MPRational> & cnum)
    {
      do {
        Ret ret = parse_poly_term(num, cnum);
        if (ret == FAILED)
          return FAILED;
      } while (!curr_token.is_end_of_poly());
      return SUCCESS;
    }

    Ret RatFunParser::parse_poly_term(std::vector<Monomial> & num,
                                          std::vector<MPRational> & cnum)
    {
      bool got_coeff = false;
      std::string cstr;
      num.push_back(Monomial(nvars));
      Monomial & m = num.back();
      m.coeff() = 1;

      if (curr_token.type == RatFunToken::COEFF) {
        got_coeff = true;
        if (curr_token.val == 1) {
          cnum.push_back(MPRational(1));
        } else if (curr_token.val == 2) {
          cnum.push_back(MPRational(-1));
        } else {
          cstr.assign(curr_token.start, curr_token.size);
          cstr.erase(std::remove_if(cstr.begin(), cstr.end(), isspace),
                     cstr.end());
          cnum.push_back(MPRational(cstr.c_str()));
        }
        DO_NEXT_TOKEN();
      } else {
        cnum.push_back(MPRational(1));
      }

      bool first_var = true;
      unsigned tot_deg = 0;
      while (curr_token.type == RatFunToken::VAR ||
             curr_token.type == RatFunToken::TIMES) {
        if (!first_var || got_coeff)
          if (curr_token.type == RatFunToken::TIMES)
            DO_NEXT_TOKEN();
        unsigned cur_var = curr_token.val;
        DO_NEXT_TOKEN();
        if (curr_token.type == RatFunToken::EXPONENT) {
          m.exponent(cur_var) += curr_token.val;
          tot_deg += curr_token.val;
          DO_NEXT_TOKEN();
        } else {
          m.exponent(cur_var) += 1;
          tot_deg += 1;
        }
        first_var = false;
      }
      m.degree() = tot_deg;

      if (!got_coeff && first_var)
        return FAILED;

      // Check if there's an integer denominator, which we absorb in
      // the coefficient
      if (curr_token.type != RatFunToken::DIV)
        return SUCCESS;
      {
        while (std::isspace(*tokenizer.cur))
          ++tokenizer.cur;
        const char * cur = tokenizer.cur;
        if (std::isdigit(*cur)) {
          const char * cdenstart = cur++;
          while(std::isdigit(*cur))
            ++cur;
          cstr.assign(cdenstart, cur - cdenstart);
          MPRational cden(cstr.c_str());
          mpq_div(cnum.back().get(), cnum.back().get(), cden.get());
          tokenizer.cur = cur;
          DO_NEXT_TOKEN();
        }
      }

      return SUCCESS;
    }

    void RatFunParser::set_denom_to_one(std::vector<Monomial> & num,
                                        std::vector<MPRational> & cnum)
    {
      num.push_back(Monomial(nvars));
      Monomial & m = num.back();
      m.coeff() = 1;
      cnum.push_back(MPRational(1));
    }

  } // namespace


  Ret parse_rat_fun(const std::string vars[], unsigned nvars,
                    const char * start, const char * end,
                    std::vector<Monomial> & num,
                    std::vector<MPRational> & cnum,
                    std::vector<Monomial> & den,
                    std::vector<MPRational> & cden)
  {
    RatFunParser parser;
    parser.vars = vars;
    parser.nvars = nvars;
    return parser.parse(start, end, num, cnum, den, cden);
  }


  Ret parse_rat_fun_mod(const std::string vars[], unsigned nvars,
                        const char * start, const char * end,
                        UInt mod,
                        std::vector<Monomial> & num,
                        std::vector<Monomial> & den)
  {
    std::vector<MPRational> cnum, cden;
    Ret ret = parse_rat_fun(vars, nvars, start, end, num, cnum, den, cden);
    if (ret == FAILED)
      return FAILED;

    MPInt mmod(mod);
    MPInt cc;
    unsigned i=0;
    for (auto & coeff : cnum) {
      if (!rat_mod(coeff,mmod,cc))
        return FAILED;
      num[i].coeff() = cc.to_uint();
      ++i;
    }

    i=0;
    for (auto & coeff : cden) {
      if (!rat_mod(coeff,mmod,cc))
        return FAILED;
      den[i].coeff() = cc.to_uint();
      ++i;
    }

    return SUCCESS;
  }


} // namespace fflow
