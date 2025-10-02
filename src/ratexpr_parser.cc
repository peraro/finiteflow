#include <cctype>
#include <string>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <memory>
#include <fflow/ratexpr_parser.hh>
#include "spooky_v2.hh"

namespace fflow {

  namespace {

    const UInt HASH_SEED = 4680351337364691545ULL;

    struct HashVarPow {
      inline std::size_t operator()(AnalyticExpression::VarPow x) const
      {
        return SpookyHash::Hash64(&x, sizeof(x), HASH_SEED);
      }
    };

    struct EqVarPow {
      inline bool operator()(AnalyticExpression::VarPow v1,
                             AnalyticExpression::VarPow v2) const
      {
        return v1.var == v2.var && v1.exponent == v2.exponent;
      }
    };

    struct HashMPRational {
      inline std::size_t operator()(const MPRational & q) const
      {
        std::size_t hash = HASH_SEED;
        if (q.get()->_mp_num._mp_size < 0)
          hash += 134;
        hash = SpookyHash::Hash64(q.get()->_mp_den._mp_d,
                                  iabs(q.get()->_mp_den._mp_size),
                                  hash);
        hash = SpookyHash::Hash64(q.get()->_mp_num._mp_d,
                                  iabs(q.get()->_mp_num._mp_size),
                                  hash);
        return hash;
      }
    };

    struct EqVarMPRational {
      inline bool operator()(const MPRational & v1,
                             const MPRational & v2) const
      {
        return v1.cmp(v2) == 0;
      }
    };

    typedef std::unordered_map<AnalyticExpression::VarPow,
                               unsigned,HashVarPow, EqVarPow> VarPowMap;
    typedef std::unordered_map<MPRational,unsigned,
                               HashMPRational,EqVarMPRational> MPRationalMap;


    struct RatExprToken {

      enum {
        PLUS,
        TIMES,
        MINUS,
        DIV,
        POW,
        VAR,
        VARPOW,
        MACHINT,
        BIGINT,
        OPENPAR,
        CLOSEPAR,
        END,
        INVALID
      };

      const char * start;
      unsigned size;
      unsigned type;
      union {
        Int val;
        struct {
          unsigned var;
          int exponent;
        };
      };
    };

    struct RatExprTokenizer {

      typedef AnalyticExpression::VarPow VarPow;

      static RatExprToken end_token()
      {
        return RatExprToken{0,0,RatExprToken::END,{0}};
      }

      static RatExprToken invalid_token()
      {
        return RatExprToken{0,0,RatExprToken::INVALID,{0}};
      }

      void init();

      RatExprToken parse_var_or_varpow();
      RatExprToken parse_int();
      RatExprToken next_token();

      const char * cur;
      const char * end;
      std::string var_prefix;
      const std::string * vars;
      std::unique_ptr<unsigned[]> var_idxs;
      unsigned nvars;
      VarPowMap * varspow_map;
    };

    void RatExprTokenizer::init()
    {
      if (!vars)
        return;
      // sort variable indexes from longest to shortest name, so that
      // the tokenizer is "greedy" when matching variable names
      var_idxs.reset(new unsigned[nvars]);
      std::iota(var_idxs.get(), var_idxs.get()+nvars, 0);
      std::sort(var_idxs.get(), var_idxs.get()+nvars,
                [this](unsigned i, unsigned j) -> bool
                {
                  int cmp = compare(this->vars[i].size(),
                                    this->vars[j].size());
                  if (cmp)
                    return cmp > 0;
                  else
                    return i < j;
                });
    }

    RatExprToken RatExprTokenizer::next_token()
    {
      if (cur > end)
        return invalid_token();

      if (cur == end) {
        ++cur;
        return end_token();
      }

      while (std::isspace(*cur))
        ++cur;

      if (cur == end || *cur == '\0')
        return end_token();

      if (std::isdigit(*cur))
        return parse_int();

      if (*cur == '+')
        return RatExprToken{cur++,1,RatExprToken::PLUS,{0}};

      if (*cur == '-')
        return RatExprToken{cur++,1,RatExprToken::MINUS,{0}};

      if (*cur == '*')
        return RatExprToken{cur++,1,RatExprToken::TIMES,{0}};

      if (*cur == '/')
        return RatExprToken{cur++,1,RatExprToken::DIV,{0}};

      if (*cur == '(')
        return RatExprToken{cur++,1,RatExprToken::OPENPAR,{0}};

      if (*cur == ')')
        return RatExprToken{cur++,1,RatExprToken::CLOSEPAR,{0}};

      if (std::isalpha(*cur) || *cur == '_')
        return parse_var_or_varpow();

      if (*cur == '^')
        return RatExprToken{cur++,1,RatExprToken::POW,{0}};

      return invalid_token();
    }

    RatExprToken RatExprTokenizer::parse_var_or_varpow()
    {
      std::size_t max_length = end-cur;
      RatExprToken res = invalid_token();
      std::size_t pref_len = var_prefix.length();

      // Check variables
      if (pref_len && max_length >= pref_len+1
          && std::equal(var_prefix.c_str(),
                        var_prefix.c_str()+pref_len,
                        cur)
          && std::isdigit(cur[pref_len])) {
        const char * start = cur;
        unsigned val = cur[pref_len]-'0';
        cur = cur+pref_len+1;
        while (std::isdigit(*cur)) {
          val = 10*val + (*cur - '0');
          ++cur;
        }
        if (val < nvars) {
          res = RatExprToken{start, unsigned(cur-start),
                             RatExprToken::VAR, {val}};
        } else {
          logerr("Variable index out of bounds");
          return invalid_token();
        }
      }

      if (res.type == RatExprToken::INVALID && vars) {
        for (unsigned j=0; j<nvars; ++j) {
          unsigned i = var_idxs[j];
          if (vars[i].size() > max_length)
            continue;
          if (vars[i].compare(0, vars[i].size(), cur, vars[i].size()) == 0) {
            const char * start = cur;
            cur += vars[i].size();
            res = RatExprToken{start, unsigned(cur-start),
                               RatExprToken::VAR, {i}};
            break;
          }
        }
        if (res.type == RatExprToken::INVALID)
          return res;
      }

      while (std::isspace(*cur))
        ++cur;

      if (*cur == '^') {
        // For powers of variables we parse the exponent directly here
        // and consider var^exponent as a single token
        ++cur;
        while (std::isspace(*cur))
          ++cur;
        bool open_par = false;
        if (*cur == '(') {
          open_par = true;
          ++cur;
          while (std::isspace(*cur))
            ++cur;
        }
        int sign = 1;
        if (*cur == '-') {
          sign = -1;
          ++cur;
        }
        if (!std::isdigit(*cur)) {
          logerr("Non-numeric exponent found");
          return invalid_token();
        }
        int exponent = *cur-'0';
        ++cur;
        while (std::isdigit(*cur)) {
          exponent = 10*exponent + (*cur - '0');
          ++cur;
        }
        if (open_par) {
          while (std::isspace(*cur))
            ++cur;
          if (*cur != ')')
            return invalid_token();
          ++cur;
        }
        exponent *= sign;
        (*varspow_map)[VarPow{static_cast<unsigned int>(res.val),exponent}] = 0;
        res = RatExprToken{res.start, unsigned(cur-res.start),
                           RatExprToken::VARPOW, {res.val}};
        res.var = res.val;
        res.exponent = exponent;
      }

      return res;
    }

    RatExprToken RatExprTokenizer::parse_int()
    {
      // The abs of any integer up to 18 digits is smaller than 2^62,
      // which is smaller than our smallest prime, hence we can just
      // save it as a machine integer (with no mod required).  From 19
      // digits and above, we use multiprecision numbers.

      const char * start = cur;
      ++cur;

      while (std::isdigit(*cur))
        ++cur;

      const unsigned ndigits = cur - start;

      if (ndigits <= 18) {
        Int val = start[0]-'0';
        for (unsigned j=1; j<ndigits; ++j)
          val = val*10 + (start[j]-'0');
        return RatExprToken{start, ndigits, RatExprToken::MACHINT, {val}};
      }

      return RatExprToken{start, ndigits, RatExprToken::BIGINT, {0}};
    }


    struct RatExprParser {

      typedef AnalyticExpression::VarPow VarPow;

      MPRationalMap bignums_map;

      const RatExprToken * tokens_beg;
      VarPowMap * varspow_map;
      std::vector<Instruction> * bytecode;
      std::vector<MPRational> * bignums;

      const RatExprToken * parse();

      const RatExprToken * parse_expr(const RatExprToken * start);
      const RatExprToken * parse_term(const RatExprToken * start, int sign);
      const RatExprToken * parse_factor(const RatExprToken * start,
                                        int expsign, MPRational & coeff,
                                        std::size_t & n_factors);
      const RatExprToken * parse_exponent(const RatExprToken * start,
                                          int expsign);
      const RatExprToken * get_exponent(const RatExprToken * start,
                                        Int & expval);
      void push_exponent(Int expval);

      void push_machinenum(Int num);
      void push_bignum(const MPRational & q);

    };

    void RatExprParser::push_machinenum(Int num)
    {
      int sign = num >= 0 ? 1 : -1;
      UInt absnum = num*sign;

      if (absnum < 256) {
        (*bytecode).push_back(AnalyticExpression::SMALLNUM);
        (*bytecode).push_back(absnum);
      } else {
        (*bytecode).push_back(AnalyticExpression::MEDNUM);
        bool init_num = false;
        for (int j=sizeof(UInt)/sizeof(Instruction)-1; j>=0; --j) {
          Instruction instr = (absnum >> j*8); // mod 8 implicit in conversion
          if (!init_num && instr) {
            (*bytecode).push_back(j+1);
            init_num = true;
          }
          if (init_num)
            (*bytecode).push_back(instr);
        }
      }

      if (sign < 0)
        (*bytecode).push_back(AnalyticExpression::NEG);
    }

    void RatExprParser::push_bignum(const MPRational & q)
    {
      auto found = bignums_map.find(q);
      if (found != bignums_map.end()) {
        push_machinenum(found->second);
        bytecode->push_back(AnalyticExpression::BIGNUM);
        return;
      }

      unsigned idx = bignums_map.size();
      bignums_map[q] = idx;
      push_machinenum(idx);
      bytecode->push_back(AnalyticExpression::BIGNUM);
      return;
    }

    const RatExprToken * RatExprParser::parse()
    {
      const RatExprToken * tk = parse_expr(tokens_beg);
      if (!tk)
        return 0;

      if (tk->type != RatExprToken::END) {
        logerr("Unexpected token");
        return 0;
      }

      (*bytecode).push_back(AnalyticExpression::END);

      return tk+1;
    }

    const RatExprToken * RatExprParser::parse_expr(const RatExprToken * start)
    {
      const RatExprToken * tk = parse_term(start,1);
      std::size_t n_terms = 1;
      while (true) {
        if (!tk)
          return 0;
        if (tk->type == RatExprToken::PLUS)
          tk = parse_term(tk+1, 1);
        else if (tk->type == RatExprToken::MINUS)
          tk = parse_term(tk+1, -1);
        else
          break;
        ++n_terms;
      }
      if (n_terms > 1) {
        push_machinenum(n_terms);
        bytecode->push_back(AnalyticExpression::ADD);
      }
      return tk;
    }

    static bool implicit_times_allowed(const RatExprToken * tk)
    {
      // range of tokens than can be the start of a factor after an implicit '*'
      return tk->type >= RatExprToken::VAR && tk->type <= RatExprToken::OPENPAR;
    }

    const RatExprToken * RatExprParser::parse_term(const RatExprToken * start,
                                                   int sign)
    {
      const Int BIG_NUM_POS_THRESH = UInt(1) << 62;
      const Int BIG_NUM_NEG_THRESH = -Int(UInt(1) << 62);

      MPRational coeff(sign);
      std::size_t n_factors = 0;
      const RatExprToken * tk = parse_factor(start,1,coeff,n_factors);
      while (true) {
        if (!tk)
          return 0;
        if (tk->type == RatExprToken::TIMES)
          tk = parse_factor(tk+1, 1, coeff, n_factors);
        else if (tk->type == RatExprToken::DIV)
          tk = parse_factor(tk+1, -1, coeff, n_factors);
        else if (implicit_times_allowed(tk))
          tk = parse_factor(tk, 1, coeff, n_factors);
        else
          break;
      }


      // Handle numerical coefficient

      bool to_neg = false;

      if (coeff.cmp(1) == 0) {

        if (!n_factors) {
          push_machinenum(1);
          ++n_factors;
        }

      } else if (coeff.cmp(-1) == 0) {

        if (!n_factors) {
          push_machinenum(-1);
          ++n_factors;
        } else {
          to_neg = true;
        }

      } else {

        if (mpz_cmp_ui(mpq_denref(coeff.get()),1) != 0 // has denominator
            || coeff.cmp(BIG_NUM_NEG_THRESH) < 0
            || coeff.cmp(BIG_NUM_POS_THRESH) > 0) {

          // big rational number
          push_bignum(coeff);
          ++n_factors;

        } else {

          // machine number
          Int num, den;
          coeff.to_int(num, den);
          // den==1 is guaranteed here
          push_machinenum(num);
          ++n_factors;

        }

      }

      if (n_factors > 1) {
        push_machinenum(n_factors);
        bytecode->push_back(AnalyticExpression::MUL);
      }
      if (to_neg)
        bytecode->push_back(AnalyticExpression::NEG);
      return tk;
    }

    const RatExprToken * RatExprParser::parse_factor(const RatExprToken * start,
                                                     int exponentsign,
                                                     MPRational & coeff,
                                                     std::size_t & n_factors)
    {
      const RatExprToken * tk = start;
      const unsigned base_type = tk->type;

      if (base_type == RatExprToken::VAR) {
        push_machinenum(tk->val);
        bytecode->push_back(AnalyticExpression::VAR);
        ++n_factors;
        return tk+1;
      }

      if (base_type == RatExprToken::VARPOW) {
        auto key = VarPow{tk->var,tk->exponent};
        unsigned idx = (*varspow_map)[key];
        push_machinenum(idx);
        bytecode->push_back(AnalyticExpression::VARPOW);
        ++n_factors;
        return tk+1;
      }

      if (base_type == RatExprToken::MACHINT) {

        Int val = tk->val;
        Int exponent;
        tk = get_exponent(tk+1, exponent);
        if (!tk)
          return 0;
        exponent *= exponentsign;

        if (exponent == 1) {

          MPRational qval(val);
          mpq_mul(coeff.get(), coeff.get(), qval.get());
          return tk;

        } else if (exponent == -1) {

          MPRational qval(val);
          mpq_div(coeff.get(), coeff.get(), qval.get());
          return tk;

        } else {

          // Integers raised to non-trivial exponents are inserted as
          // additional factors and exponentiation is handled during
          // evaluation.
          push_machinenum(val);
          push_exponent(exponent);
          ++n_factors;
          return tk;

        }
      }

      if (base_type == RatExprToken::OPENPAR) {
        tk = parse_expr(tk+1);
        if (!tk)
          return 0;
        if (tk->type != RatExprToken::CLOSEPAR) {
          logerr("Missing ')'");
          return 0;
        }
        ++n_factors;
        return parse_exponent(tk+1, exponentsign);
      }


      if (base_type == RatExprToken::BIGINT) {

        std::string qstr(tk->start, tk->size);
        MPRational qval(qstr.c_str());

        // logic similar to MACHINT
        Int exponent;
        tk = get_exponent(tk+1, exponent);
        if (!tk)
          return 0;
        exponent *= exponentsign;

        if (exponent == 1) {

          mpq_mul(coeff.get(), coeff.get(), qval.get());
          return tk;

        } else if (exponent == -1) {

          mpq_div(coeff.get(), coeff.get(), qval.get());
          return tk;

        } else {

          // Non-trivial exponents handled as additional factors
          push_bignum(qval);
          push_exponent(exponent);
          ++n_factors;
          return tk;

        }

        return 0;
      }

      // handle prefix operators

      if (base_type == RatExprToken::PLUS)
        return tk+1;

      if (base_type == RatExprToken::MINUS) {
        mpq_neg(coeff.get(), coeff.get());
        return tk+1;
      }

      logerr(format("Unexpected token: {}",base_type));

      return 0;
    }

    const RatExprToken *
    RatExprParser::parse_exponent(const RatExprToken * start,
                                  int exponentsign)
    {
      Int exponent;
      const RatExprToken * tk = get_exponent(start, exponent);
      if (!tk)
        return 0;

      push_exponent(exponent*exponentsign);

      return tk;
    }

    void RatExprParser::push_exponent(Int exponent)
    {
      if (exponent == 1)
        return;
      if (exponent >= 0) {
        push_machinenum(exponent);
        bytecode->push_back(AnalyticExpression::POW);
      } else {
        push_machinenum(-exponent);
        bytecode->push_back(AnalyticExpression::NEGPOW);
      }
    }


    const RatExprToken *
    RatExprParser::get_exponent(const RatExprToken * start, Int & expval)
    {
      if (start->type != RatExprToken::POW) {
        expval = 1;
        return start;
      }

      const RatExprToken * tk = start + 1;

      bool open_par = false;
      if (tk->type == RatExprToken::OPENPAR) {
        open_par = true;
        ++tk;
      }

      int sign = 1;
      if (tk->type == RatExprToken::MINUS) {
        sign = -1;
        ++tk;
      }

      if (tk->type != RatExprToken::MACHINT) {
        logerr("Invalid exponent");
        return 0;
      }
      expval = tk->val*sign;
      ++tk;

      if (open_par) {
        if (tk->type != RatExprToken::CLOSEPAR)
          return 0;
        ++tk;
      }

      return tk;
    }


  } // namespace


  Ret parse_ratexpr_list(unsigned npars,
                         const std::string & vars_prefix,
                         const std::string vars[],
                         const ExprCStr inputs[], const unsigned len[],
                         unsigned nfunctions,
                         std::vector<std::vector<Instruction>> & bytecode,
                         std::vector<MPRational> & bignums,
                         std::vector<AnalyticExpression::VarPow> & varpows)
  {
    if (!nfunctions) {
      bytecode.clear();
      bignums.clear();
      varpows.clear();
      return SUCCESS;
    }

    RatExprTokenizer tokenizer;
    tokenizer.var_prefix = vars_prefix;
    tokenizer.vars = vars;
    tokenizer.nvars = npars;
    VarPowMap varspow_map;
    tokenizer.varspow_map = &varspow_map;
    tokenizer.init();
    std::vector<RatExprToken> tokens;
    for (unsigned j=0; j<nfunctions; ++j) {
      tokenizer.cur = inputs[j];
      tokenizer.end = inputs[j] + len[j];
      while (true) {
        RatExprToken tk = tokenizer.next_token();
        if (tk.type == RatExprToken::INVALID)
          return FAILED;
        tokens.push_back(tk);
        if (tk.type == RatExprToken::END)
          break;
      }
    }
    tokens.push_back(RatExprTokenizer::invalid_token());

    varpows.clear();
    varpows.reserve(varspow_map.size());
    for (const auto & el : varspow_map)
      varpows.push_back(el.first);
    std::sort(varpows.begin(), varpows.end(),
              [](AnalyticExpression::VarPow v1,
                 AnalyticExpression::VarPow v2) -> bool
              {
                int cmp = compare(v1.var,v2.var);
                if (cmp)
                  return cmp < 0;
                else
                  return v1.exponent < v2.exponent;
              });
    unsigned idx=0;
    for (const auto el : varpows)
      varspow_map[el] = idx++;

    RatExprParser parser;
    parser.varspow_map = & varspow_map;
    parser.tokens_beg = &tokens[0];
    bytecode.clear();
    bytecode.resize(nfunctions);
    for (unsigned j=0; j<nfunctions; ++j) {
      parser.bytecode = &bytecode[j];
      const RatExprToken * ret = parser.parse();
      if (!ret)
        return FAILED;
      parser.tokens_beg = ret;
    }

    bignums.clear();
    bignums.resize(parser.bignums_map.size());
    for (const auto & el : parser.bignums_map)
      bignums[el.second] = el.first;

    return SUCCESS;
  }

} // namespace fflow
