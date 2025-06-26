#include <fstream>
#include <numeric>
#include <unordered_map>
#include <fflow/graph.hh>
#include <fflow/analytic_solver.hh>
#include <fflow/alg_functions.hh>
#include <fflow/json.hh>
#include "ffuhash.h"

namespace {

  struct HashFFUHash {
    inline std::size_t operator()(const FFUHashVal & val) const
    {
      return ((std::size_t)val.val[0] << 32) | val.val[1];
    }
  };

  struct EqFFUHash {
    inline bool operator()(const FFUHashVal & v1, const FFUHashVal & v2) const
    {
      return (v1.val[0] == v2.val[0])
        && (v1.val[1] == v2.val[1])
        && (v1.val[2] == v2.val[2])
        && (v1.val[3] == v2.val[3])
        && (v1.val[4] == v2.val[4]);
    }
  };

}

namespace fflow {

  std::string read_file(const char * filename)
  {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (file) {
      std::string ret;
      file.seekg(0, std::ios::end);
      ret.resize(file.tellg());
      file.seekg(0, std::ios::beg);
      file.read(&ret[0], ret.size());
      file.close();
      return ret;
    }
    return "";
  }

  // Right now it only handles "\/" -> "/" and "\\" -> "\"
  void unescape_string(std::string & s)
  {
    const char * read = s.c_str();
    const char * end = read + s.size();
    auto write = s.begin();

    while (read+1 < end) {
      if (*read == '\\') {
        if (read[1] == '/')
          *write = '/';
        if (read[1] == '\\')
          *write = '\\';
        read += 2;
        ++write;
      } else {
        *write = *read;
        ++read;
        ++write;
      }
    }
    if (read < end) {
      *write = *read;
      ++write;
      ++read;
    }

    s.resize(write - s.begin());
  }

  namespace  {

    struct JSONToken {

      enum {
        OPEN_SQUARE_PAR,
        CLOSE_SQUARE_PAR,
        COMMA,
        STRING,
        INTEGER,
        END,
        INVALID
      };

      bool is_invalid() const
      {
        return type == INVALID;
      }

      unsigned type;
      unsigned size;
      const char * start;
      UInt val;
    };

    struct JSONTokenizer {

      static JSONToken end_token()
      {
        return JSONToken{JSONToken::END, 0, nullptr, 0};
      }

      static JSONToken invalid()
      {
        return JSONToken{JSONToken::INVALID, 0, nullptr, 0};
      }

      JSONToken next_token();

      const char * cur, * end;
    };


    JSONToken JSONTokenizer::next_token()
    {
      if (cur == end)
        return end_token();

      while (std::isspace(*cur))
        ++cur;

      if (cur == end || *cur == '\0')
        return end_token();

      if (*cur == '[') {
        auto ret = JSONToken{JSONToken::OPEN_SQUARE_PAR,1,cur,0};
        ++cur;
        return ret;
      }

      if (*cur == ']') {
        auto ret = JSONToken{JSONToken::CLOSE_SQUARE_PAR,1,cur,0};
        ++cur;
        return ret;
      }

      if (*cur == ',') {
        auto ret = JSONToken{JSONToken::COMMA,1,cur,0};
        ++cur;
        return ret;
      }

      if (*cur == '"') {
        ++cur;
        const char * start = cur;
        while (cur != end && *cur != '"')
          ++cur;
        if (cur == end)
          return invalid();
        auto ret = JSONToken{JSONToken::STRING, unsigned(cur-start), start, 0};
        ++cur;
        return ret;
      }

      if (std::isdigit(*cur)) {
        const char * start = cur;
        UInt val = std::strtol(cur, const_cast<char**>(&cur), 10);
        auto ret =JSONToken{JSONToken::INTEGER,
                            unsigned(cur-start), start, val};
        return ret;
      }

      return invalid();
    }


    struct JSONParser {

      void next_token()
      {
        curr_token = tokenizer.next_token();
      }

      Ret parse_sparse_system(const char * filename,
                              AnalyticSparseSolver & sys,
                              AnalyticSparseSolverData & data,
                              bool has_info = false);
      Ret parse_sparse_system_data(const char * filename,
                                   AnalyticSparseSolver & sys,
                                   AnalyticSparseSolverData & data,
                                   bool has_info);
      Ret skip_list();
      Ret parse_eq(unsigned i,
                   AnalyticSparseSolver & sys,
                   AnalyticSparseSolverData & data);
      Ret check_ls_ratfun(std::vector<HornerRatFunPtr> & f,
                          std::vector<MPHornerRatFunMap> & map,
                          std::size_t & idx);
      Ret parse_ratfun(HornerRatFunPtr & f, MPHornerRatFunMap & map);
      Ret parse_ratfun_list(const char * filename,
                            AnalyticFunction & fun,
                            AnalyticFunctionData & data);
      Ret parse_integer_list(const char * filename,
                             std::vector<unsigned> & out);
      Ret parse_poly(Monomial * monomials, MPRational* coeffs,
                     unsigned n_terms);

      void resize_monomials(std::vector<Monomial> & m, std::size_t n)
      {
        if (m.size() < n) {
          std::size_t curr_size = m.size();
          m.resize(n);
          for (unsigned j=curr_size; j<n; ++j) {
            m[j] = Monomial(npars);
            m[j].coeff() = 1;
          }
        }
      }

#define DO_NEXT_TOKEN()                                         \
      do {                                                      \
        next_token();                                           \
        if (FF_ERRCOND(curr_token.type == JSONToken::INVALID))  \
          return FAILED;                                        \
      } while (0)
#define EXPECT_TOKEN(etype)                       \
      do {                                        \
        next_token();                             \
        if (FF_ERRCOND(curr_token.type != etype)) \
          return FAILED;                          \
      } while (0)
#define EXPECT_TOKEN_EX(etype,logmsg)               \
      do {                                          \
        next_token();                               \
        if (FF_ERRCOND(curr_token.type != etype)) { \
          logerr(logmsg);                           \
          return FAILED;                            \
        }                                           \
      } while (0)

      JSONTokenizer tokenizer;
      JSONToken curr_token;
      unsigned neqs, nvars, npars, this_eq;
      std::size_t needed_workspace_size = 0;
      std::vector<unsigned> needed_vars;
      std::vector<std::string> files;
      std::string str_buff;

      std::vector<Monomial> num_monomials, den_monomials;

      std::unique_ptr<bool[]> eq_is_needed;

      std::unordered_map<FFUHashVal, std::size_t, HashFFUHash, EqFFUHash> cmap;
    };


    // The info on the system should have the form
    //
    // [neqs, nvars, npars, n_needed, neededvars, nfiles, files]
    //
    // where everything is an integer except the entries:
    //
    //   neededvars = [var1, var2, ...] | 0
    //   files = ["file1", "file2", ...]
    Ret JSONParser::parse_sparse_system(const char * filename,
                                        AnalyticSparseSolver & sys,
                                        AnalyticSparseSolverData & data,
                                        bool has_info)
    {
      std::string json_data = read_file(filename);
      if (json_data.size() == 0) {
        logerr(format("JSON: file '{}' not found or empty", filename));
        return FAILED;
      }

      tokenizer = JSONTokenizer{json_data.c_str(),
                                json_data.c_str() + json_data.size()};

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);

      EXPECT_TOKEN(JSONToken::INTEGER);
      neqs = curr_token.val;

      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::INTEGER);
      nvars = curr_token.val;

      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::INTEGER);
      npars = curr_token.val;

      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::INTEGER);
      needed_vars.resize(curr_token.val);

      EXPECT_TOKEN(JSONToken::COMMA);

      DO_NEXT_TOKEN();
      if (curr_token.type == JSONToken::INTEGER && curr_token.val == 0) {
        if (needed_vars.size() != nvars)
          return FAILED;
        std::iota(needed_vars.begin(), needed_vars.end(), 0);
      } else if (curr_token.type == JSONToken::OPEN_SQUARE_PAR) {
        unsigned n_needed = needed_vars.size();
        for (unsigned j=0; j<n_needed; ++j) {
          if (j)
            EXPECT_TOKEN(JSONToken::COMMA);
          EXPECT_TOKEN(JSONToken::INTEGER);
          needed_vars[j] = curr_token.val;
        }
        EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);
      } else {
        return FAILED;
      }

      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::INTEGER);
      unsigned n_files = curr_token.val;
      files.resize(n_files);

      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);

      for (unsigned j=0; j<n_files; ++j) {
        if (j)
          EXPECT_TOKEN(JSONToken::COMMA);
        EXPECT_TOKEN(JSONToken::STRING);
        files[j].assign(curr_token.start, curr_token.size);
        unescape_string(files[j]);
      }

      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);
      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      sys.rinfo.resize(neqs);
      data.c.resize(neqs);
      sys.cmap.resize(neqs);
      if (!has_info) {
        sys.nparsin.resize(1);
        sys.nparsin[0] = npars;
        sys.init(neqs, nvars, needed_vars.data(), needed_vars.size(), data);
      }

      if (has_info) {
        eq_is_needed.reset(new bool[neqs]());
        const unsigned * needed_eqs = sys.indep_eqs();
        std::size_t n_needed_eqs = sys.n_indep_eqs();
        for (unsigned j=0; j<n_needed_eqs; ++j)
          eq_is_needed[needed_eqs[j]] = true;
      }

      needed_workspace_size = 0;
      this_eq = 0;
      for (auto & f : files) {
        Ret ret = parse_sparse_system_data(f.c_str(), sys, data, has_info);
        if (ret != SUCCESS)
          return FAILED;
      }

      return SUCCESS;
    }


    // Format is:
    //
    //   [neqs_in_file, [col1, col2, ...], [coeff1, coeff2, ...]]
    //
    // where the first sublist are integers, and the second are
    // rational functions
    Ret JSONParser::parse_sparse_system_data(const char * filename,
                                             AnalyticSparseSolver & sys,
                                             AnalyticSparseSolverData & data,
                                             bool has_info)
    {
      std::string json_data = read_file(filename);
      if (json_data.size() == 0) {
        logerr(format("JSON: file '{}' not found or empty", filename));
        return FAILED;
      }

      tokenizer = JSONTokenizer{json_data.c_str(),
                                json_data.c_str() + json_data.size()};

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
      EXPECT_TOKEN(JSONToken::INTEGER);
      unsigned eq_in_file = curr_token.val;

      if (this_eq + eq_in_file > neqs)
        return FAILED;

      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
      for (unsigned j=0; j<eq_in_file; ++j) {

        if (j)
          EXPECT_TOKEN(JSONToken::COMMA);

        Ret ret;
        if (has_info && !eq_is_needed[this_eq])
          ret = skip_list();
        else
          ret = parse_eq(this_eq, sys, data);
        if (ret != SUCCESS)
          return FAILED;

        ++this_eq;
      }
      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      return SUCCESS;
    }


    Ret JSONParser::skip_list()
    {
      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
      int balanced_pars = 1;
      while (1) {
        DO_NEXT_TOKEN();
        if (curr_token.type == JSONToken::OPEN_SQUARE_PAR)
          balanced_pars += 1;
        else if (curr_token.type == JSONToken::CLOSE_SQUARE_PAR)
          balanced_pars -= 1;
        if (balanced_pars == 0)
          return SUCCESS;
      };
      return FAILED;
    }


    Ret JSONParser::parse_eq(unsigned i,
                             AnalyticSparseSolver & sys,
                             AnalyticSparseSolverData & data)
    {
      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);

      EXPECT_TOKEN(JSONToken::INTEGER);
      AnalyticSparseSolver::RowInfo & rinf = sys.rinfo[i];
      unsigned csize = curr_token.val;
      rinf.size = csize;
      rinf.cols.reset(new unsigned[csize]);
      rinf.idx.reset(new std::size_t[csize]);

      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
      unsigned * rinfptr = rinf.cols.get();
      for (unsigned j=0; j<csize; ++j) {
        if (j)
          EXPECT_TOKEN(JSONToken::COMMA);
        EXPECT_TOKEN(JSONToken::INTEGER);
        if (curr_token.val > nvars+1) {
          logerr("JSON: Column index of sparse equation is out of bounds");
          return FAILED;
        }
        rinfptr[j] = curr_token.val;
      }
      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
      for (unsigned j=0; j<csize; ++j) {
        if (j)
          EXPECT_TOKEN(JSONToken::COMMA);
        Ret ret = check_ls_ratfun(data.c, sys.cmap, rinf.idx[j]);
        if (ret != SUCCESS)
            return ret;
      }
      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      return SUCCESS;
    }


    // Parse rational function for a Sparse Solver.
    //
    // We try to quickly get to the end of the substring for the
    // function.  If check with our hashmap whether (the "unique" hash
    // of) the same substring has already been parsed and in that case
    // we add it to the list of coefficients.  Otherwise we proceed
    // via normal parsing.
    Ret JSONParser::check_ls_ratfun(std::vector<HornerRatFunPtr> & f,
                                    std::vector<MPHornerRatFunMap> & map,
                                    std::size_t & idx)
    {
      const char * c = tokenizer.cur;
      const char * cend = tokenizer.end;
      while (isspace(*c) && c < cend)
        ++c;

      // NOTE: zero functions are not allowed in sparse solvers, so we
      // don't check for that
      if (*c != '[')
        return FAILED;

      const char * c2 = c+1;
      int par_count = 1;
      for (;c2 < cend; ++c2) {
        if (*c2 == '[') {
          ++par_count;
        } else if (*c2 == ']') {
          --par_count;
          if (par_count == 0)
            break;
        }
      }
      ++c2;

      if (par_count != 0)
        return FAILED;

      FFUHashVal hval = ffUHash(static_cast<const void *>(c), c2-c);
      auto found = cmap.find(hval);
      if (found != cmap.end()) {
        idx = found->second;
        tokenizer.cur = c2;
        return SUCCESS;
      } else {
        // insert a new element at the end
        idx = map.size();
        cmap[hval] = idx;
        f.push_back(HornerRatFunPtr());
        map.push_back(MPHornerRatFunMap());
        auto & fel = f.back();
        auto & mapel = map.back();
        return parse_ratfun(fel, mapel);
      }

      return SUCCESS;
    }


    // List of rational functions:
    //
    //   [nfunctions, nvars, [func1, func2, ...]]
    //
    Ret JSONParser::parse_ratfun_list(const char * filename,
                                      AnalyticFunction & fun,
                                      AnalyticFunctionData & data)
    {
      std::string json_data = read_file(filename);
      if (json_data.size() == 0) {
        logerr(format("JSON: file '{}' not found or empty", filename));
        return FAILED;
      }

      tokenizer = JSONTokenizer{json_data.c_str(),
                                json_data.c_str() + json_data.size()};

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);

      EXPECT_TOKEN(JSONToken::INTEGER);
      unsigned nfunctions = curr_token.val;
      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::INTEGER);
      unsigned nvars = curr_token.val;
      EXPECT_TOKEN(JSONToken::COMMA);

      npars = nvars;
      data.f.reset(new HornerRatFunPtr[nfunctions]);
      fun.fmap.reset(new MPHornerRatFunMap[nfunctions]);

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
      for (unsigned j=0; j<nfunctions; ++j) {
        if (j)
          EXPECT_TOKEN(JSONToken::COMMA);
        Ret ret = parse_ratfun(data.f[j], fun.fmap[j]);
        if (ret != SUCCESS)
          return ret;
      }
      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      fun.init(nvars, nfunctions, data);

      return SUCCESS;
    }


    // A rational function has the form:
    //
    //   [[n_num_terms, monomial_num_list], [n_den_terms, monomial_den_list]]
    //
    // Alternatively, it can be just the number 0 which stands for a
    // vanishing function.
    Ret JSONParser::parse_ratfun(HornerRatFunPtr & f, MPHornerRatFunMap & map)
    {
      DO_NEXT_TOKEN();

      if (curr_token.type == JSONToken::INTEGER && curr_token.val == 0)
        return SUCCESS;

      if (curr_token.type != JSONToken::OPEN_SQUARE_PAR)
        return FAILED;

      unsigned num_terms, den_terms;

      // num
      {
        EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
        EXPECT_TOKEN(JSONToken::INTEGER);
        num_terms = curr_token.val;
        map.num_map.resize(num_terms);
        resize_monomials(num_monomials, num_terms);
        EXPECT_TOKEN(JSONToken::COMMA);
        Ret ret = parse_poly(num_monomials.data(), map.num_map.coeff.get(),
                             num_terms);
        if (ret != SUCCESS)
          return FAILED;
        EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);
      }

      EXPECT_TOKEN(JSONToken::COMMA);

      // den
      {
        EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
        EXPECT_TOKEN(JSONToken::INTEGER);
        den_terms = curr_token.val;
        map.den_map.resize(den_terms);
        resize_monomials(den_monomials, den_terms);
        EXPECT_TOKEN(JSONToken::COMMA);
        Ret ret = parse_poly(den_monomials.data(), map.den_map.coeff.get(),
                             den_terms);
        if (ret != SUCCESS)
          return FAILED;
        EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);
      }

      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      f.from_sparse_poly(num_monomials.data(), num_terms,
                         den_monomials.data(), den_terms,
                         npars, 0,
                         map.num_map.pos.get(), map.den_map.pos.get());

      std::size_t wspace = std::max(horner_required_workspace(f.num()),
                                    horner_required_workspace(f.den()));
      needed_workspace_size = std::max(wspace, needed_workspace_size);

      return SUCCESS;
    }


    // A list of monomials in the form
    //
    //  [["rational_coeff", [exponent1, exponent2,...]], ...]
    Ret JSONParser::parse_poly(Monomial * monomials, MPRational* coeffs,
                               unsigned n_terms)
    {
      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);

      for (unsigned j=0; j<n_terms; ++j) {
        if (j)
          EXPECT_TOKEN(JSONToken::COMMA);

        EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
        EXPECT_TOKEN(JSONToken::STRING);
        str_buff.assign(curr_token.start, curr_token.size);
        unescape_string(str_buff);
        coeffs[j].set(str_buff.c_str());
        Monomial & m = monomials[j];
        unsigned tot_deg = 0;
        EXPECT_TOKEN(JSONToken::COMMA);
        EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
        for (unsigned ex=0; ex<npars; ++ex) {
          if (ex)
            EXPECT_TOKEN_EX(JSONToken::COMMA,
                            "JSON: unexpected number of free parameters");
          EXPECT_TOKEN(JSONToken::INTEGER);
          m.exponent(ex) = curr_token.val;
          tot_deg += curr_token.val;
        }
        m.degree() = tot_deg;
        EXPECT_TOKEN_EX(JSONToken::CLOSE_SQUARE_PAR,
                        "JSON: unexpected number of free parameters");
        EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);
      }

      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      return SUCCESS;
    }


    // Format is:
    //
    //   [listlen, [i1,i2,i3,...]]
    //
    Ret JSONParser::parse_integer_list(const char * filename,
                                       std::vector<unsigned> & out)
    {
      std::string json_data = read_file(filename);
      if (json_data.size() == 0) {
        logerr(format("JSON: file '{}' not found or empty", filename));
        return FAILED;
      }

      tokenizer = JSONTokenizer{json_data.c_str(),
                                json_data.c_str() + json_data.size()};

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
      EXPECT_TOKEN(JSONToken::INTEGER);
      out.resize(curr_token.val);
      EXPECT_TOKEN(JSONToken::COMMA);

      EXPECT_TOKEN(JSONToken::OPEN_SQUARE_PAR);
      bool first = true;
      for (unsigned & i : out) {
        if (!first)
          EXPECT_TOKEN(JSONToken::COMMA);
        else
          first = false;
        EXPECT_TOKEN(JSONToken::INTEGER);
        i = curr_token.val;
      }
      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      EXPECT_TOKEN(JSONToken::CLOSE_SQUARE_PAR);

      return SUCCESS;
    }

  } // namespace

  Ret json_sparse_system(const char * info_file,
                         AnalyticSparseSolver & sys,
                         AnalyticSparseSolverData & data,
                         unsigned & required_workspace)
  {
    JSONParser parser;
    Ret ret = parser.parse_sparse_system(info_file, sys, data);
    required_workspace = parser.needed_workspace_size;
    return ret;
  }


  Ret json_sparse_system_with_info(const char * json_info_file,
                                   const char * binary_info_file,
                                   AnalyticSparseSolver & sys,
                                   AnalyticSparseSolverData & data,
                                   unsigned & required_workspace)
  {
    JSONParser parser;
    Ret ret = sys.load_info(&data, binary_info_file);
    if (ret == FAILED)
      return FAILED;
    ret = parser.parse_sparse_system(json_info_file, sys, data, true);
    required_workspace = parser.needed_workspace_size;
    return ret;
  }


  Ret json_sparse_ratfun(const char * file,
                         AnalyticFunction & fun,
                         AnalyticFunctionData & data,
                         unsigned & required_workspace)
  {
    JSONParser parser;
    Ret ret = parser.parse_ratfun_list(file, fun, data);
    required_workspace = parser.needed_workspace_size;
    return ret;
  }


  Ret json_integer_list(const char * filename,
                        std::vector<unsigned> & list)
  {
    JSONParser parser;
    return parser.parse_integer_list(filename, list);
  }


  namespace  {

    void json_write_poly(std::ofstream & out,
                         const MPReconstructedPoly & poly,
                         unsigned nvars)
    {
      out << "[";

      unsigned size = poly.size();
      for (unsigned j=0; j<size; ++j) {
        if (j)
          out << ",";

        out << "[";

        out << "\"";
        out << poly.coeff(j);
        out << "\"";

        out << ",";

        const Monomial::VarExponent * exps = poly.monomial(j).exponents();
        out << "[";
        for (unsigned k=0; k<nvars; ++k) {
          if (k)
            out << ",";
          out << exps[k];
        }
        out << "]";

        out << "]";
      }

      out << "]";
    }

    void json_write_ratfun_impl(std::ofstream & out,
                                const MPReconstructedRatFun & fun,
                                unsigned nvars)
    {
      if (fun.numerator().size() == 0) {
        out << "0";
        return;
      }

      out << "[";

      {
        out << "[";
        out << fun.numerator().size();
        out << ",";
        json_write_poly(out, fun.numerator(), nvars);
        out << "]";
      }

      out << ",";

      {
        out << "[";
        out << fun.denominator().size();
        out << ",";
        json_write_poly(out, fun.denominator(), nvars);
        out << "]";
      }

      out << "]";
    }

  } // namespace

  Ret json_write_ratfun(const char * file,
                        const MPReconstructedRatFun * fun,
                        unsigned len,
                        unsigned nvars)
  {
    std::ofstream out(file);
    if (out.fail())
      return FAILED;

    out << "[";
    out << len;
    out << ",";
    out << nvars;
    out << ",";

    out << "[";
    for (unsigned j=0; j<len; ++j) {
      if (j)
        out << ",";
      json_write_ratfun_impl(out, fun[j], nvars);
    }
    out << "]";

    out << "]";

    return SUCCESS;
  }

  Ret json_write_integer_list(const char * filename,
                              const unsigned * list,
                              unsigned len)
  {
    std::ofstream out(filename);
    if (out.fail())
      return FAILED;

    out << "[";
    {
      out << len;
      out << ",";
      out << "[";
      for (unsigned j=0; j<len; ++j) {
        if (j)
          out << ",";
        out << list[j];
      }
      out << "]";
    }
    out << "]";

    return SUCCESS;
  }

} // namespace fflow

