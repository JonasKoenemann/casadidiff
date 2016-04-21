/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */



#include "code_generator.hpp"
#include "function_internal.hpp"
#include <iomanip>
#include "casadi/core/runtime/runtime_embedded.hpp"

using namespace std;
namespace casadi {

  CodeGenerator::CodeGenerator(const Dict& opts) {
    // Default options
    this->verbose = false;
    this->mex = false;
    this->cpp = false;
    this->main = false;
    this->real_t = "double";
    this->codegen_scalars = false;
    this->with_header = false;

    // Read options
    for (auto&& e : opts) {
      if (e.first=="verbose") {
        this->verbose = e.second;
      } else if (e.first=="mex") {
        this->mex = e.second;
      } else if (e.first=="cpp") {
        this->cpp = e.second;
      } else if (e.first=="main") {
        this->main = e.second;
      } else if (e.first=="real_t") {
        this->real_t = e.second.to_string();
      } else if (e.first=="codegen_scalars") {
        this->codegen_scalars = e.second;
      } else if (e.first=="with_header") {
        this->with_header = e.second;
      } else {
        casadi_error("Unrecongnized option: " << e.first);
      }
    }

    // Includes needed
    if (this->main) addInclude("stdio.h");

    // Mex and main needs string.h
    if (this->mex || this->main) {
      addInclude("string.h");
    }

    // Mex file?
    if (this->mex) {
      addInclude("mex.h", false, "MATLAB_MEX_FILE");
      // Define printf (note file should be compilable both with and without mex)
      this->auxiliaries
        << "#ifdef MATLAB_MEX_FILE" << endl
        << "#define PRINTF mexPrintf" << endl
        << "#else" << endl
        << "#define PRINTF printf" << endl
        << "#endif" << endl;
    } else {
      // Define printf as standard printf from stdio.h
      this->auxiliaries << "#define PRINTF printf" << endl;
    }
  }

  void CodeGenerator::add(const Function& f) {
    f->generateFunction(*this, f.name(), false);
    if (this->with_header) {
      if (this->cpp) this->header << "extern \"C\" " ; // C linkage
      this->header << f->signature(f.name()) << ";" << endl;
    }
    f->generateMeta(*this, f.name());
    this->exposed_fname.push_back(f.name());
  }

  std::string CodeGenerator::generate() const {
    stringstream s;
    generate(s);
    return s.str();
  }

  void CodeGenerator::file_open(std::ofstream& f, const std::string& name) const {
    // Open a file for writing
    f.open(name);

    // Print header
    f << "/* This function was automatically generated by CasADi */" << endl;

    // C linkage
    if (!this->cpp) {
      f << "#ifdef __cplusplus" << endl
        << "extern \"C\" {" << endl
        << "#endif" << endl << endl;
    }
  }

  void CodeGenerator::file_close(std::ofstream& f) const {
    // C linkage
    if (!this->cpp) {
      f << "#ifdef __cplusplus" << endl
        << "} /* extern \"C\" */" << endl
        << "#endif" << endl;
    }

    // Close file(s)
    f.close();
  }

  void CodeGenerator::define_real_t(std::ostream &s) const {
    s << "#ifndef real_t" << endl
      << "#define real_t " << this->real_t << endl
      << "#define to_double(x) "
      << (this->cpp ? "static_cast<double>(x)" : "(double) x") << endl
      << "#define to_int(x) "
      << (this->cpp ? "static_cast<int>(x)" : "(int) x") << endl
      << "#endif /* real_t */" << endl << endl;
  }

  void CodeGenerator::generate(const std::string& name) const {
    // Divide name into base and suffix (.c by default)
    string basename, suffix;
    string::size_type dotpos = name.rfind('.');
    if (dotpos==string::npos) {
      basename = name;
      suffix = this->cpp ? ".cpp" : ".c";
    } else {
      basename = name.substr(0, dotpos);
      suffix = name.substr(dotpos);
    }

    // File(s) being generated, header is optional
    ofstream s;

    // Make sure that the base name is sane
    casadi_assert(Function::check_name(basename));

    // Create files
    file_open(s, basename + suffix);

    // Prefix internal symbols to avoid symbol collisions
    s << "#ifdef CODEGEN_PREFIX" << endl
         << "  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)" << endl
         << "  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID" << endl
         << "  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)" << endl
         << "#else /* CODEGEN_PREFIX */" << endl
         << "  #define CASADI_PREFIX(ID) " << basename << "_ ## ID" << endl
         << "#endif /* CODEGEN_PREFIX */" << endl << endl;

    s << this->includes.str();
    s << endl;

    // Real type (usually double)
    define_real_t(s);

    // External function declarations
    if (!added_externals_.empty()) {
      s << "/* External functions */" << endl;
      for (auto&& i : added_externals_) {
        s << i << endl;
      }
      s << endl << endl;
    }

    // Pre-C99
    s << "/* Pre-c99 compatibility */" << endl
         << "#if __STDC_VERSION__ < 199901L" << endl
         << "real_t CASADI_PREFIX(fmin)(real_t x, real_t y) { return x<y ? x : y;}" << endl
         << "#define fmin(x,y) CASADI_PREFIX(fmin)(x,y)" << endl
         << "real_t CASADI_PREFIX(fmax)(real_t x, real_t y) { return x>y ? x : y;}" << endl
         << "#define fmax(x,y) CASADI_PREFIX(fmax)(x,y)" << endl
         << "#endif" << endl << endl;

    // Generate the actual function
    generate(s);

    // Mex gateway
    if (this->mex) {
      // Begin conditional compilation
      s << "#ifdef MATLAB_MEX_FILE" << endl;

      // Function prototype
      if (this->cpp) s << "extern \"C\"" << endl; // C linkage
      s << "void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {"
           << endl;

      // Create a buffer
      size_t buf_len = 0;
      for (int i=0; i<exposed_fname.size(); ++i) {
        buf_len = std::max(buf_len, exposed_fname[i].size());
      }
      s << "  char buf[" << (buf_len+1) << "];" << endl;

      // Read string argument
      s << "  int buf_ok = --argc >= 0 && !mxGetString(*argv++, buf, sizeof(buf));" << endl;

      // Create switch
      s << "  if (!buf_ok) {" << endl
           << "    /* name error */" << endl;
      for (int i=0; i<exposed_fname.size(); ++i) {
        s << "  } else if (strcmp(buf, \"" << exposed_fname[i] << "\")==0) {" << endl
             << "    return mex_" << exposed_fname[i] << "(resc, resv, argc, argv);" << endl;
      }
      s << "  }" << endl;

      // Error
      s << "  mexErrMsgTxt(\"First input should be a command string. Possible values:";
      for (int i=0; i<exposed_fname.size(); ++i) {
        s << " '" << exposed_fname[i] << "'";
      }
      s << "\");" << endl;

      // End conditional compilation and function
      s << "}" << endl
           << "#endif" << endl;
    }

    // Generate main
    if (this->main) {
      s << "int main(int argc, char* argv[]) {" << endl;

      // Create switch
      s << "  if (argc<2) {" << endl
           << "    /* name error */" << endl;
      for (int i=0; i<exposed_fname.size(); ++i) {
        s << "  } else if (strcmp(argv[1], \"" << exposed_fname[i] << "\")==0) {" << endl
             << "    return main_" << exposed_fname[i] << "(argc-2, argv+2);" << endl;
      }
      s << "  }" << endl;

      // Error
      s << "  fprintf(stderr, \"First input should be a command string. Possible values:";
      for (int i=0; i<exposed_fname.size(); ++i) {
        s << " '" << exposed_fname[i] << "'";
      }
      s << "\\n\");" << endl;

      // End main
      s << "  return 1;" << endl
           << "}" << endl;
    }

    // Finalize file
    file_close(s);

    // Generate header
    if (this->with_header) {
      // Create a header file
      file_open(s, basename + ".h");

      // Define the real_t type (typically double)
      define_real_t(s);

      // Add declarations
      s << this->header.str();

      // Finalize file
      file_close(s);
    }
  }

  void CodeGenerator::generate(std::ostream& s) const {
    // Codegen auxiliary functions
    s << this->auxiliaries.str();

    // Print integer constants
    stringstream name;
    for (int i=0; i<integer_constants_.size(); ++i) {
      name.str(string());
      name << "CASADI_PREFIX(s" << i << ")";
      print_vector(s, name.str(), integer_constants_[i]);
      s << "#define s" << i << " CASADI_PREFIX(s" << i << ")" << endl;
    }

    // Print double constants
    for (int i=0; i<double_constants_.size(); ++i) {
      name.str(string());
      name << "CASADI_PREFIX(c" << i << ")";
      print_vector(s, name.str(), double_constants_[i]);
      s << "#define c" << i << " CASADI_PREFIX(c" << i << ")" << endl;
    }

    // Codegen body
    s << this->body.str();

    // End with new line
    s << endl;
  }

  std::string CodeGenerator::to_string(int n) {
    stringstream ss;
    ss << n;
    return ss.str();
  }

  std::string CodeGenerator::work(int n, int sz) const {
    if (n<0 || sz==0) {
      return "0";
    } else if (sz==1 && !this->codegen_scalars) {
      return "(&w" + to_string(n) + ")";
    } else {
      return "w" + to_string(n);
    }
  }

  std::string CodeGenerator::workel(int n) const {
    if (n<0) return "0";
    stringstream s;
    if (this->codegen_scalars) s << "*";
    s << "w" << n;
    return s.str();
  }

  void CodeGenerator::assign(std::ostream &s, const std::string& lhs, const std::string& rhs) {
    s << "  " << lhs << " = " << rhs << ";" << endl;
  }

  void CodeGenerator::print_vector(std::ostream &s, const std::string& name, const vector<int>& v) {
    s << "static const int " << name << "[] = {";
    for (int i=0; i<v.size(); ++i) {
      if (i!=0) s << ", ";
      s << v[i];
    }
    s << "};" << endl;
  }

  void CodeGenerator::print_vector(std::ostream &s, const std::string& name,
                                  const vector<double>& v) {
    s << "static const real_t " << name << "[] = {";
    for (int i=0; i<v.size(); ++i) {
      if (i!=0) s << ", ";
      s << constant(v[i]);
    }
    s << "};" << endl;
  }

  void CodeGenerator::addInclude(const std::string& new_include, bool relative_path,
                                 const std::string& use_ifdef) {
    // Register the new element
    bool added = added_includes_.insert(new_include).second;

    // Quick return if it already exists
    if (!added) return;

    // Ifdef opening
    if (!use_ifdef.empty()) this->includes << "#ifdef " << use_ifdef << endl;

    // Print to the header section
    if (relative_path) {
      this->includes << "#include \"" << new_include << "\"" << endl;
    } else {
      this->includes << "#include <" << new_include << ">" << endl;
    }

    // Ifdef closing
    if (!use_ifdef.empty()) this->includes << "#endif" << endl;
  }

  bool CodeGenerator::simplifiedCall(const Function& f) {
    return f->simplifiedCall();
  }

  std::string CodeGenerator::
  operator()(const Function& f, const std::string& arg,
             const std::string& res, const std::string& iw,
             const std::string& w, const std::string& mem) const {
    return f->codegen_name(*this) + "(" + arg + ", " + res + ", "
      + iw + ", " + w + ", " + mem + ")";
  }

  std::string CodeGenerator::operator()(const Function& f,
                                        const std::string& arg, const std::string& res) const {
    return f->codegen_name(*this) + "(" + arg + ", " + res + ")";
  }

  void CodeGenerator::addExternal(const std::string& new_external) {
    added_externals_.insert(new_external);
  }

  int CodeGenerator::addSparsity(const Sparsity& sp) {
    // Get the current number of patterns before looking for it
    size_t num_patterns_before = added_sparsities_.size();

    // Get index of the pattern
    const void* h = static_cast<const void*>(sp.get());
    int& ind = added_sparsities_[h];

    // Generate it if it does not exist
    if (added_sparsities_.size() > num_patterns_before) {

      // Compact version of the sparsity pattern
      std::vector<int> sp_compact = sp.compress();

      // Codegen vector
      ind = getConstant(sp_compact, true);
    }

    return ind;
  }

  std::string CodeGenerator::sparsity(const Sparsity& sp) {
    return "s" + to_string(addSparsity(sp));
  }

  int CodeGenerator::get_sparsity(const Sparsity& sp) const {
    const void* h = static_cast<const void*>(sp.get());
    PointerMap::const_iterator it=added_sparsities_.find(h);
    casadi_assert(it!=added_sparsities_.end());
    return it->second;
  }

  size_t CodeGenerator::hash(const std::vector<double>& v) {
    // Calculate a hash value for the vector
    std::size_t seed=0;
    if (!v.empty()) {
      casadi_assert(sizeof(double) % sizeof(size_t)==0);
      const int int_len = v.size()*(sizeof(double)/sizeof(size_t));
      const size_t* int_v = reinterpret_cast<const size_t*>(&v.front());
      for (size_t i=0; i<int_len; ++i) {
        hash_combine(seed, int_v[i]);
      }
    }
    return seed;
  }

  size_t CodeGenerator::hash(const std::vector<int>& v) {
    size_t seed=0;
    hash_combine(seed, v);
    return seed;
  }

  int CodeGenerator::getConstant(const std::vector<double>& v, bool allow_adding) {
    // Hash the vector
    size_t h = hash(v);

    // Try to locate it in already added constants
    pair<multimap<size_t, size_t>::iterator, multimap<size_t, size_t>::iterator> eq =
      added_double_constants_.equal_range(h);
    for (multimap<size_t, size_t>::iterator i=eq.first; i!=eq.second; ++i) {
      if (equal(v, double_constants_[i->second])) return i->second;
    }

    if (allow_adding) {
      // Add to constants
      int ind = double_constants_.size();
      double_constants_.push_back(v);
      added_double_constants_.insert(pair<size_t, size_t>(h, ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  int CodeGenerator::getConstant(const std::vector<int>& v, bool allow_adding) {
    // Hash the vector
    size_t h = hash(v);

    // Try to locate it in already added constants
    pair<multimap<size_t, size_t>::iterator, multimap<size_t, size_t>::iterator> eq =
      added_integer_constants_.equal_range(h);
    for (multimap<size_t, size_t>::iterator i=eq.first; i!=eq.second; ++i) {
      if (equal(v, integer_constants_[i->second])) return i->second;
    }

    if (allow_adding) {
      // Add to constants
      int ind = integer_constants_.size();
      integer_constants_.push_back(v);
      added_integer_constants_.insert(pair<size_t, size_t>(h, ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  void CodeGenerator::addAuxiliary(Auxiliary f) {
    // Register the new auxiliary
    bool added = added_auxiliaries_.insert(f).second;

    // Quick return if it already exists
    if (!added) return;

    // Add the appropriate function
    switch (f) {
    case AUX_COPY:
      this->auxiliaries
        << codegen_str_copy
        << codegen_str_copy_define << endl
        << endl;
      break;
    case AUX_SWAP:
      this->auxiliaries << codegen_str_swap
        << codegen_str_swap_define
        << endl;
      break;
    case AUX_SCAL:
      this->auxiliaries << codegen_str_scal
        << codegen_str_scal_define
        << endl;
      break;
    case AUX_AXPY:
      this->auxiliaries << codegen_str_axpy
        << codegen_str_axpy_define
        << endl;
      break;
    case AUX_DOT:
      this->auxiliaries
        << codegen_str_dot
        << codegen_str_dot_define << endl
        << endl;
      break;
    case AUX_BILIN:
      this->auxiliaries
        << codegen_str_bilin
        << codegen_str_bilin_define << endl
        << endl;
      break;
    case AUX_RANK1:
      this->auxiliaries
        << codegen_str_rank1
        << codegen_str_rank1_define << endl
        << endl;
      break;
    case AUX_ASUM:
      this->auxiliaries << codegen_str_asum
        << codegen_str_asum_define
        << endl;
      break;
    case AUX_IAMAX:
      this->auxiliaries << codegen_str_iamax
        << codegen_str_iamax_define
        << endl;
      break;
    case AUX_NRM2:
      this->auxiliaries << codegen_str_nrm2
        << codegen_str_nrm2_define
        << endl;
      break;
    case AUX_FILL:
      this->auxiliaries
        << codegen_str_fill
        << codegen_str_fill_define << endl
        << endl;
      break;
    case AUX_MTIMES:
      this->auxiliaries << codegen_str_mtimes
        << codegen_str_mtimes_define
        << endl;
      break;
    case AUX_SQ:
      auxSq();
      break;
    case AUX_SIGN:
      auxSign();
      break;
    case AUX_PROJECT:
      this->auxiliaries
        << codegen_str_project
        << codegen_str_project_define
        << endl << endl;
      break;
    case AUX_TRANS:
      this->auxiliaries << codegen_str_trans
        << "#define trans(x, sp_x, y, sp_y, tmp) CASADI_PREFIX(trans)(x, sp_x, y, sp_y, tmp)"
        << endl << endl;
      break;
    case AUX_TO_MEX:
      this->auxiliaries
        << "#ifdef MATLAB_MEX_FILE" << endl
        << "mxArray* CASADI_PREFIX(to_mex)(const int* sp, const real_t* x) {" << endl
        << "  int nrow = *sp++, ncol = *sp++, nnz = sp[ncol];" << endl
        << "  mxArray* p = mxCreateSparse(nrow, ncol, nnz, mxREAL);" << endl
        << "  int i;" << endl
        << "  mwIndex* j;" << endl
        << "  for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *sp++;" << endl
        << "  for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *sp++;" << endl
        << "  if (x) {" << endl
        << "    double* d = (double*)mxGetData(p);" << endl
        << "    for (i=0; i<nnz; ++i) *d++ = to_double(*x++);" << endl
        << "  }" << endl
        << "  return p;" << endl
        << "}" << endl
        << "#define to_mex(sp, x) CASADI_PREFIX(to_mex)(sp, x)" << endl
        << "#endif" << endl << endl;
      break;
    case AUX_FROM_MEX:
      addAuxiliary(AUX_FILL);
      this->auxiliaries
        << "#ifdef MATLAB_MEX_FILE" << endl
        << "real_t* CASADI_PREFIX(from_mex)(const mxArray *p, "
        << "real_t* y, const int* sp, real_t* w) {" << endl
        << "  if (!mxIsDouble(p) || mxGetNumberOfDimensions(p)!=2)" << endl
        << "    mexErrMsgIdAndTxt(\"Casadi:RuntimeError\",\"\\\"from_mex\\\" failed: "
        << "Not a two-dimensional matrix of double precision.\");" << endl
        << "  int nrow = *sp++, ncol = *sp++, nnz = sp[ncol];" << endl
        << "  const int *colind=sp, *row=sp+ncol+1;" << endl
        << "  size_t p_nrow = mxGetM(p), p_ncol = mxGetN(p);" << endl
        << "  const double* p_data = (const double*)mxGetData(p);" << endl
        << "  bool is_sparse = mxIsSparse(p);" << endl
        << "  mwIndex *Jc = is_sparse ? mxGetJc(p) : 0;" << endl
        << "  mwIndex *Ir = is_sparse ? mxGetIr(p) : 0;" << endl
        << "  if (p_nrow==1 && p_ncol==1) {" << endl
        << "    double v = is_sparse && Jc[1]==0 ? 0 : *p_data;" << endl
        << "    fill(y, nnz, v);" << endl
        << "  } else {" << endl
        << "    bool tr = false;" << endl
        << "    if (nrow!=p_nrow || ncol!=p_ncol) {" << endl
        << "      tr = nrow==p_ncol && ncol==p_nrow && (nrow==1 || ncol==1);" << endl
        << "      if (!tr) mexErrMsgIdAndTxt(\"Casadi:RuntimeError\",\"\\\"from_mex\\\""
        << " failed: Dimension mismatch.\");" << endl
        << "    }" << endl
        << "    int r,c,k;" << endl
        << "    if (is_sparse) {" << endl
        << "      if (tr) {" << endl
        << "        for (c=0; c<ncol; ++c)" << endl
        << "          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]+c*nrow]=0;" << endl
        << "        for (c=0; c<p_ncol; ++c)" << endl
        << "          for (k=Jc[c]; k<Jc[c+1]; ++k) w[c+Ir[k]*p_ncol] = p_data[k];" << endl
        << "        for (c=0; c<ncol; ++c)" << endl
        << "          for (k=colind[c]; k<colind[c+1]; ++k) y[k] = w[row[k]+c*nrow];" << endl
        << "      } else {" << endl
        << "        for (c=0; c<ncol; ++c) {" << endl
        << "          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]]=0;" << endl
        << "          for (k=Jc[c]; k<Jc[c+1]; ++k) w[Ir[k]]=p_data[k];" << endl
        << "          for (k=colind[c]; k<colind[c+1]; ++k) y[k]=w[row[k]];" << endl
        << "        }" << endl
        << "      }" << endl
        << "    } else {" << endl
        << "      for (c=0; c<ncol; ++c) {" << endl
        << "        for (k=colind[c]; k<colind[c+1]; ++k) {" << endl
        << "          y[k] = p_data[row[k]+c*nrow];" << endl
        << "        }" << endl
        << "      }" << endl
        << "    }" << endl
        << "  }" << endl
        << "  return y;" << endl
        << "}" << endl
        << "#define from_mex(p, y, sp, w) CASADI_PREFIX(from_mex)(p, y, sp, w)" << endl
        << "#endif" << endl << endl;
      break;
    }
  }

  std::string CodeGenerator::to_mex(const Sparsity& sp, const std::string& arg) {
    addAuxiliary(AUX_TO_MEX);
    stringstream s;
    s << "to_mex(" << sparsity(sp) << ", " << arg << ");";
    return s.str();
  }

  std::string CodeGenerator::from_mex(std::string& arg,
                                      const std::string& res, std::size_t res_off,
                                      const Sparsity& sp_res, const std::string& w) {
    // Handle offset with recursion
    if (res_off!=0) return from_mex(arg, res+"+"+to_string(res_off), 0, sp_res, w);

    addAuxiliary(AUX_FROM_MEX);
    stringstream s;
    s << "from_mex(" << arg
      << ", " << res << ", " << sparsity(sp_res) << ", " << w << ");";
    return s.str();
  }

  void CodeGenerator::auxSq() {
    this->auxiliaries
      << "real_t CASADI_PREFIX(sq)(real_t x) "
      << "{ return x*x;}" << endl
      << "#define sq(x) CASADI_PREFIX(sq)(x)" << endl << endl;
  }

  void CodeGenerator::auxSign() {
    this->auxiliaries
      << "real_t CASADI_PREFIX(sign)(real_t x) "
      << "{ return x<0 ? -1 : x>0 ? 1 : x;}" << endl
      << "#define sign(x) CASADI_PREFIX(sign)(x)" << endl << endl;
  }

  std::string CodeGenerator::constant(double v) {
    stringstream s;
    if (isnan(v)) {
      s << "NAN";
    } else if (isinf(v)) {
      if (v<0) s << "-";
      s << "INFINITY";
    } else {
      int v_int(v);
      if (v_int==v) {
        // Print integer
        s << v_int << ".";
      } else {
        // Print real
        std::ios_base::fmtflags fmtfl = s.flags(); // get current format flags
        s << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1) << v;
        s.flags(fmtfl); // reset current format flags
      }
    }
    return s.str();
  }

  std::string CodeGenerator::copy(const std::string& arg,
                                  std::size_t n, const std::string& res) {
    stringstream s;
    // Perform operation
    addAuxiliary(AUX_COPY);
    s << "copy(" << arg << ", " << n << ", " << res << ");";
    return s.str();
  }

  std::string CodeGenerator::fill(const std::string& res,
                                  std::size_t n, const std::string& v) {
    stringstream s;
    // Perform operation
    addAuxiliary(AUX_FILL);
    s << "fill(" << res << ", " << n << ", " << v << ");";
    return s.str();
  }

  std::string CodeGenerator::dot(int n, const std::string& x,
                                        const std::string& y) {
    addAuxiliary(AUX_DOT);
    stringstream s;
    s << "dot(" << n << ", " << x << ", " << y << ")";
    return s.str();
  }

  std::string CodeGenerator::bilin(const std::string& A, const Sparsity& sp_A,
                                   const std::string& x, const std::string& y) {
    addAuxiliary(AUX_BILIN);
    stringstream s;
    s << "bilin(" << A << ", " << sparsity(sp_A) << ", " << x << ", " << y << ")";
    return s.str();
  }

  std::string CodeGenerator::rank1(const std::string& A, const Sparsity& sp_A,
                                   const std::string& alpha, const std::string& x,
                                   const std::string& y) {
    addAuxiliary(AUX_RANK1);
    stringstream s;
    s << "rank1(" << A << ", " << sparsity(sp_A) << ", "
      << alpha << ", " << x << ", " << y << ");";
    return s.str();
  }

  std::string CodeGenerator::declare(std::string s) {
    // Add C linkage?
    if (this->cpp) {
      s = "extern \"C\" " + s;
    }

    // To header file
    if (this->with_header) {
      this->header << s << ";" << endl;
    }

    return s;
  }

  std::string
  CodeGenerator::project(const std::string& arg, const Sparsity& sp_arg,
                         const std::string& res, const Sparsity& sp_res,
                         const std::string& w) {
    // If sparsity match, simple copy
    if (sp_arg==sp_res) return copy(arg, sp_arg.nnz(), res);

    // Create call
    addAuxiliary(CodeGenerator::AUX_PROJECT);
    stringstream s;
    s << "  project(" << arg << ", " << sparsity(sp_arg) << ", " << res << ", "
      << sparsity(sp_res) << ", " << w << ");";
    return s.str();
  }

  std::string CodeGenerator::printf(const std::string& str, const std::vector<std::string>& arg) {
    addInclude("stdio.h");
    stringstream s;
    s << "PRINTF(\"" << str << "\"";
    for (int i=0; i<arg.size(); ++i) s << ", " << arg[i];
    s << ");";
    return s.str();
  }

  std::string CodeGenerator::printf(const std::string& str, const std::string& arg1) {
    std::vector<std::string> arg;
    arg.push_back(arg1);
    return printf(str, arg);
  }

  std::string CodeGenerator::printf(const std::string& str, const std::string& arg1,
                                    const std::string& arg2) {
    std::vector<std::string> arg;
    arg.push_back(arg1);
    arg.push_back(arg2);
    return printf(str, arg);
  }

  std::string CodeGenerator::printf(const std::string& str, const std::string& arg1,
                                    const std::string& arg2, const std::string& arg3) {
    std::vector<std::string> arg;
    arg.push_back(arg1);
    arg.push_back(arg2);
    arg.push_back(arg3);
    return printf(str, arg);
  }

  std::string CodeGenerator::compile(const std::string& name,
                                     const std::string& compiler) {
    // Flag to get a DLL
#ifdef __APPLE__
    string dlflag = " -dynamiclib";
#else // __APPLE__
    string dlflag = " -shared";
#endif // __APPLE__

    // File names
    string cname = name + (this->cpp ? ".cpp" : ".c"), dlname = "./" + name + ".so";

    // Remove existing files, if any
    string rm_command = "rm -rf " + cname + " " + dlname;
    int flag = system(rm_command.c_str());
    casadi_assert_message(flag==0, "Failed to remove old source");

    // Codegen it
    generate(name);

    // Compile it
    string compile_command = compiler + " " + dlflag + " " + cname + " -o " + dlname;
    flag = system(compile_command.c_str());
    casadi_assert_message(flag==0, "Compilation failed");

    // Return name of compiled function
    return dlname;
  }

  std::string CodeGenerator::mtimes(const std::string& x, const Sparsity& sp_x,
                                    const std::string& y, const Sparsity& sp_y,
                                    const std::string& z, const Sparsity& sp_z,
                                    const std::string& w, bool tr) {
    addAuxiliary(CodeGenerator::AUX_MTIMES);
    stringstream s;
    s << "mtimes(" << x << ", " << sparsity(sp_x) << ", " << y << ", " << sparsity(sp_y) << ", "
      << z << ", " << sparsity(sp_z) << ", " << w << ", " <<  (tr ? "1" : "0") << ");";
    return s.str();
  }

} // namespace casadi
