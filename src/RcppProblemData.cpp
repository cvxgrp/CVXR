#include <climits>
#include "Utils.h"
#include "ProblemData.h"

// Create a new ProblemData object.
//
// @return an external ptr (Rcpp::XPtr) to a ProblemData object instance.
// [[Rcpp::export(.ProblemData__new)]]
SEXP ProblemData__new() {
  Rcpp::XPtr<ProblemData> ptr( new ProblemData( ), true );
  return ptr;
}

// Set the param_id extraction pointer.
//
// @param xp the ProblemData Object XPtr
// @param pid the param_id to set
// [[Rcpp::export(.ProblemData__set_param_id)]]
void ProblemData__set_param_id(SEXP xp, int pid) {
  Rcpp::XPtr<ProblemData> ptr(xp);
  ptr->param_id = pid;
}

// Set the vec_idx extraction pointer.
//
// @param xp the ProblemData Object XPtr
// @param idx the vec_idx to set
// [[Rcpp::export(.ProblemData__set_vec_idx)]]
void ProblemData__set_vec_idx(SEXP xp, int idx) {
  Rcpp::XPtr<ProblemData> ptr(xp);
  ptr->vec_idx = idx;
}

// Get the length of V for current (param_id, vec_idx).
//
// @param xp the ProblemData Object XPtr
// @return integer length
// [[Rcpp::export(.ProblemData__getLen)]]
int ProblemData__getLen(SEXP xp) {
  Rcpp::XPtr<ProblemData> ptr(xp);
  return ptr->getLen();
}

// Get the V data for current (param_id, vec_idx).
//
// @param xp the ProblemData Object XPtr
// @return numeric vector
// [[Rcpp::export(.ProblemData__getV)]]
std::vector<double> ProblemData__getV(SEXP xp) {
  Rcpp::XPtr<ProblemData> ptr(xp);
  return ptr->getV();
}

// Get the I data for current (param_id, vec_idx).
//
// @param xp the ProblemData Object XPtr
// @return integer vector
// [[Rcpp::export(.ProblemData__getI)]]
std::vector<int> ProblemData__getI(SEXP xp) {
  Rcpp::XPtr<ProblemData> ptr(xp);
  return ptr->getI();
}

// Get the J data for current (param_id, vec_idx).
//
// @param xp the ProblemData Object XPtr
// @return integer vector
// [[Rcpp::export(.ProblemData__getJ)]]
std::vector<int> ProblemData__getJ(SEXP xp) {
  Rcpp::XPtr<ProblemData> ptr(xp);
  return ptr->getJ();
}

// Get all param_ids present in the tensor.
//
// @param xp the ProblemData Object XPtr
// @return integer vector of param_ids
// [[Rcpp::export(.ProblemData__get_param_ids)]]
std::vector<int> ProblemData__get_param_ids(SEXP xp) {
  Rcpp::XPtr<ProblemData> ptr(xp);
  return ptr->get_param_ids();
}

// Get the number of inner vectors for a given param_id.
//
// @param xp the ProblemData Object XPtr
// @param pid the param_id
// @return integer count
// [[Rcpp::export(.ProblemData__get_num_vecs)]]
int ProblemData__get_num_vecs(SEXP xp, int pid) {
  Rcpp::XPtr<ProblemData> ptr(xp);
  return ptr->get_num_vecs(pid);
}

// Flatten every (param_id, element) slice of the tensor in one call.
//
// This is the extraction-and-concatenation block of CVXPY's
// cppbackend.py build_matrix(), lines 66-91, moved into C++: CVXPY walks the
// slices with three getters per parameter element from Python, and the R
// mirror of that loop cost four .Call()s per element (about 1 ms per call at
// 500 parameter elements, independent of nnz; claims audit decision 30,
// CVXR 1.9.2.9165). The three input vectors are parallel, in the caller's
// iteration order; the result is the tensor's flat COO triple with 0-based
// indices, in the same entry order the R loop produced:
//   I = row + col * constr_length          (cppbackend.py:88)
//   J = element index + column offset      (cppbackend.py:89)
// Every id in `param_ids` must have a slot of at least `param_sizes` inner
// vectors (init_data_tensor() in CVXcanon.cpp creates exactly that); anything
// else is a programming error and stops.
//
// @param xp the ProblemData Object XPtr
// @param param_ids integer vector of param ids (CONSTANT_ID = -1 included)
// @param param_sizes integer vector, number of elements per param id
// @param col_offsets integer vector, tensor column offset per param id
// @param constr_length total number of constraint rows across the LinOps
// @return list(V = numeric, I = integer, J = integer), 0-based I and J
// [[Rcpp::export(.ProblemData__get_tensor)]]
Rcpp::List ProblemData__get_tensor(SEXP xp, Rcpp::IntegerVector param_ids,
                                   Rcpp::IntegerVector param_sizes,
                                   Rcpp::IntegerVector col_offsets,
                                   int constr_length) {
  Rcpp::XPtr<ProblemData> ptr(xp);
  const R_xlen_t n_ids = param_ids.size();
  if (param_sizes.size() != n_ids || col_offsets.size() != n_ids) {
    Rcpp::stop("get_tensor: param_ids, param_sizes and col_offsets must have the same length");
  }

  // First pass: total number of entries, so each output is allocated once.
  R_xlen_t total = 0;
  for (R_xlen_t k = 0; k < n_ids; ++k) {
    std::map<int, std::vector<std::vector<double> > >::const_iterator it =
      ptr->TensorV.find(param_ids[k]);
    if (it == ptr->TensorV.end() || (int) it->second.size() < param_sizes[k]) {
      Rcpp::stop("get_tensor: param id %d has no tensor slot of size %d",
                 param_ids[k], param_sizes[k]);
    }
    for (int i = 0; i < param_sizes[k]; ++i) {
      total += (R_xlen_t) it->second[i].size();
    }
  }

  Rcpp::NumericVector V(total);
  Rcpp::IntegerVector I(total);
  Rcpp::IntegerVector J(total);
  const long long row_limit = (long long) INT_MAX;
  R_xlen_t pos = 0;
  for (R_xlen_t k = 0; k < n_ids; ++k) {
    const int pid = param_ids[k];
    const std::vector<std::vector<double> >& tV = ptr->TensorV[pid];
    const std::vector<std::vector<int> >& tI = ptr->TensorI[pid];
    const std::vector<std::vector<int> >& tJ = ptr->TensorJ[pid];
    for (int i = 0; i < param_sizes[k]; ++i) {
      const std::vector<double>& v = tV[i];
      const std::vector<int>& ii = tI[i];
      const std::vector<int>& jj = tJ[i];
      const int col = i + col_offsets[k];
      for (size_t e = 0; e < v.size(); ++e) {
        const long long row = (long long) ii[e] + (long long) jj[e] * (long long) constr_length;
        if (row > row_limit) {
          Rcpp::stop("get_tensor: flattened tensor row exceeds the integer limit");
        }
        V[pos] = v[e];
        I[pos] = (int) row;
        J[pos] = col;
        ++pos;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("V") = V,
                            Rcpp::Named("I") = I,
                            Rcpp::Named("J") = J);
}
