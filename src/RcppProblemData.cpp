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
