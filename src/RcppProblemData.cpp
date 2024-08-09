#include "Utils.hpp"
#include "ProblemData.hpp"

//' Create a new ProblemData object.
//'
//' @return an external ptr (Rcpp::XPtr) to a ProblemData object instance.
// [[Rcpp::export(.ProblemData__new)]]
SEXP ProblemData__new() {
  // create a pointer to an Uniform object and wrap it
  // as an external pointer
  Rcpp::XPtr<ProblemData> ptr( new ProblemData( ), true );

  // return the external pointer to the R side
  return ptr;
}

//' Get the vec_idx field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return an integer vector of the field vec_idx from the ProblemData Object
// [[Rcpp::export(.ProblemData__get_vec_idx)]]
int ProblemData__get_vec_idx(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);
  return ptr->vec_idx;
}

//' Set the vec_idx field in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param idx an int value for field vec_idx of the ProblemData object
// [[Rcpp::export(.ProblemData__set_vec_idx)]]
void ProblemData__set_vec_idx(SEXP xp, int idx) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);
  // Set the result
  ptr->vec_idx = idx;
}

//' Get the param_id field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return an integer vector of the field param_id from the ProblemData Object
// [[Rcpp::export(.ProblemData__get_param_id)]]
int ProblemData__get_param_id(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);
  return ptr->param_id;
}

//' Set the param_id field in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param idx an int value for field param_id of the ProblemData object
// [[Rcpp::export(.ProblemData__set_param_id)]]
void ProblemData__set_param_id(SEXP xp, int idx) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);
  // Set the result
  ptr->param_id = idx;
}

//' Get the length of V, I, J
//'
//' @param xp the ProblemData Object XPtr
//' @return the length of V, I, J 
// [[Rcpp::export(.ProblemData__getLen)]]
std::vector<double> ProblemData__getLen(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // return the result to R
  return ptr->getLen();
}

//' Get the V field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param num_values the number of values to return
//' @return a numeric vector of doubles (the field V) from the ProblemData Object
// [[Rcpp::export(.ProblemData__get_V)]]
Rcpp::NumericVector ProblemData__get_V(SEXP xp, int num_values) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);
  Rcpp::NumericVector result(num_values);
  ptr->getV(result.begin(), num_values);
  // return the result to R
  return result;
}

//' Get the I field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param num_values the number of values to return
//' @return an integer vector of the field I from the ProblemData Object
// [[Rcpp::export(.ProblemData__get_I)]]
Rcpp::IntVector ProblemData__get_I(SEXP xp, int num_values) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);
  Rcpp::IntVector result(num_values);
  ptr->getI(result.begin(), num_values);
  // return the result to R
  return result;
}

//' Get the J field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param num_values the number of values to return
//' @return an integer vector of the field J from the ProblemData Object
// [[Rcpp::export(.ProblemData__get_J)]]
Rcpp::IntVector ProblemData__get_J(SEXP xp, int num_values) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);
  Rcpp::IntVector result(num_values);
  ptr->getJ(result.begin(), num_values);
  // return the result to R
  return result;
}


