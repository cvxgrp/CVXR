#include "ProblemData.hpp"

//' Create a new ProblemData object.
//'
//' @return an external ptr (Rcpp::XPtr) to a ProblemData object instance.
// [[Rcpp::export]]
SEXP ProblemData__new() {
  // create a pointer to an Uniform object and wrap it
  // as an external pointer
  Rcpp::XPtr<ProblemData> ptr( new ProblemData( ), true );

  // return the external pointer to the R side
  return ptr;
}

//' Invoke the getV method on the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return a std::vector of doubles (the field V) from the ProblemData Object
// [[Rcpp::export]]
std::vector<double> ProblemData__getV(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // return the result to R
  std::vector<double> res = ptr->V;
  return res;
}

//' Set the V vector in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param vp the values for V
// [[Rcpp::export]]
void ProblemData__set_V(SEXP xp, std::vector<double> vp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->V = vp;
}

//' Return the const_vec field value from the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return a numeric vector of the field const_vec from the ProblemData Object
// [[Rcpp::export]]
std::vector<double> ProblemData__getConstVec(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // invoke the function
  std::vector<double> res = ptr->const_vec;

  // return the result to R
  return res;
}

//' Set the const_vec in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param cvp the values for const_vec
// [[Rcpp::export]]
void ProblemData__set_const_vec(SEXP xp, std::vector<double> cvp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->const_vec = cvp;
}


//' Invoke the getI method on the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return an integer vector of the field I from the ProblemData Object
// [[Rcpp::export]]
std::vector<int> ProblemData__getI(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // invoke the function
  std::vector<int> res = ptr->I;

  // return the result to R
  return res;
}

//' Set the I in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param ip the values for I
// [[Rcpp::export]]
void ProblemData__set_I(SEXP xp, std::vector<int> ip) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->I = ip;
}

//' Invoke the getJ method on the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return an integer vector of the field J from the ProblemData Object
// [[Rcpp::export]]
std::vector<int> ProblemData__getJ(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // invoke the function
  std::vector<int> res = ptr->J;

  // return the result to R
  return res;
}

//' Set the J in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param jp the values for J
// [[Rcpp::export]]
void ProblemData__set_J(SEXP xp, std::vector<int> jp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->J = jp;
}

//' Return the id_to_col map of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return the id_to_col map
// [[Rcpp::export]]
std::map<int, int> ProblemData__get_id_to_col(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Get the result
  return ptr->id_to_col;
}

//' Set the id_to_col map of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param iv the named integer vector
//' @return the id_to_col map
// [[Rcpp::export]]
void ProblemData__set_id_to_col(SEXP xp, Rcpp::IntegerVector iv) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  std::map<int, int> id_to_col_map = ptr->id_to_col;

  Rcpp::CharacterVector iNames = iv.names();
  Rcpp::CharacterVector::iterator ct;
  Rcpp::IntegerVector::iterator it;

  id_to_col_map.clear();
  for (ct = iNames.begin(), it = iv.begin(); ct != iNames.end(); ++it, ++ct) {
    id_to_col_map[ atoi(*ct) ] = *it;
  }
}

//' Return the const_to_row map of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return the const_to_row map
// [[Rcpp::export]]
std::map<int, int> ProblemData__get_const_to_row(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Get the result
  return ptr->const_to_row;
}

//' Set the const_to_row map of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param iv the named integer vector
//' @return the const_to_row map
// [[Rcpp::export]]
void ProblemData__set_const_to_row(SEXP xp, Rcpp::IntegerVector iv) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  std::map<int, int> const_to_row_map = ptr->const_to_row;

  Rcpp::CharacterVector iNames = iv.names();
  Rcpp::CharacterVector::iterator ct;
  Rcpp::IntegerVector::iterator it;

  const_to_row_map.clear();
  for (ct = iNames.begin(), it = iv.begin(); ct != iNames.end(); ++it, ++ct) {
    const_to_row_map[ atoi(*ct) ] = *it;
  }
}

//' Return the number of constraints of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return the number of constraints
// [[Rcpp::export]]
int ProblemData__get_num_constaints(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Get the result
  return ptr->num_constraints;
}

//' Get the vals value from the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return a std::vector of doubles (the field vals) from the ProblemData Object
// [[Rcpp::export]]
std::vector<double> ProblemData__get_vals(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // return the result to R
  std::vector<double> res = ptr->vals;
  return res;
}

//' Set the vals vector in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param vp the values for vals
// [[Rcpp::export]]
void ProblemData__set_vals(SEXP xp, std::vector<double> vp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->vals = vp;
}

//' Get the row_idxs vector from the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return an integer vector of the field row_idxs from the ProblemData Object
// [[Rcpp::export]]
std::vector<int> ProblemData__get_row_idxs(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // invoke the function
  std::vector<int> res = ptr->row_idxs;

  // return the result to R
  return res;
}

//' Set the row_idxs field in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param idxp the values for idxs
// [[Rcpp::export]]
void ProblemData__set_row_idxs(SEXP xp, std::vector<int> idxp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->row_idxs = idxp;
}

//' Get the col_ptrs field from the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return an integer vector of the field col_ptrs from the ProblemData Object
// [[Rcpp::export]]
std::vector<int> ProblemData__get_col_ptrs(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // invoke the function
  std::vector<int> res = ptr->col_ptrs;

  // return the result to R
  return res;
}

//' Set the col_ptrs field in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param icp the values for col_ptrs
// [[Rcpp::export]]
void ProblemData__set_col_ptrs(SEXP xp, std::vector<int> icp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->col_ptrs = icp;
}
