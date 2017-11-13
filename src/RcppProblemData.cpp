#include "Utils.h"
#include "ProblemData.h"

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

//' Get the V field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return a numeric vector of doubles (the field V) from the ProblemData Object
// [[Rcpp::export(.ProblemData__get_V)]]
std::vector<double> ProblemData__get_V(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // return the result to R
  return ptr->V;
}

//' Set the V field in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param vp a numeric vector of values for field V
// [[Rcpp::export(.ProblemData__set_V)]]
void ProblemData__set_V(SEXP xp, std::vector<double> vp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->V = vp;
}

//' Get the I field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return an integer vector of the field I from the ProblemData Object
// [[Rcpp::export(.ProblemData__get_I)]]
std::vector<int> ProblemData__get_I(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  return ptr->I;
}

//' Set the I field in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param ip an integer vector of values for field I of the ProblemData object
// [[Rcpp::export(.ProblemData__set_I)]]
void ProblemData__set_I(SEXP xp, std::vector<int> ip) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->I = ip;
}

//' Get the J field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return an integer vector of the field J from the ProblemData Object
// [[Rcpp::export(.ProblemData__get_J)]]
std::vector<int> ProblemData__get_J(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // return the result to R
  return ptr ->J;
}

//' Set the J field in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param jp an integer vector of the values for field J of the ProblemData object
// [[Rcpp::export(.ProblemData__set_J)]]
void ProblemData__set_J(SEXP xp, std::vector<int> jp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->J = jp;
}

//' Get the const_vec field from the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return a numeric vector of the field const_vec from the ProblemData Object
// [[Rcpp::export(.ProblemData__get_const_vec)]]
std::vector<double> ProblemData__get_const_vec(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  return ptr->const_vec;
}

//' Set the const_vec field in the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param cvp a numeric vector of values for const_vec field of the ProblemData object
// [[Rcpp::export(.ProblemData__set_const_vec)]]
void ProblemData__set_const_vec(SEXP xp, std::vector<double> cvp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Set the result
  ptr->const_vec = cvp;
}


//' Get the id_to_col field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return the id_to_col field as a named integer vector where the names are integers converted to characters
// [[Rcpp::export(.ProblemData__get_id_to_col)]]
std::map<int, int> ProblemData__get_id_to_col(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Get the result
  return ptr->id_to_col;
}

//' Set the id_to_col field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param iv a named integer vector with names being integers converted to characters
// [[Rcpp::export(.ProblemData__set_id_to_col)]]
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

//' Get the const_to_row field of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @return the const_to_row field as a named integer vector where the names are integers converted to characters
// [[Rcpp::export(.ProblemData__get_const_to_row)]]
std::map<int, int> ProblemData__get_const_to_row(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ProblemData> ptr(xp);

  // Get the result
  return ptr->const_to_row;
}

//' Set the const_to_row map of the ProblemData Object
//'
//' @param xp the ProblemData Object XPtr
//' @param iv a named integer vector with names being integers converted to characters
// [[Rcpp::export(.ProblemData__set_const_to_row)]]
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

