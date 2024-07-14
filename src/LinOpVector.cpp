#include "cvxcore.hpp"

//' Create a new LinOpVector object.
//'
//' @return an external ptr (Rcpp::XPtr) to a LinOp object instance.
// [[Rcpp::export(.LinOpVector__new)]]
SEXP LinOpVector__new() {
  // create a pointer to an LinopVector object and wrap it
  // as an external pointer
  Rcpp::XPtr<LinOpVector> ptr( new LinOpVector( ), true );

  // return the external pointer to the R side
  return ptr;
}


//' Perform a push back operation on the \code{args} field of LinOp
//'
//' @param xp the LinOpVector Object XPtr
//' @param yp the LinOp Object XPtr to push
// [[Rcpp::export(.LinOpVector__push_back)]]
void LinOpVector__push_back(SEXP xp, SEXP yp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOpVector> ptrX(xp);
  Rcpp::XPtr<LinOp> ptrY(yp);

  ptrX->linop_vec.push_back(ptrY);
}


//' Return the LinOp element at index i (0-based)
//'
//' @param lvec the LinOpVector Object XPtr
//' @param i the index
// [[Rcpp::export(.LinOp_at_index)]]
SEXP LinOp_at_index(SEXP lvec, int i) {
  Rcpp::XPtr<LinOpVector> vPtr(lvec);
  return Rcpp::XPtr<LinOp>((vPtr->linop_vec)[i]);
}


//' Create a new ConstLinOpVector object.
//'
//' @return an external ptr (Rcpp::XPtr) to a LinOp object instance.
// [[Rcpp::export(.ConstLinOpVector__new)]]
SEXP ConstLinOpVector__new() {
  // create a pointer to an LinopVector object and wrap it
  // as an external pointer
  Rcpp::XPtr<ConstLinOpVector> ptr( new ConstLinOpVector( ), true );

  // return the external pointer to the R side
  return ptr;
}


//' Perform a push back operation on the \code{args} field of LinOp
//'
//' @param xp the ConstLinOpVector Object XPtr
//' @param yp the LinOp Object XPtr to push
// [[Rcpp::export(.ConstLinOpVector__push_back)]]
void ConstLinOpVector__push_back(SEXP xp, SEXP yp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ConstLinOpVector> ptrX(xp);
  Rcpp::XPtr<LinOp> ptrY(yp);

  ptrX->linop_vec.push_back(ptrY);
}


//' Return the LinOp element at index i (0-based)
//'
//' @param lvec the ConstLinOpVector Object XPtr
//' @param i the index
// [[Rcpp::export(.ConstLinOp_at_index)]]
SEXP ConstLinOp_at_index(SEXP lvec, int i) {
  Rcpp::XPtr<ConstLinOpVector> vPtr(lvec);
  return Rcpp::XPtr<LinOp>((vPtr->linop_vec)[i]);
}
