#include "Solution.hpp"

//' Create a new Solution object.
//'
//' @return an external ptr (Rcpp::XPtr) to a Solution object instance.
// [[Rcpp::export]]
SEXP Solution__new() {
  // create a pointer to an Uniform object and wrap it
  // as an external pointer
  Rcpp::XPtr<Solution> ptr( new Solution( ), true );

  // return the external pointer to the R side
  return ptr;
}

//' Return the solver status of the Solution Object
//'
//' @param xp the Solution Object XPtr
//' @return the value of the status field of the Solution Object
// [[Rcpp::export]]
int Solution__get_status(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<Solution> ptr(xp);

  // Get the result
  solverStatus code = ptr->status;
  int numericCode;
  switch (code) {
  case OPTIMAL:
    numericCode = 1;
    break;
  case INFEASIBLE:
    numericCode = 2;
    break;
  case UNBOUNDED:
    numericCode = 3;
    break;
  case OPTIMAL_INACCURATE:
    numericCode = 4;
    break;
  case INFEASIBLE_INACCURATE:
    numericCode = 5;
    break;
  case UNBOUNDED_INACCURATE:
    numericCode = 6;
    break;
  default:
    numericCode = 7;
  }
  return numericCode;
}


//' Return the optimal value of the Solution Object
//'
//' @param xp the Solution Object XPtr
//' @return the value of the optimal_value field of the Solution Object
// [[Rcpp::export]]
double Solution__get_optimal_value(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<Solution> ptr(xp);

  // Get the result
  return ptr->optimal_value;
}

//' Return the primal_values map of the Solution Object
//'
//' @param xp the Solution Object XPtr
//' @return the primal_values map
// [[Rcpp::export]]
std::map<int, Eigen::MatrixXd> Solution__get_primal_values(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<Solution> ptr(xp);

  // Get the result
  return ptr->primal_values;
}

//' Return the dual_values map of the Solution Object
//'
//' @param xp the Solution Object XPtr
//' @return the dual_values map
// [[Rcpp::export]]
std::map<int, Eigen::MatrixXd> Solution__get_dual_values(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<Solution> ptr(xp);

  // Get the result
  return ptr->dual_values;
}

