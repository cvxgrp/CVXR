#include "CVXR.h"

// Source:  Rcpp Introduction

// [[Rcpp::export(.cpp_convolve)]]
Rcpp::NumericVector cpp_convolve(Rcpp::NumericVector xa, Rcpp::NumericVector xb) {
  int n_xa = xa.size(), n_xb = xb.size();
  Rcpp::NumericVector xab(n_xa + n_xb - 1);
  typedef Rcpp::NumericVector::iterator vec_iterator;
  vec_iterator ia = xa.begin(), ib = xb.begin();
  vec_iterator iab = xab.begin();
  for (int i = 0; i < n_xa; i++)
    for (int j = 0; j < n_xb; j++)
      iab[i + j] += ia[i] * ib[j];
  return xab;
}


// [[Rcpp::export(.sweep_in_place)]]
void sweep_in_place(Rcpp::NumericMatrix P, Rcpp::NumericVector c_part) {
  int nrow = P.nrow();
  int ncol = P.ncol();
  if (ncol != c_part.size()) {
    Rcpp::stop("sweep_in_place: Incompatible dimensions");
  }
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      P(i, j) = P(i, j) * c_part[j];
    }
  }
}
