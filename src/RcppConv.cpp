#include <Rcpp.h>
#include <RcppEigen.h>

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

// [[Rcpp::export(.sweep_dgCmat_in_place)]]
void multiply_dgCMatrix_vector(Rcpp::S4 A, Rcpp::NumericVector c_part) {
  // Get the slot values from the dgCMatrix object
  Rcpp::IntegerVector i = A.slot("i");
  Rcpp::IntegerVector p = A.slot("p");
  Rcpp::NumericVector x_values = A.slot("x");
  
  // Get the number of columns in the matrix
  int n = p.length() - 1;
  int k = c_part.length();
  
  // Check if the length of the vector matches the number of columns or for a scalar
  if (k != n && k != 1) {
    Rcpp::stop("mutiply_dgCMatrix_vector: Incompatible dimensions");
  }

  if (k == 1) { // c_part is a scalar
    for (int i = 0; i < x_values.length(); ++i) {
      x_values[i] *= c_part[0];
    }
  } else {
    // Perform in-place multiplication
    for (int col = 0; col < n; ++col) {
      int start = p[col];
      int end = p[col + 1];
      
      for (int idx = start; idx < end; ++idx) {
	x_values[idx] *= c_part[col];
      }
    }
  }
  
}

// [[Rcpp::export(.sweep_in_place)]]
void sweep_in_place(Rcpp::NumericMatrix P, Rcpp::NumericVector c_part) {
  int nrow = P.nrow();
  int ncol = P.ncol();
  int k = c_part.length();
  if (ncol != k && k != 1) {
    Rcpp::stop("sweep_in_place: Incompatible dimensions");
  }

  if (k == 1) { // x is a scalar
    for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < nrow; i++) {
	P(i, j) = P(i, j) * c_part[0];
      }
    }
  } else {
    for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < nrow; i++) {
	P(i, j) = P(i, j) * c_part[j];
      }
    }
  }
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::SparseMatrix<double> upper_tri_to_full(int n) {
  if (n == 0) {
    Eigen::SparseMatrix<double> result(0, 0);
    return result;
  } else {
    int entries = floor((double) n * ((double) n + 1.0) / 2.0);
    std::vector<Eigen::Triplet<double>> triplets;  
    int count = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = i; j < n; ++j) {
	triplets.emplace_back(j * n + i, count, 1.0);      
	if (i != j) {
	  triplets.emplace_back(i * n + j, count, 1.0);      
	}
	count++;
      }
    }
    Eigen::SparseMatrix<double> result(n * n, entries);
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
  }
}
