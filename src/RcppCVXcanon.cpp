//    This file is part of cvxr

#ifdef _R_INTERFACE_
#include "CVXR.h"
#include "CVXcanon.h"

static const char hexArray[] = { '0', '1', '2', '3', '4',
				 '5', '6', '7', '8', '9',
				 'a', 'b', 'c', 'd', 'e', 'f' };
// Generate a random 16-hex digit id. grouped by four
std::string genRandomId() {
  char result[20];
  Rcpp::NumericVector unif = Rcpp::runif(16, 0.0, 16.0);
  for (int i = 0, j = 0; i < 19; ++i) {
    if (i == 4 || i == 9 || i == 14) {
      result[i] = '-';
    } else {
      result[i] = hexArray[static_cast<int>(unif[j])];
      j++;
    }
  }
  result[19] = '\0';
  return(std::string(result));
}

// Make a map out of an Rcpp List. Caller's responsibility
// to ensure proper names etc.
std::map<std::string, double>  makeMap(Rcpp::List L) {
  std::map<std::string, double> result;
  Rcpp::StringVector s = L.names();
  for (int i = 0; i < s.size(); i++) {
    result[Rcpp::as<std::string>(s[i])] = L[i];
  }
  return(result);
}

// Helper: parse a named IntegerVector into map<int,int>
static std::map<int, int> parse_named_intvec(Rcpp::IntegerVector v) {
  std::map<int, int> result;
  if (v.size() > 0) {
    Rcpp::StringVector s = v.names();
    for (int i = 0; i < s.size(); i++) {
      result[std::stoi(Rcpp::as<std::string>(s[i]))] = v[i];
    }
  }
  return result;
}

// Build the problem matrix (no constraint offsets).
//
// @param xp the LinOpVector Object XPtr
// @param v the \code{id_to_col} named int vector in R with integer names
// @param var_length total number of variable columns
// @param param_to_size_v named int vector: param_id -> param_size (empty for non-DPP)
// @return a XPtr to ProblemData Object
// [[Rcpp::export(.build_matrix_0)]]
SEXP build_matrix_0(SEXP xp, Rcpp::IntegerVector v, int var_length,
                    Rcpp::IntegerVector param_to_size_v) {
#ifdef _R_DEBUG_
  Rcpp::Rcout << "In Build Matrix 0" <<std::endl;
#endif
  Rcpp::XPtr<LinOpVector> ptrX(xp);
  std::map<int, int> id_to_col = parse_named_intvec(v);
  std::map<int, int> param_to_size = parse_named_intvec(param_to_size_v);

  Rcpp::XPtr<ProblemData> resPtr(new ProblemData(), true);
#ifdef _R_DEBUG_
  Rcpp::Rcout << "After resPtr" <<std::endl;
#endif
  build_matrix_2(ptrX->linvec, var_length, id_to_col, param_to_size, resPtr);
#ifdef _R_DEBUG_
  Rcpp::Rcout << "After constructing external ptr" <<std::endl;
#endif
  return resPtr;
}

// Build the problem matrix (with constraint offsets).
//
// @param xp the LinOpVector Object XPtr
// @param v1 the \code{id_to_col} named int vector in R with integer names
// @param var_length total number of variable columns
// @param param_to_size_v named int vector: param_id -> param_size (empty for non-DPP)
// @param v2 the \code{constr_offsets} vector of offsets (an int vector in R)
// @return a XPtr to ProblemData Object
// [[Rcpp::export(.build_matrix_1)]]
SEXP build_matrix_1(SEXP xp, Rcpp::IntegerVector v1, int var_length,
                    Rcpp::IntegerVector param_to_size_v,
                    Rcpp::IntegerVector v2) {
  Rcpp::XPtr<LinOpVector> ptrX(xp);
  std::map<int, int> id_to_col = parse_named_intvec(v1);
  std::map<int, int> param_to_size = parse_named_intvec(param_to_size_v);

  std::vector<int> constr_offsets;
  for (int i = 0; i < v2.size(); i++) {
    constr_offsets.push_back(v2[i]);
  }

#ifdef _R_DEBUG_
  Rcpp::Rcout << "Before Build Matrix 1" <<std::endl;
#endif
  Rcpp::XPtr<ProblemData> resPtr(new ProblemData(), true);
  build_matrix_3(ptrX->linvec, var_length, id_to_col, param_to_size,
                 constr_offsets, resPtr);
#ifdef _R_DEBUG_
  Rcpp::Rcout << "After constructing external ptr" <<std::endl;
#endif
  return resPtr;
}

#endif
