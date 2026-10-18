//    This file is part of cvxr

#ifdef _R_INTERFACE_
#include "CVXR.h"
#include "CVXcanon.h"

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
  Rcpp::XPtr<LinOpVector> ptrX(xp);
  std::map<int, int> id_to_col = parse_named_intvec(v);
  std::map<int, int> param_to_size = parse_named_intvec(param_to_size_v);

  Rcpp::XPtr<ProblemData> resPtr(new ProblemData(), true);
  build_matrix_2(ptrX->linvec, var_length, id_to_col, param_to_size, resPtr);
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

  Rcpp::XPtr<ProblemData> resPtr(new ProblemData(), true);
  build_matrix_3(ptrX->linvec, var_length, id_to_col, param_to_size,
                 constr_offsets, resPtr);
  return resPtr;
}

#endif
