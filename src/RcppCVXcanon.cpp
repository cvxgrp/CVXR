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

//' Get the \code{sparse} flag field for the LinOp object
//'
//' @param xp the LinOpVector Object XPtr
//' @param v the \code{id_to_col} named int vector in R with integer names
//' @return a XPtr to ProblemData Object
// [[Rcpp::export(.build_matrix_0)]]
SEXP build_matrix_0(SEXP xp, Rcpp::IntegerVector v) {
  // grab the object as a XPtr (smart pointer)
#ifdef _R_DEBUG_
  Rcpp::Rcout << "In Build Matrix 0" <<std::endl;
#endif
  Rcpp::XPtr<LinOpVector> ptrX(xp);
  std::map<int, int> id_to_col;
#ifdef _R_DEBUG_
  Rcpp::Rcout << "Build Matrix 0: v.size() is " << v.size() <<std::endl;
#endif
  if (v.size() > 0) {
    Rcpp::StringVector s = v.names();
    for (int i = 0; i < s.size(); i++) {
#ifdef _R_DEBUG_
      Rcpp::Rcout << "Build_matrix_0 loop: i is " << i<<std::endl;
      Rcpp::Rcout << "Build_matrix_0 loop: s[i] is " << s[i]<<std::endl;
      Rcpp::Rcout << "Build_matrix_0 loop: v[i] is " << v[i]<<std::endl;
#endif
      id_to_col[atoi(s[i])] = v[i];
    }
  }
  //  ProblemData res = build_matrix(ptrX->linvec, id_to_col);
  //  Rcpp::Rcout << "After Build Matrix" <<std::endl;  
  //  Rcpp::XPtr<ProblemData> resPtr(&res, true);

  Rcpp::XPtr<ProblemData> resPtr(new ProblemData(), true);
#ifdef _R_DEBUG_
  Rcpp::Rcout << "After resPtr" <<std::endl;
#endif    
  build_matrix_2(ptrX->linvec, id_to_col, resPtr);
#ifdef _R_DEBUG_
  Rcpp::Rcout << "After constructing external ptr" <<std::endl;
#endif  
  return resPtr;
}

//' Get the \code{sparse} flag field for the LinOp object
//'
//' @param xp the LinOpVector Object XPtr
//' @param v1 the \code{id_to_col} named int vector in R with integer names
//' @param v2 the \code{constr_offsets} vector of offsets (an int vector in R)
//' @return a XPtr to ProblemData Object
// [[Rcpp::export(.build_matrix_1)]]
SEXP build_matrix_1(SEXP xp, Rcpp::IntegerVector v1, Rcpp::IntegerVector v2) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOpVector> ptrX(xp);
  std::map<int, int> id_to_col;
  if (v1.size() > 0 ) {
    Rcpp::StringVector s = v1.names();
    for (int i = 0; i < s.size(); i++) {
      id_to_col[atoi(s[i])] = v1[i];
    }
  }
  std::vector<int> constr_offsets;
  for (int i = 0; i < v2.size(); i++) {
    constr_offsets.push_back(v2[i]);
  }

#ifdef _R_DEBUG_  
  Rcpp::Rcout << "Before Build Matrix 1" <<std::endl;
#endif
  // ProblemData res = build_matrix(ptrX->linvec, id_to_col, constr_offsets);
  // Rcpp::Rcout << "After Build Matrix" <<std::endl;    
  // Rcpp::XPtr<ProblemData> resPtr(&res, true);
  Rcpp::XPtr<ProblemData> resPtr(new ProblemData(), true);
  build_matrix_3(ptrX->linvec, id_to_col, constr_offsets, resPtr);
#ifdef _R_DEBUG_  
  Rcpp::Rcout << "After constructing external ptr" <<std::endl;
#endif
  return resPtr;
}

#endif
