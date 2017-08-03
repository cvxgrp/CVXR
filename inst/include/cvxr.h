#ifndef RCPP_LINOP_VECTOR_H
#define RCPP_LINOP_VECTOR_H

#include "LinOp.hpp"

class LinOpVector {
public:
  /* The vector */
  std::vector<LinOp *> linvec;

  /* uuid */
  boost::uuids::uuid id;

  LinOpVector() {
    id = boost::uuids::random_generator()();
#ifdef _R_DEBUG
    Rcpp::Rcout << "LinOpVector id " << id << " Created!" << std::endl;
#endif
  }
  ~LinOpVector() {
#ifdef _R_DEBUG
    Rcpp::Rcout << "LinOpVector id " << id << " Destroyed!" << std::endl;
#endif

  }

};

#endif
