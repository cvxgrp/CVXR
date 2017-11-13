#ifndef _CVXR_H_
#define _CVXR_H_

#include "Utils.h"
#include "LinOp.h"
#include "ProblemData.h"

class LinOpVector {
public:
  /* The vector */
  std::vector<LinOp *> linvec;

  /* almost uuid */
  std::string id;

  LinOpVector() {
    id = genRandomId();
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

class IntVector {
public:
  /* The vector */
  std::vector<int> int_vector;

  IntVector() {
  }

};

class Slice {
public:
  /* The slice */
  std::vector<std::vector<int> > slice;

  Slice() {
  }
  
};

#endif
