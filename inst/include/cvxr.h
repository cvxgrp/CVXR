#ifndef _CVXR_H_
#define _CVXR_H_

#include "Utils.hpp"
#include "LinOp.hpp"
#include "ProblemData.hpp"

class LinOpVector {
public:
  /* The vector */
  std::vector<LinOp *> linvec;

  /* uuid */
  boost::uuids::uuid id;

  LinOpVector() {
    boost::uuids::basic_random_generator<boost::mt19937> gen;
    // id = boost::uuids::random_generator()();
    id = gen();
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
