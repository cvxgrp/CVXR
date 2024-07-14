#ifndef _LINOP_VECTOR_H_
#define _LINOP_VECTOR_H_

#include "Utils.h"
#include "LinOp.h"
#include "ProblemData.h"

class LinOpVector {
public:
  /* The vector */
  std::vector<LinOp *> linop_vec;
  
#ifdef _R_INTERFACE_
#ifdef _R_DEBUG
  /* almost uuid */
  std::string id;
#endif
#endif
  
  LinOpVector() {
#ifdef _R_INTERFACE_
#ifdef _R_DEBUG
    id = genRandomId();
    Rcpp::Rcout << "LinOpVector id " << id << " Created!" << std::endl;
#endif
#endif
  }

  ~LinOpVector() {
#ifdef _R_INTERFACE_    
#ifdef _R_DEBUG
    Rcpp::Rcout << "LinOpVector id " << id << " Destroyed!" << std::endl;
#endif
#endif
  }
  
};

class ConstLinOpVector {
public:
  /* The vector */
  std::vector<const LinOp *> linop_vec;
  
#ifdef _R_INTERFACE_
#ifdef _R_DEBUG
  /* almost uuid */
  std::string id;
#endif
#endif
  
  ConstLinOpVector() {
#ifdef _R_INTERFACE_
#ifdef _R_DEBUG
    id = genRandomId();
    Rcpp::Rcout << "LinOpVector id " << id << " Created!" << std::endl;
#endif
#endif
  }

  ~ConstLinOpVector() {
#ifdef _R_INTERFACE_    
#ifdef _R_DEBUG
    Rcpp::Rcout << "LinOpVector id " << id << " Destroyed!" << std::endl;
#endif
#endif
  }
  
};



#endif
