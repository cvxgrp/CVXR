#ifndef _CVXR_H_
#define _CVXR_H_

#include "Utils.h"
#include "LinOp.h"
#include "ProblemData.h"

class LinOpVector {
public:
  /* The vector */
  std::vector<LinOp *> linvec;

  LinOpVector() {
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
