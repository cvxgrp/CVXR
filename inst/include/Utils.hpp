//   Copyright 2017 Steven Diamond
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

// Some useful defines for Matricies/etc.

#include "../include/Eigen/Sparse"
#include "../include/Eigen/Core"

#ifdef _R_INTERFACE_

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/lexical_cast.hpp>         // lexical cast
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

#include <Rcpp.h>
#include <RcppEigen.h>

#endif

#define NULL_MATRIX Eigen::SparseMatrix<double>(0,0)

typedef Eigen::Matrix<int, Eigen::Dynamic, 1> Vector;
typedef Eigen::SparseMatrix<double> Matrix;
typedef std::map<int, Matrix> CoeffMap;
typedef Eigen::Triplet<double> Triplet;
