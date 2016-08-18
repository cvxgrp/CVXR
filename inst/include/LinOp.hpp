//    This file is part of CVXcanon.
//
//    CVXcanon is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    CVXcanon is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with CVXcanon.  If not, see <http://www.gnu.org/licenses/>.

#ifndef LINOP_H
#define LINOP_H

#include <vector>
#include <cassert>
#include <iostream>
#include "Utils.hpp"

/* ID for all coefficient matrices associated with linOps of CONSTANT_TYPE */
static const int CONSTANT_ID = -1;

/* TYPE of each LinOP */
enum operatortype {
  VARIABLE,
  PROMOTE,
  MUL,
  RMUL,
  MUL_ELEM,
  DIV,
  SUM,
  NEG,
  INDEX,
  TRANSPOSE,
  SUM_ENTRIES,
  TRACE,
  RESHAPE,
  DIAG_VEC,
  DIAG_MAT,
  UPPER_TRI,
  CONV,
  HSTACK,
  VSTACK,
  SCALAR_CONST,
  DENSE_CONST,
  SPARSE_CONST,
  NO_OP,
  KRON,

  /* Used only to specify constraint type */
  EQ,     // equality constraint
  LEQ,    // non-negative orthant
  SOC,    // second-order cone
  EXP,    // exponential cone
  SDP,    // semi-definite cone
};

/* linOp TYPE */
typedef operatortype OperatorType;

/* LinOp Class mirrors the CVXPY linOp class. Data fields are determined
   by the TYPE of LinOp. No error checking is performed on the data fields,
   and the semantics of SIZE, ARGS, and DATA depends on the linop TYPE. */
class LinOp {
public:
  OperatorType type;
  std::vector<int> size;

  /* Children LinOps in the tree */
  std::vector<LinOp*> args;

  /* Sparse Data Fields */
  bool sparse; // True only if linOp has sparse_data
  Matrix sparse_data;

  /* Dense Data Field */
  Eigen::MatrixXd dense_data;

  /* Slice Data: stores slice data as (row_slice, col_slice)
   * where slice = (start, end, step_size) */
  std::vector<std::vector<int> > slice;

  /* uuid */
  boost::uuids::uuid id;

  /* Constructor */
  LinOp() {
    id = boost::uuids::random_generator()();
#ifdef CVXCANON_DEBUG
    Rcpp::Rcout << "LinOp id " << id << " Created!" << std::endl;
#endif
    sparse = false; // dense by default
  }

  ~LinOp() {
#ifdef CVXCANON_DEBUG
    Rcpp::Rcout << "LinOp id " << id << " ; type " << type << " Destroyed!!" << std::endl;
#endif
  }

  /* Checks if LinOp is constant type */
  bool has_constant_type() {
    return  type == SCALAR_CONST || type == DENSE_CONST
            || type == SPARSE_CONST;
  }

  /* Initializes DENSE_DATA. MATRIX is a pointer to the data of a 2D
   * numpy array, ROWS and COLS are the size of the ARRAY.
   *
   * MATRIX must be a contiguous array of doubles aligned in fortran
   * order.
   *
   * NOTE: The function prototype must match the type-map in CVXCanon.i
   * exactly to compile and run properly.
   */
  void set_dense_data(double* matrix, int rows, int cols) {
    dense_data = Eigen::Map<Eigen::MatrixXd> (matrix, rows, cols);
  }

  /* Initializes SPARSE_DATA from a sparse matrix in COO format.
   * DATA, ROW_IDXS, COL_IDXS are assumed to be contiguous 1D numpy arrays
   * where (DATA[i], ROW_IDXS[i], COLS_IDXS[i]) is a (V, I, J) triplet in
   * the matrix. ROWS and COLS should refer to the size of the matrix.
   *
   * NOTE: The function prototype must match the type-map in CVXCanon.i
   * exactly to compile and run properly.
   */
  void set_sparse_data(double *data, int data_len, double *row_idxs,
                       int rows_len, double *col_idxs, int cols_len,
                       int rows, int cols) {

    assert(rows_len == data_len && cols_len == data_len);
    sparse = true;
    Matrix sparse_coeffs(rows, cols);
    std::vector<Triplet> tripletList;
    tripletList.reserve(data_len);
    for (int idx = 0; idx < data_len; idx++) {
      tripletList.push_back(Triplet(int(row_idxs[idx]), int(col_idxs[idx]),
                                    data[idx]));
    }
    sparse_coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
    sparse_coeffs.makeCompressed();
    sparse_data = sparse_coeffs;
  }
};

struct Variable {
  int id;
  std::vector<int> size;
  OperatorType type; // dual variables only
  bool operator < (const Variable &other) const { return id < other.id; }
};
#endif
