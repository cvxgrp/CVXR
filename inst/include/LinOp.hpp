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

#ifndef LINOP_H
#define LINOP_H

#include "Utils.hpp"
#include <cassert>
#include <iostream>
#include <vector>

/* TYPE of each LinOP */
enum operatortype {
  VARIABLE,
  PARAM,
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
  KRON, // for backwards compatibility; equivalent to KRON_R (which is preferred)
  KRON_R,
  KRON_L
};

/* linOp TYPE */
typedef operatortype OperatorType;

#ifdef _R_INTERFACE_
OperatorType to_optype(int typeValue) {
  // typeValue starts at 1!
  OperatorType result;
  switch (typeValue - 1) {
  case 0:
    result = VARIABLE;
    break;
  case 1:
    result = PARAM;
    break;
  case 2:
    result = PROMOTE;
    break;
  case 3:
    result = MUL;
    break;
  case 4:
    result = RMUL;
    break;
  case 5:
    result = MUL_ELEM;
    break;
  case 6:
    result = DIV;
    break;
  case 7:
    result = SUM;
    break;
  case 8:
    result = NEG;
    break;
  case 9:
    result = INDEX;
    break;
  case 10:
    result = TRANSPOSE;
    break;
  case 11:
    result = SUM_ENTRIES;
    break;
  case 12:
    result = TRACE;
    break;
  case 13:
    result = RESHAPE;
    break;
  case 14:
    result = DIAG_VEC;
    break;
  case 15:
    result = DIAG_MAT;
    break;
  case 16:
    result = UPPER_TRI;
    break;
  case 17:
    result = CONV;
    break;
  case 18:
    result = HSTACK;
    break;
  case 19:
    result = VSTACK;
    break;
  case 20:
    result = SCALAR_CONST;
    break;
  case 21:
    result = DENSE_CONST;
    break;
  case 22:
    result = SPARSE_CONST;
    break;
  case 23:
    result = NO_OP;
    break;
  case 24:
    result = KRON_R;  // KRON and KRON_R are same, latter preferred!
    break;
  case 25:
    result = KRON_R;
    break;
  case 26:
    result = KRON_L;
    break;
  default:
    Rcpp::stop("Invalid operator coding specified");
  }
  return(result);
}

int from_optype(OperatorType type) {
  int result = 0;
  switch (type) {
  case VARIABLE:
    result = 1;
    break;
  case PARAM:
    result = 2;
    break;
  case PROMOTE:
    result = 3;
    break;
  case MUL:
    result = 4;
    break;
  case RMUL:
    result = 5;
    break;
  case MUL_ELEM:
    result = 6;
    break;
  case DIV:
    result = 7;
    break;
  case SUM:
    result = 8;
    break;
  case NEG:
    result = 9;
    break;
  case INDEX:
    result = 10;
    break;
  case TRANSPOSE:
    result = 11;
    break;
  case SUM_ENTRIES:
    result = 12;
    break;
  case TRACE:
    result = 13;
    break;
  case RESHAPE:
    result = 14;
    break;
  case DIAG_VEC:
    result = 15;
    break;
  case DIAG_MAT:
    result = 16;
    break;
  case UPPER_TRI:
    result = 17;
    break;
  case CONV:
    result = 18;
    break;
  case HSTACK:
    result = 19;
    break;
  case VSTACK:
    result = 20;
    break;
  case SCALAR_CONST:
    result = 21;
    break;
  case DENSE_CONST:
    result = 22;
    break;
  case SPARSE_CONST:
    result = 23;
    break;
  case NO_OP:
    result = 24;
    break;
  case KRON:
    result = 25;
    break;
  case KRON_R:
    result = 26;
    break;
  case KRON_L:
    result = 27;
    break;
  default:
    Rcpp::stop("Invalid operator coding specified");
  }
  return(result);
}

#endif  

/* LinOp Class mirrors the CVXPY linOp class. Data fields are determined
         by the TYPE of LinOp. No error checking is performed on the data
   fields,
         and the semantics of SIZE, ARGS, and DATA depends on the linOp TYPE. */
class LinOp {
public:
  LinOp(OperatorType type,  const std::vector<int> &shape,
        const std::vector<const LinOp *> &args)
    : type_(type), shape_(shape), args_(args), sparse_(false),
      data_has_been_set_(false) {
#ifdef _R_INTERFACE_
#ifdef _R_DEBUG
    id_ = genRandomId();
    Rcpp::Rcout << "New LinOp id " << id_ << std::endl;
#endif
#endif
  }

#ifdef _R_INTERFACE_  
#ifdef _R_DEBUG 
  ~LinOp() {
    Rcpp::Rcout << "LinOp id " << id_ << "; type " << type_ << " destroyed!!" << std::endl;
  }
#endif
#endif

  OperatorType get_type() const { return type_; }
  bool is_constant() const {
    return type_ == SCALAR_CONST || type_ == DENSE_CONST ||
           type_ == SPARSE_CONST;
  }
  std::vector<int> get_shape() const { return shape_; }
  const std::vector<const LinOp *> get_args() const { return args_; }
  const std::vector<std::vector<int> > get_slice() const { return slice_; }
  void push_back_slice_vec(const std::vector<int> &slice_vec) {
    slice_.push_back(slice_vec);
  }

  bool has_numerical_data() const { return data_has_been_set_; }
  const LinOp *get_linOp_data() const { return linOp_data_; }
  void set_linOp_data(const LinOp *tree) {
    assert(!data_has_been_set_);
    linOp_data_ = tree;
    data_has_been_set_ = true;
  }
  int get_data_ndim() const { return data_ndim_; }
  void set_data_ndim(int ndim) { data_ndim_ = ndim; }
  bool is_sparse() const { return sparse_; }
  const Matrix &get_sparse_data() const { return sparse_data_; }
  const Eigen::MatrixXd &get_dense_data() const { return dense_data_; }

#ifdef _R_INTERFACE_
  /* In R we directly have the dense matrix, so no need to process things.
   */
  void set_dense_data(Eigen::Map<Eigen::MatrixXd> dense_data) {
    assert(!data_has_been_set_);
    dense_data_ = dense_data;
    sparse_ = false;
    data_has_been_set_ = true;
  }     
#else
  /* Initializes DENSE_DATA. MATRIX is a pointer to the data of a 2D
   * numpy array, ROWS and COLS are the size of the ARRAY.
   *
   * MATRIX must be a contiguous array of doubles aligned in fortran
   * order.
   *
   * NOTE: The function prototype must match the type-map in CVXCanon.i
   * exactly to compile and run properly.
   */
  void set_dense_data(double *matrix, int rows, int cols) {
    assert(!data_has_been_set_);
    dense_data_ = Eigen::Map<Eigen::MatrixXd>(matrix, rows, cols);
    sparse_ = false;
    data_has_been_set_ = true;
  }
#endif
  
#ifdef _R_INTERFACE_
  /* In R we directly have the sparse matrix, so no need to process things.
   */
  void set_sparse_data(Eigen::Map<Eigen::Matrix> sparse_data) {
    assert(!data_has_been_set_);
    sparse_data_ = sparse_data;
    sparse_ = true;
    data_ndim_ = 2;
    data_has_been_set_ = true;
  }     
#else  
  /* Initializes SPARSE_DATA from a sparse matrix in COO format.
   * DATA, ROW_IDXS, COL_IDXS are assumed to be contiguous 1D numpy arrays
   * where (DATA[i], ROW_IDXS[i], COLS_IDXS[i]) is a (V, I, J) triplet in
   * the matrix. ROWS and COLS should refer to the size of the matrix.
   *
   * NOTE: The function prototype must match the type-map in CVXCanon.i
   * exactly to compile and run properly.
   */
  void set_sparse_data(double *data, int data_len, double *row_idxs,
                       int rows_len, double *col_idxs, int cols_len, int rows,
                       int cols) {
    assert(!data_has_been_set_);
    assert(rows_len == data_len && cols_len == data_len);
    sparse_ = true;
    Matrix sparse_coeffs(rows, cols);
    std::vector<Triplet> tripletList;
    tripletList.reserve(data_len);
    for (int idx = 0; idx < data_len; idx++) {
      tripletList.push_back(
          Triplet(int(row_idxs[idx]), int(col_idxs[idx]), data[idx]));
    }
    sparse_coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
    sparse_coeffs.makeCompressed();
    sparse_data_ = sparse_coeffs;
    data_ndim_ = 2;
    data_has_been_set_ = true;
  }
#endif
  
private:
  const OperatorType type_;  
  std::vector<int> shape_;
  // children of this LinOp
  std::vector<const LinOp *> args_;
  // stores slice data as (row_slice, col_slice), where slice = (start, end,
  // step_size)
  std::vector<std::vector<int> > slice_;
  // numerical data acted upon by this linOp
  const LinOp *linOp_data_;
  int data_ndim_;
  // true iff linOp has sparse_data; at most one of sparse_data and
  // dense_data should be set
  bool sparse_;
  Matrix sparse_data_;
  Eigen::MatrixXd dense_data_;
  bool data_has_been_set_;
#ifdef _R_INTERFACE_
#ifdef _R_DEBUG
  /* almost uuid */
  std::string id_;
#endif
#endif
};
#endif
