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

// DPP tensor-based LinOp operations.
// Mirrors CVXPY cvxcore/src/LinOpOperations.cpp.
// All get_*_mat functions now return Tensor instead of vector<Matrix>.
// For non-DPP (no PARAM nodes), Tensor degenerates to {CONSTANT_ID: {var_id: [matrix]}}.

#include "LinOp.h"
#include "LinOpOperations.h"
#include "Utils.h"
#include <Rcpp.h>
#include <map>
#include <iostream>

/***********************
 * FUNCTION PROTOTYPES *
 ***********************/
Tensor build_tensor(Matrix &mat);
Tensor get_sum_coefficients(LinOp &lin, int arg_idx);
Tensor get_sum_entries_mat(LinOp &lin, int arg_idx);
Tensor get_trace_mat(LinOp &lin, int arg_idx);
Tensor get_neg_mat(LinOp &lin, int arg_idx);
Tensor get_div_mat(LinOp &lin, int arg_idx);
Tensor get_promote_mat(LinOp &lin, int arg_idx);
Tensor get_mul_mat(LinOp &lin, int arg_idx);
Tensor get_mul_elemwise_mat(LinOp &lin, int arg_idx);
Tensor get_rmul_mat(LinOp &lin, int arg_idx);
Tensor get_index_mat(LinOp &lin, int arg_idx);
Tensor get_transpose_mat(LinOp &lin, int arg_idx);
Tensor get_reshape_mat(LinOp &lin, int arg_idx);
Tensor get_diag_vec_mat(LinOp &lin, int arg_idx);
Tensor get_diag_matrix_mat(LinOp &lin, int arg_idx);
Tensor get_upper_tri_mat(LinOp &lin, int arg_idx);
Tensor get_conv_mat(LinOp &lin, int arg_idx);
Tensor get_hstack_mat(LinOp &lin, int arg_idx);
Tensor get_vstack_mat(LinOp &lin, int arg_idx);
Tensor get_kron_mat(LinOp &lin, int arg_idx);
Tensor get_kron_l_mat(LinOp &lin, int arg_idx);
Tensor get_variable_coeffs(LinOp &lin, int arg_idx);
Tensor get_const_coeffs(LinOp &lin, int arg_idx);
Tensor get_param_coeffs(LinOp &lin, int arg_idx);

/**
 * Dispatch to the appropriate coefficient function based on LinOp type.
 * Returns a Tensor for arg_idx.
 */
Tensor get_node_coeffs(LinOp &lin, int arg_idx) {
	Tensor coeffs;
	switch (lin.type) {
	case VARIABLE:
		coeffs = get_variable_coeffs(lin, arg_idx);
		break;
	case SCALAR_CONST:
	case DENSE_CONST:
	case SPARSE_CONST:
		coeffs = get_const_coeffs(lin, arg_idx);
		break;
	case PARAM:
		coeffs = get_param_coeffs(lin, arg_idx);
		break;
	case PROMOTE:
		coeffs = get_promote_mat(lin, arg_idx);
		break;
	case MUL_EXPR:
		coeffs = get_mul_mat(lin, arg_idx);
		break;
	case RMUL_EXPR:
		coeffs = get_rmul_mat(lin, arg_idx);
		break;
	case MUL_ELEM:
		coeffs = get_mul_elemwise_mat(lin, arg_idx);
		break;
	case DIV:
		coeffs = get_div_mat(lin, arg_idx);
		break;
	case SUM:
		coeffs = get_sum_coefficients(lin, arg_idx);
		break;
	case NEG:
		coeffs = get_neg_mat(lin, arg_idx);
		break;
	case INDEX:
		coeffs = get_index_mat(lin, arg_idx);
		break;
	case TRANSPOSE:
		coeffs = get_transpose_mat(lin, arg_idx);
		break;
	case SUM_ENTRIES:
		coeffs = get_sum_entries_mat(lin, arg_idx);
		break;
	case TRACE:
		coeffs = get_trace_mat(lin, arg_idx);
		break;
	case RESHAPE:
		coeffs = get_reshape_mat(lin, arg_idx);
		break;
	case DIAG_VEC:
		coeffs = get_diag_vec_mat(lin, arg_idx);
		break;
	case DIAG_MAT:
		coeffs = get_diag_matrix_mat(lin, arg_idx);
		break;
	case UPPER_TRI:
		coeffs = get_upper_tri_mat(lin, arg_idx);
		break;
	case CONV:
		coeffs = get_conv_mat(lin, arg_idx);
		break;
	case HSTACK:
		coeffs = get_hstack_mat(lin, arg_idx);
		break;
	case VSTACK:
		coeffs = get_vstack_mat(lin, arg_idx);
		break;
	case KRON:
		coeffs = get_kron_mat(lin, arg_idx);
		break;
	case KRON_L:
		coeffs = get_kron_l_mat(lin, arg_idx);
		break;
	default:
		Rcpp::stop("Error: linOp type invalid");
	}
	return coeffs;
}

/**
 * Recursive tensor computation — mirrors CVXPY lin_to_tensor.
 * Leaf nodes: return get_node_coeffs directly.
 * Composite nodes: multiply node coefficients with recursive arg tensors.
 */
static const int LIN_TO_TENSOR_MAX_DEPTH = 5000;

Tensor lin_to_tensor(LinOp &lin, int depth) {
	if (depth > LIN_TO_TENSOR_MAX_DEPTH) {
		Rcpp::stop("lin_to_tensor: recursion depth %d exceeds limit %d; "
		           "flatten the expression tree before canonicalization",
		           depth, LIN_TO_TENSOR_MAX_DEPTH);
	}
	if (lin.args.size() == 0) {
		return get_node_coeffs(lin, 0);
	} else {
		Tensor result;
		for (unsigned i = 0; i < lin.args.size(); ++i) {
			Tensor lh_coeff = get_node_coeffs(lin, i);
			Tensor rh_coeff = lin_to_tensor(*lin.args[i], depth + 1);
			Tensor prod = tensor_mul(lh_coeff, rh_coeff);
			acc_tensor(result, prod);
		}
		return result;
	}
}

/*******************
 * HELPER FUNCTIONS
 *******************/

/**
 * Wraps a single Matrix into a Tensor: {CONSTANT_ID: {CONSTANT_ID: [mat]}}
 */
Tensor build_tensor(Matrix &mat) {
	Tensor ten;
	ten[CONSTANT_ID] = DictMat();
	DictMat* dm = &(ten[CONSTANT_ID]);
	(*dm)[CONSTANT_ID] = std::vector<Matrix>();

	std::vector<Matrix>* mat_vec = &(*dm)[CONSTANT_ID];
	mat_vec->push_back(Matrix());
	(*mat_vec)[0].swap(mat);
	return ten;
}

/**
 * Returns an N x N sparse identity matrix.
 */
Matrix sparse_eye (int n) {
	Matrix eye_n(n, n);
	eye_n.setIdentity();
	return eye_n;
}

/**
 * Returns a sparse ROWS x COLS matrix of all ones.
 */
Matrix sparse_ones(int rows, int cols) {
	Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(rows, cols);
	Matrix mat = ones.sparseView();
	return mat;
}

/**
 * Returns a sparse rows x cols matrix with matrix[row_sel, col_sel] = 1.
 * Used for PARAM coefficient construction.
 */
Matrix sparse_selector(int rows, int cols, int row_sel, int col_sel) {
	Matrix selector(rows * cols, 1);
	selector.insert(row_sel + rows * col_sel, 0) = 1.0;
	return selector;
}

/**
 * Reshapes the input matrix into a single column vector.
 */
Matrix sparse_reshape_to_vec(Matrix &mat) {
	int rows = mat.rows();
	int cols = mat.cols();
	Matrix out(rows * cols, 1);
	std::vector<Triplet> tripletList;
	tripletList.reserve(rows * cols);
	for ( int k = 0; k < mat.outerSize(); ++k) {
		for (Matrix::InnerIterator it(mat, k); it; ++it) {
			tripletList.push_back(Triplet(it.col() * rows + it.row(), 0,
			                              it.value()));
		}
	}
	out.setFromTriplets(tripletList.begin(), tripletList.end());
	out.makeCompressed();
	return out;
}

/******************
 * Data access helpers
 ******************/

/**
 * Returns the constant data from a LinOp as a sparse matrix.
 */
Matrix get_constant_data(LinOp &lin, bool column) {
	Matrix coeffs;
	if (lin.sparse) {
		if (column) {
			coeffs = sparse_reshape_to_vec(lin.sparse_data);
		} else {
			coeffs = lin.sparse_data;
		}
	} else {
		if (column) {
			Eigen::Map<Eigen::MatrixXd> column(lin.dense_data.data(),
			                                   lin.dense_data.rows() *
			                                   lin.dense_data.cols(), 1);
			coeffs = column.sparseView();
		} else {
			coeffs = lin.dense_data.sparseView();
		}
	}
	coeffs.makeCompressed();
	return coeffs;
}

/**
 * Interface for the INDEX linOp to retrieve slice data.
 */
std::vector<std::vector<int> > get_slice_data(LinOp &lin, int rows, int cols) {
	if (lin.type != INDEX) Rcpp::stop("get_slice_data: expected INDEX LinOp type");
	std::vector<int> row_slice = lin.slice[0];
	std::vector<int> col_slice = lin.slice[1];
	if (row_slice.size() != 3) Rcpp::stop("get_slice_data: row_slice must have 3 elements (start, end, step)");
	if (col_slice.size() != 3) Rcpp::stop("get_slice_data: col_slice must have 3 elements (start, end, step)");

	std::vector<std::vector<int> > slices;
	slices.push_back(row_slice);
	slices.push_back(col_slice);
	return slices;
}

/**
 * Returns the divisor from a DIV LinOp.
 */
double get_divisor_data(LinOp &lin) {
	if (lin.type != DIV) Rcpp::stop("get_divisor_data: expected DIV LinOp type");
	return lin.dense_data(0, 0);
}

/**
 * Returns the ID from a VARIABLE or PARAM LinOp.
 */
int get_id_data(LinOp &lin) {
	if (lin.type != VARIABLE && lin.type != PARAM)
		Rcpp::stop("get_id_data: expected VARIABLE or PARAM LinOp type");
	return int(lin.dense_data(0, 0));
}

/*****************************
 * LinOP -> Tensor FUNCTIONS
 *****************************/

/**
 * KRON: variable on right, constant on left.
 */
Tensor get_kron_mat(LinOp &lin, int arg_idx) {
	if (lin.type != KRON) Rcpp::stop("get_kron_mat: expected KRON LinOp type");
	Matrix constant = get_constant_data(lin, false);
	int lh_rows = constant.rows();
	int lh_cols = constant.cols();
	int rh_rows =  lin.args[0]->size[0];
	int rh_cols =  lin.args[0]->size[1];

	int rows = rh_rows * rh_cols * lh_rows * lh_cols;
	int cols = rh_rows * rh_cols;
	Matrix coeffs(rows, cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(rh_rows * rh_cols * constant.nonZeros());
	for ( int k = 0; k < constant.outerSize(); ++k ) {
		for ( Matrix::InnerIterator it(constant, k); it; ++it ) {
			int row = (rh_rows * rh_cols * (lh_rows * it.col())) + (it.row() * rh_rows);
			int col = 0;
			for(int j = 0; j < rh_cols; j++){
				for(int i = 0; i < rh_rows; i++) {
					tripletList.push_back(Triplet(row + i, col, it.value()));
					col++;
				}
				row += lh_rows * rh_rows;
			}
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * KRON_L: variable on left, constant on right.
 */
Tensor get_kron_l_mat(LinOp &lin, int arg_idx) {
	if (lin.type != KRON_L) Rcpp::stop("get_kron_l_mat: expected KRON_L LinOp type");
	Matrix rh = get_constant_data(lin, false);
	int rh_rows = rh.rows();
	int rh_cols = rh.cols();
	int lh_rows = lin.args[0]->size[0];
	int lh_cols = lin.args[0]->size[1];

	int kron_rows = lh_rows * rh_rows;
	int rh_nnz = rh.nonZeros();
	std::vector<int> base_row_indices;
	std::vector<double> vec_rh;
	base_row_indices.reserve(rh_nnz);
	vec_rh.reserve(rh_nnz);
	int row_offset = 0;
	for (int k = 0; k < rh.outerSize(); ++k) {
		for (Matrix::InnerIterator it(rh, k); it; ++it) {
			int cur_row = it.row() + row_offset;
			base_row_indices.push_back(cur_row);
			vec_rh.push_back(it.value());
		}
		row_offset += kron_rows;
	}

	int lh_size = lh_rows * lh_cols;
	int rh_size = rh_rows * rh_cols;
	Matrix mat(lh_size * rh_size, lh_size);
	std::vector<Triplet> tripletList;
	tripletList.reserve(lh_size * rh_nnz);

	int row, col;
	int outer_row_offset = 0;
	for (int j = 0; j < lh_cols; ++j) {
		row_offset = outer_row_offset;
		for (int i = 0; i < lh_rows; ++i) {
			col = i + j * lh_rows;
			for (int ell = 0; ell < rh_nnz; ++ell) {
				row = base_row_indices[ell] + row_offset;
				tripletList.push_back(Triplet(row, col, vec_rh[ell]));
			}
			row_offset += rh_rows;
		}
		outer_row_offset += lh_rows * rh_size;
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	mat.makeCompressed();
	return build_tensor(mat);
}

/**
 * VSTACK: vertical stacking. Returns coefficient for arg_idx only.
 */
Tensor get_vstack_mat(LinOp &lin, int arg_idx) {
	if (lin.type != VSTACK) Rcpp::stop("get_vstack_mat: expected VSTACK LinOp type");

	/* Compute row offset for this specific arg */
	int row_offset = 0;
	for (int idx = 0; idx < arg_idx; ++idx) {
		row_offset += lin.args[idx]->size[0];
	}

	LinOp &arg = *lin.args[arg_idx];
	int arg_rows = arg.size[0];
	int arg_cols = arg.size[1];
	int column_offset = lin.size[0];

	std::vector<Triplet> tripletList;
	tripletList.reserve(arg_rows * arg_cols);
	for (int i = 0; i < arg_rows; i++) {
		for (int j = 0; j < arg_cols; j++) {
			int row_idx = i + (j * column_offset) + row_offset;
			int col_idx = i + (j * arg_rows);
			tripletList.push_back(Triplet(row_idx, col_idx, 1));
		}
	}

	Matrix coeff(lin.size[0] * lin.size[1], arg_rows * arg_cols);
	coeff.setFromTriplets(tripletList.begin(), tripletList.end());
	coeff.makeCompressed();
	return build_tensor(coeff);
}

/**
 * HSTACK: horizontal stacking. Returns coefficient for arg_idx only.
 */
Tensor get_hstack_mat(LinOp &lin, int arg_idx) {
	if (lin.type != HSTACK) Rcpp::stop("get_hstack_mat: expected HSTACK LinOp type");

	/* Compute row offset for this specific arg */
	int row_offset = 0;
	for (int idx = 0; idx < arg_idx; ++idx) {
		LinOp &prev_arg = *lin.args[idx];
		row_offset += prev_arg.size[0] * prev_arg.size[1];
	}

	LinOp &arg = *lin.args[arg_idx];
	int arg_rows = arg.size[0];
	int arg_cols = arg.size[1];
	int column_offset = arg_rows;

	std::vector<Triplet> tripletList;
	tripletList.reserve(arg_rows * arg_cols);
	for (int i = 0; i < arg_rows; i++) {
		for (int j = 0; j < arg_cols; j++) {
			int row_idx = i + (j * column_offset) + row_offset;
			int col_idx = i + (j * arg_rows);
			tripletList.push_back(Triplet(row_idx, col_idx, 1));
		}
	}

	Matrix coeff(lin.size[0] * lin.size[1], arg_rows * arg_cols);
	coeff.setFromTriplets(tripletList.begin(), tripletList.end());
	coeff.makeCompressed();
	return build_tensor(coeff);
}

/**
 * CONV: convolution via toeplitz matrix.
 */
Tensor get_conv_mat(LinOp &lin, int arg_idx) {
	if (lin.type != CONV) Rcpp::stop("get_conv_mat: expected CONV LinOp type");
	Matrix constant = get_constant_data(lin, false);
	int rows = lin.size[0];
	int nonzeros = constant.rows();
	int cols = lin.args[0]->size[0];

	Matrix toeplitz(rows, cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(nonzeros * cols);
	for (int col = 0; col < cols; col++) {
		int row_start = col;
		for ( int k = 0; k < constant.outerSize(); ++k ) {
			for ( Matrix::InnerIterator it(constant, k); it; ++it ) {
				int row_idx = row_start + it.row();
				tripletList.push_back(Triplet(row_idx, col, it.value()));
			}
		}
	}
	toeplitz.setFromTriplets(tripletList.begin(), tripletList.end());
	toeplitz.makeCompressed();
	return build_tensor(toeplitz);
}

/**
 * UPPER_TRI: extract upper triangular elements.
 */
Tensor get_upper_tri_mat(LinOp &lin, int arg_idx) {
	if (lin.type != UPPER_TRI) Rcpp::stop("get_upper_tri_mat: expected UPPER_TRI LinOp type");
	int rows = lin.args[0]->size[0];
	int cols = lin.args[0]->size[1];

	int entries = lin.size[0];
	Matrix coeffs(entries, rows * cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(entries);
	int count = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (j > i) {
				int row_idx = count;
				count++;
				int col_idx = j * rows + i;
				tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
			}
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * DIAG_MAT: diagonal matrix to vector (with offset k).
 */
Tensor get_diag_matrix_mat(LinOp &lin, int arg_idx) {
	if (lin.type != DIAG_MAT) Rcpp::stop("get_diag_matrix_mat: expected DIAG_MAT LinOp type");
	int n = lin.size[0];
	int k = (lin.dense_data.size() > 0) ? int(lin.dense_data(0, 0)) : 0;
	int ak = (k >= 0) ? k : -k;
	int m = n + ak;

	Matrix coeffs(n, m * m);
	std::vector<Triplet> tripletList;
	tripletList.reserve(n);
	for (int i = 0; i < n; i++) {
		int col_idx;
		if (k >= 0) {
			col_idx = (i + k) * m + i;
		} else {
			col_idx = i * m + (i + ak);
		}
		tripletList.push_back(Triplet(i, col_idx, 1.0));
	}

	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * DIAG_VEC: vector to diagonal matrix (with offset k).
 */
Tensor get_diag_vec_mat(LinOp &lin, int arg_idx) {
	if (lin.type != DIAG_VEC) Rcpp::stop("get_diag_vec_mat: expected DIAG_VEC LinOp type");
	int M = lin.size[0];
	int k = (lin.dense_data.size() > 0) ? int(lin.dense_data(0, 0)) : 0;
	int ak = (k >= 0) ? k : -k;
	int n = M - ak;

	Matrix coeffs(M * M, n);
	std::vector<Triplet> tripletList;
	tripletList.reserve(n);
	for (int i = 0; i < n; i++) {
		int row_idx;
		if (k >= 0) {
			row_idx = (i + k) * M + i;
		} else {
			row_idx = i * M + (i + ak);
		}
		tripletList.push_back(Triplet(row_idx, i, 1.0));
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * TRANSPOSE: permutation matrix.
 */
Tensor get_transpose_mat(LinOp &lin, int arg_idx) {
	if (lin.type != TRANSPOSE) Rcpp::stop("get_transpose_mat: expected TRANSPOSE LinOp type");
	int rows = lin.size[0];
	int cols = lin.size[1];

	Matrix coeffs(rows * cols, rows * cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(rows * cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int row_idx = rows * j + i;
			int col_idx = i * cols + j;
			tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * INDEX: selection matrix.
 */
Tensor get_index_mat(LinOp &lin, int arg_idx) {
	if (lin.type != INDEX) Rcpp::stop("get_index_mat: expected INDEX LinOp type");
	int rows = lin.args[0]->size[0];
	int cols = lin.args[0]->size[1];
	Matrix coeffs (lin.size[0] * lin.size[1], rows * cols);

	if (coeffs.rows () == 0 ||  coeffs.cols() == 0) {
		return build_tensor(coeffs);
	}

#ifdef _R_INTERFACE_
	std::vector<Triplet> tripletList;
	std::vector<int> col_slice = lin.slice[1];
	std::vector<int> row_slice = lin.slice[0];
	int counter = 0;
	for (unsigned j = 0; j < col_slice.size(); j++) {
	  for (unsigned i = 0; i < row_slice.size(); i++) {
	    int row_idx = counter;
	    int col_idx = col_slice[j] * rows + row_slice[i];
	    tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
	    counter++;
	  }
	}
#else
	std::vector<std::vector<int> > slices = get_slice_data(lin, rows, cols);
	int row_start = slices[0][0];
	int row_end = slices[0][1];
	int row_step = slices[0][2];
	int col_start = slices[1][0];
	int col_end = slices[1][1];
	int col_step = slices[1][2];

	std::vector<Triplet> tripletList;
	int col = col_start;
	int counter = 0;
	while (true) {
		if (col < 0 || col >= cols) break;
		int row = row_start;
		while (true) {
			if (row < 0 || row >= rows) break;
			int row_idx = counter;
			int col_idx = col * rows + row;
			tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
			counter++;
			row += row_step;
			if ((row_step > 0 && row >= row_end) || (row_step < 0 && row <= row_end)) break;
		}
		col += col_step;
		if ((col_step > 0 && col >= col_end) || (col_step < 0 && col <= col_end)) break;
	}
#endif
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * MUL_ELEM: element-wise multiplication (diagonal matrix).
 * For non-DPP, wraps constant data in build_tensor then diagonalizes.
 */
Tensor get_mul_elemwise_mat(LinOp &lin, int arg_idx) {
	if (lin.type != MUL_ELEM) Rcpp::stop("get_mul_elemwise_mat: expected MUL_ELEM LinOp type");

	/* DPP path: use lin_to_tensor on data sub-tree (matches CVXPY).
	 * This handles both constant and parametric data uniformly. */
	if (lin.linOp_data != nullptr) {
		Tensor mul_ten = lin_to_tensor(*lin.linOp_data);
		/* Diagonalize all matrices in the tensor */
		for (auto it = mul_ten.begin(); it != mul_ten.end(); ++it) {
			int param_id = it->first;
			DictMat &var_map = it->second;
			for (auto jit = var_map.begin(); jit != var_map.end(); ++jit) {
				int var_id = jit->first;
				std::vector<Matrix> &mat_vec = jit->second;
				for (unsigned i = 0; i < mat_vec.size(); ++i) {
					mat_vec[i] = diagonalize(mat_vec[i]);
				}
			}
		}
		return mul_ten;
	}

	/* Non-DPP fallback: inline constant data */
	Matrix constant = get_constant_data(lin, true);
	int n = constant.rows();

	std::vector<Triplet> tripletList;
	tripletList.reserve(n);
	for ( int k = 0; k < constant.outerSize(); ++k ) {
		for ( Matrix::InnerIterator it(constant, k); it; ++it ) {
			tripletList.push_back(Triplet(it.row(), it.row(), it.value()));
		}
	}
	Matrix coeffs(n, n);
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * RMUL_EXPR: right multiplication (kronecker product with identity).
 */
Tensor get_rmul_mat(LinOp &lin, int arg_idx) {
	if (lin.type != RMUL_EXPR) Rcpp::stop("get_rmul_mat: expected RMUL_EXPR LinOp type");

	/* DPP path: use lin_to_tensor on data sub-tree.
	 * lin_to_tensor returns flattened column vectors for constant data,
	 * so we must use modular arithmetic to un-flatten indices.
	 * Matches CVXPY LinOpOperations.cpp lines 801-892. */
	if (lin.linOp_data != nullptr) {
		Tensor rmul_ten = lin_to_tensor(*lin.linOp_data);

		/* Interpret as row or column vector as needed (CVXPY lines 808-822).
		 * When data_ndim==1 and first arg is not a row vector, transpose
		 * all coefficient matrices. In R, shapes are always 2D so
		 * data_ndim is always >= 2, making this a no-op in practice. */
		if (lin.data_ndim == 1 && lin.args[0]->size[0] != 1) {
			for (auto it = rmul_ten.begin(); it != rmul_ten.end(); ++it) {
				DictMat &var_map = it->second;
				for (auto jit = var_map.begin(); jit != var_map.end(); ++jit) {
					std::vector<Matrix> &mat_vec = jit->second;
					for (unsigned i = 0; i < mat_vec.size(); ++i) {
						mat_vec[i] = mat_vec[i].transpose();
					}
				}
			}
		}

		/* Get data shape. Matches CVXPY lines 823-840. */
		int data_rows = (lin.linOp_data->size.size() >= 1)
			? lin.linOp_data->size[0] : 1;
		int data_cols = (lin.linOp_data->size.size() >= 2)
			? lin.linOp_data->size[1] : 1;
		/* n = result_rows: for 2D shapes, equivalent to lin.size[0]
		 * (CVXPY uses 3-way shape analysis for 0D/1D/2D args) */
		int n = (lin.args[0]->size.size() >= 2)
			? lin.args[0]->size[0] : 1;
		/* Kronecker product with identity for each tensor matrix */
		for (auto it = rmul_ten.begin(); it != rmul_ten.end(); ++it) {
			int param_id = it->first;
			DictMat &var_map = it->second;
			for (auto jit = var_map.begin(); jit != var_map.end(); ++jit) {
				int var_id = jit->first;
				std::vector<Matrix> &mat_vec = jit->second;
				for (unsigned i = 0; i < mat_vec.size(); ++i) {
					const Matrix &curr_matrix = mat_vec[i];
					Matrix coeffs(data_cols * n, data_rows * n);
					std::vector<Triplet> tripletList;
					tripletList.reserve(n * curr_matrix.nonZeros());
					for (int k = 0; k < curr_matrix.outerSize(); ++k) {
						for (Matrix::InnerIterator mit(curr_matrix, k); mit; ++mit) {
							double val = mit.value();
							/* Data is flattened — un-flatten using
							 * modular arithmetic (CVXPY lines 859-868) */
							int col, row;
							if (curr_matrix.rows() == 1) {
								col = mit.col() / data_rows;
								row = mit.col() % data_rows;
							} else {
								col = mit.row() / data_rows;
								row = mit.row() % data_rows;
							}
							int row_start = col * n;
							int col_start = row * n;
							for (int ii = 0; ii < n; ii++) {
								tripletList.push_back(Triplet(row_start + ii, col_start + ii, val));
							}
						}
					}
					coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
					coeffs.makeCompressed();
					mat_vec[i] = coeffs;
				}
			}
		}
		return rmul_ten;
	}

	/* Non-DPP fallback: inline constant data */
	Matrix constant = get_constant_data(lin, false);
	int rows = constant.rows();
	int cols = constant.cols();
	int n = lin.size[0];

	Matrix coeffs(cols * n, rows * n);
	std::vector<Triplet> tripletList;
	tripletList.reserve(n * constant.nonZeros());
	for ( int k = 0; k < constant.outerSize(); ++k ) {
		for ( Matrix::InnerIterator it(constant, k); it; ++it ) {
			double val = it.value();
			int row_start = it.col() * n;
			int col_start = it.row() * n;
			for (int i = 0; i < n; i++) {
				int row_idx = row_start + i;
				int col_idx = col_start + i;
				tripletList.push_back(Triplet(row_idx, col_idx, val));
			}
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * MUL_EXPR: left multiplication (block diagonal matrix).
 */
Tensor get_mul_mat(LinOp &lin, int arg_idx) {
	if (lin.type != MUL_EXPR) Rcpp::stop("get_mul_mat: expected MUL_EXPR LinOp type");

	/* DPP path: use linOp_data sub-tree.
	 * Matches CVXPY LinOpOperations.cpp lines 904-1007.
	 * Key: for SPARSE_CONST/DENSE_CONST data, use get_constant_data
	 * (keeps 2D shape) with data_flattened flag. For parametric data,
	 * use lin_to_tensor (always flattened). */
	if (lin.linOp_data != nullptr) {
		/* Get original data shape (CVXPY lines 908-914) */
		int data_rows = lin.linOp_data->size[0];
		int data_cols = (lin.linOp_data->size.size() >= 2)
			? lin.linOp_data->size[1] : 1;
		/* num_blocks from arg shape (CVXPY line 916-918).
		 * For 2D shapes, equivalent to lin.size[1]. */
		int num_blocks = (lin.args[0]->size.size() <= 1)
			? 1 : lin.args[0]->size[1];
		/* Determine block dimensions — may need swapping (CVXPY lines 920-925) */
		int block_rows = data_rows;
		int block_cols = data_cols;
		if (lin.args[0]->size[0] != data_cols) {
			block_rows = data_cols;
			block_cols = data_rows;
		}

		/* Get tensor — fast constant path or lin_to_tensor
		 * (matches CVXPY lines 928-939) */
		Tensor mul_ten;
		bool data_flattened = true;
		if (lin.linOp_data->type == SPARSE_CONST ||
		    lin.linOp_data->type == DENSE_CONST) {
			/* Fast path: keep 2D shape */
			data_flattened = (data_rows == 1 || data_cols == 1);
			Matrix coeffs = get_constant_data(*lin.linOp_data, false);
			mul_ten = build_tensor(coeffs);
		} else {
			/* Parametric or composite: lin_to_tensor flattens */
			mul_ten = lin_to_tensor(*lin.linOp_data);
		}

		/* Interpret as row or column vector as needed (CVXPY lines 941-956).
		 * When data_ndim==1 and first arg is not a row vector, transpose
		 * all coefficient matrices. In R, shapes are always 2D so
		 * data_ndim is always >= 2, making this a no-op in practice. */
		if (lin.data_ndim == 1 && lin.args[0]->size[0] != 1) {
			for (auto it = mul_ten.begin(); it != mul_ten.end(); ++it) {
				DictMat &var_map = it->second;
				for (auto jit = var_map.begin(); jit != var_map.end(); ++jit) {
					std::vector<Matrix> &mat_vec = jit->second;
					for (unsigned i = 0; i < mat_vec.size(); ++i) {
						mul_ten[it->first][jit->first][i] =
							mat_vec[i].transpose();
					}
				}
			}
		}

		/* For scalars, just return the tensor directly */
		if (block_rows == 1 && block_cols == 1) {
			return mul_ten;
		}

		/* Block-diagonalize all matrices in the tensor
		 * (matches CVXPY lines 960-1005) */
		for (auto it = mul_ten.begin(); it != mul_ten.end(); ++it) {
			int param_id = it->first;
			DictMat &var_map = it->second;
			for (auto jit = var_map.begin(); jit != var_map.end(); ++jit) {
				int var_id = jit->first;
				std::vector<Matrix> &mat_vec = jit->second;
				for (unsigned i = 0; i < mat_vec.size(); ++i) {
					Matrix block_diag(num_blocks * block_rows,
					                  num_blocks * block_cols);
					std::vector<Triplet> tripletList;
					tripletList.reserve(num_blocks * mat_vec[i].nonZeros());
					for (int curr_block = 0; curr_block < num_blocks;
					     curr_block++) {
						int start_i = curr_block * block_rows;
						int start_j = curr_block * block_cols;
						const Matrix &curr_matrix = mat_vec[i];
						for (int k = 0; k < curr_matrix.outerSize();
						     ++k) {
							for (Matrix::InnerIterator mit(
							         curr_matrix, k);
							     mit; ++mit) {
								/* Un-flatten when data was
								 * flattened by lin_to_tensor
								 * (CVXPY lines 979-993) */
								int row, col;
								if (data_flattened) {
									if (curr_matrix.rows()
									    == 1) {
										row = mit.col()
										      % block_rows;
										col = mit.col()
										      / block_rows;
									} else {
										row = mit.row()
										      % block_rows;
										col = mit.row()
										      / block_rows;
									}
								} else {
									row = mit.row();
									col = mit.col();
								}
								tripletList.push_back(
								    Triplet(start_i + row,
								            start_j + col,
								            mit.value()));
							}
						}
					}
					block_diag.setFromTriplets(
					    tripletList.begin(), tripletList.end());
					block_diag.makeCompressed();
					mat_vec[i] = block_diag;
				}
			}
		}
		return mul_ten;
	}

	/* Non-DPP fallback: inline constant data */
	Matrix block = get_constant_data(lin, false);
	int block_rows = block.rows();
	int block_cols = block.cols();

	// Don't replicate scalars
	if(block_rows == 1 && block_cols == 1){
		return build_tensor(block);
	}

	int num_blocks = lin.size[1];
	Matrix coeffs (num_blocks * block_rows, num_blocks * block_cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(num_blocks * block.nonZeros());
	for (int curr_block = 0; curr_block < num_blocks; curr_block++) {
		int start_i = curr_block * block_rows;
		int start_j = curr_block * block_cols;
		for ( int k = 0; k < block.outerSize(); ++k ) {
			for ( Matrix::InnerIterator it(block, k); it; ++it ) {
				tripletList.push_back(Triplet(start_i + it.row(), start_j + it.col(),
				                              it.value()));
			}
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * PROMOTE: column of ones.
 */
Tensor get_promote_mat(LinOp &lin, int arg_idx) {
	if (lin.type != PROMOTE) Rcpp::stop("get_promote_mat: expected PROMOTE LinOp type");
	int num_entries = lin.size[0] * lin.size[1];
	Matrix ones = sparse_ones(num_entries, 1);
	ones.makeCompressed();
	return build_tensor(ones);
}

/**
 * RESHAPE: identity matrix of size n (preserves all elements).
 */
Tensor get_reshape_mat(LinOp &lin, int arg_idx) {
	if (lin.type != RESHAPE) Rcpp::stop("get_reshape_mat: expected RESHAPE LinOp type");
	int n = lin.size[0] * lin.size[1];
	Matrix coeffs = sparse_eye(n);
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * DIV: diagonal matrix with 1/divisor entries.
 */
Tensor get_div_mat(LinOp &lin, int arg_idx) {
	if (lin.type != DIV) Rcpp::stop("get_div_mat: expected DIV LinOp type");
	double divisor = get_divisor_data(lin);
	int n = lin.size[0] * lin.size[1];
	Matrix coeffs = sparse_eye(n);
	coeffs /= divisor;
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * NEG: negative identity.
 */
Tensor get_neg_mat(LinOp &lin, int arg_idx) {
	if (lin.type != NEG) Rcpp::stop("get_neg_mat: expected NEG LinOp type");
	int n = lin.size[0] * lin.size[1];
	Matrix coeffs = sparse_eye(n);
	coeffs *= -1;
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * TRACE: row vector extracting diagonal entries.
 */
Tensor get_trace_mat(LinOp &lin, int arg_idx) {
	if (lin.type != TRACE) Rcpp::stop("get_trace_mat: expected TRACE LinOp type");
	int rows = lin.args[0]->size[0];
	Matrix coeffs (1, rows * rows);
	for (int i = 0; i < rows; i++) {
		coeffs.insert(0, i * rows + i) = 1;
	}
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * SUM_ENTRIES: row vector of ones.
 */
Tensor get_sum_entries_mat(LinOp &lin, int arg_idx) {
	if (lin.type != SUM_ENTRIES) Rcpp::stop("get_sum_entries_mat: expected SUM_ENTRIES LinOp type");
	int rows = lin.args[0]->size[0];
	int cols = lin.args[0]->size[1];
	Matrix coeffs = sparse_ones(1, rows * cols);
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * SUM: identity matrix (same coefficient for each arg).
 * In the CVXPY approach, each arg gets the same identity coefficient
 * which when accumulated produces the sum.
 */
Tensor get_sum_coefficients(LinOp &lin, int arg_idx) {
	if (lin.type != SUM) Rcpp::stop("get_sum_coefficients: expected SUM LinOp type");
	int n = lin.size[0] * lin.size[1];
	Matrix coeffs = sparse_eye(n);
	coeffs.makeCompressed();
	return build_tensor(coeffs);
}

/**
 * VARIABLE: identity matrix keyed by variable ID.
 */
Tensor get_variable_coeffs(LinOp &lin, int arg_idx) {
	if (lin.type != VARIABLE) Rcpp::stop("get_variable_coeffs: expected VARIABLE LinOp type");
	int id = get_id_data(lin);

	Tensor ten;
	DictMat &id_to_coeffs = ten[CONSTANT_ID];
	std::vector<Matrix> &mat_vec = id_to_coeffs[id];

	int n = lin.size[0] * lin.size[1];
	Matrix coeffs = sparse_eye(n);
	coeffs.makeCompressed();
	mat_vec.push_back(Matrix());
	mat_vec[0].swap(coeffs);

	return ten;
}

/**
 * PARAM: selector matrices keyed by parameter ID.
 * Each element of the parameter gets a sparse selector matrix.
 */
Tensor get_param_coeffs(LinOp &lin, int arg_idx) {
	if (lin.type != PARAM) Rcpp::stop("get_param_coeffs: expected PARAM LinOp type");
	int id = get_id_data(lin);
	unsigned m = lin.size[0];
	unsigned n = lin.size[1];

	Tensor ten;
	DictMat &dm = ten[id];
	std::vector<Matrix> &mat_vec = dm[CONSTANT_ID];

	// Make mxn matrices with one nonzero each,
	// stack them in column major order.
	for (unsigned j = 0; j < n; ++j) {
		for (unsigned i = 0; i < m; ++i) {
			mat_vec.push_back(sparse_selector(m, n, i, j));
		}
	}
	return ten;
}

/**
 * CONSTANT: column vector keyed by CONSTANT_ID.
 */
Tensor get_const_coeffs(LinOp &lin, int arg_idx) {
	if (!lin.has_constant_type()) Rcpp::stop("get_const_coeffs: expected constant LinOp type");
	Tensor ten;
	DictMat &id_to_coeffs = ten[CONSTANT_ID];
	std::vector<Matrix> &mat_vec = id_to_coeffs[CONSTANT_ID];

	Matrix coeffs = get_constant_data(lin, true);
	coeffs.makeCompressed();
	mat_vec.push_back(Matrix());
	mat_vec[0].swap(coeffs);

	return ten;
}
