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

#include "BuildMatrix.hpp"
#include <iostream>
#include <map>
#include <vector>
#include "LinOp.hpp"
#include "LinOpOperations.hpp"
#include "ProblemData.hpp"

void mul_by_const(Matrix &coeff_mat,
                  std::map<int, Matrix > &rh_coeffs,
                  std::map<int, Matrix > &result) {
  typedef std::map<int, Matrix >::iterator it_type;
  for (it_type it = rh_coeffs.begin(); it != rh_coeffs.end(); ++it) {
    int id = it->first;
    Matrix rh = it->second;
    /* Convert scalars (1x1 matrices) to primitive types */
    if (coeff_mat.rows() == 1 && coeff_mat.cols() == 1) {
      double scalar = coeff_mat.coeffRef(0, 0);
      if (result.count(id) == 0)
        result[id] = scalar * rh;
      else
        result[id] += scalar * rh;
    } else if (rh.rows() == 1 && rh.cols() == 1) {
      double scalar = rh.coeffRef(0, 0);
      if (result.count(id) == 0)
        result[id] = coeff_mat * scalar;
      else
        result[id] = coeff_mat * scalar;
    } else {
      if (result.count(id) == 0 )
        result[id] = coeff_mat * rh;
      else
        result[id] += coeff_mat * rh;
    }
  }
}

std::map<int, Matrix > get_coefficient(LinOp &lin) {
  std::map<int, Matrix > coeffs;
  if (lin.type == VARIABLE) {
    std::map<int, Matrix> new_coeffs = get_variable_coeffs(lin);
    typedef std::map<int, Matrix >::iterator it_type;
    for (it_type it = new_coeffs.begin(); it != new_coeffs.end(); ++it) {
      if (coeffs.count(it->first) == 0) {
        coeffs[it->first] = it->second;
      } else {
        coeffs[it->first] += it->second;
      }
    }
  } else if (lin.has_constant_type()) {
    /* ID will be CONSTANT_TYPE */
    std::map<int, Matrix> new_coeffs = get_const_coeffs(lin);
    typedef std::map<int, Matrix >::iterator it_type;
    for (it_type it = new_coeffs.begin(); it != new_coeffs.end(); ++it) {
      if (coeffs.count(it->first) == 0) {
        coeffs[it->first] = it->second;
      } else {
        coeffs[it->first] += it->second;
      }
    }
  } else {
    /* Multiply the arguments of the function coefficient in order */
    std::vector<Matrix> coeff_mat = get_func_coeffs(lin);
    for (unsigned i = 0; i < lin.args.size(); i++) {
      Matrix coeff = coeff_mat[i];
      std::map<int, Matrix > rh_coeffs = get_coefficient(*lin.args[i]);
      std::map<int,  Matrix > new_coeffs;
      mul_by_const(coeff, rh_coeffs, new_coeffs);

      typedef std::map<int, Matrix>::iterator it_type;
      for (it_type it = new_coeffs.begin(); it != new_coeffs.end(); ++it) {
        if (coeffs.count(it->first) == 0)
          coeffs[it->first] = it->second;
        else
          coeffs[it->first] += it->second;
      }
    }
  }
  return coeffs;
}

int get_horiz_offset(int id, std::map<int, int> &offsets,
                     int &horiz_offset, LinOp &lin) {
  if ( !offsets.count(id) ) {
    offsets[id] = horiz_offset;
    horiz_offset += lin.size[0] * lin.size[1];
  }
  return offsets[id];
}

/* function: add_matrix_to_vectors
*
* This function adds a matrix to our sparse matrix triplet
* representation, by using eigen's sparse matrix iterator
* This function takes horizontal and vertical offset, which indicate
* the offset of this block within our larger matrix.
*/
void add_matrix_to_vectors(Matrix &block, std::vector<double> &V,
                           std::vector<int>  &I, std::vector<int> &J,
                           int &vert_offset, int &horiz_offset) {
  for ( int k = 0; k < block.outerSize(); ++k ) {
    for ( Matrix::InnerIterator it(block, k); it; ++it ) {
      V.push_back(it.value());

      /* Push back current row and column indices */
      I.push_back(it.row() + vert_offset);
      J.push_back(it.col() + horiz_offset);
    }
  }
}

void extend_constant_vec(std::vector<double> &const_vec, int &vert_offset,
                         Matrix &block) {
  int rows = block.rows();
  for ( int k = 0; k < block.outerSize(); ++k ) {
    for ( Matrix::InnerIterator it(block, k); it; ++it ) {
      int idx = vert_offset + (it.col() * rows) + it.row();
      const_vec[idx] += it.value();
    }
  }
}

void process_constraint(LinOp & lin, std::vector<double> &V,
                        std::vector<int> &I, std::vector<int> &J,
                        std::vector<double> &constant_vec, int &vert_offset,
                        std::map<int, int> &id_to_col, int & horiz_offset) {
  /* Get the coefficient for the current constraint */
  std::map<int, Matrix > coeffs = get_coefficient(lin);

  typedef std::map<int, Matrix >::iterator it_type;
  for (it_type it = coeffs.begin(); it != coeffs.end(); ++it) {
    int id = it->first;                 // Horiz offset determined by the id
    Matrix block = it->second;
    if (id == CONSTANT_ID) {   // Add to CONSTANT_VEC if linop is constant
      extend_constant_vec(constant_vec, vert_offset, block);
    } else {
      int offset = get_horiz_offset(id, id_to_col, horiz_offset, lin);
      add_matrix_to_vectors(block, V, I, J, vert_offset, offset);
    }
  }
}

/* Returns the number of rows in the matrix assuming vertical stacking
   of coefficient matrices */
int get_total_constraint_length(std::vector< LinOp* > constraints) {
  int result = 0;
  for (unsigned i = 0; i < constraints.size(); i++) {
    result += constraints[i]->size[0] * constraints[i]->size[1];
  }
  return result;
}

/* Returns the number of rows in the matrix using the user provided vertical
   offsets for each constraint. */
int get_total_constraint_length(std::vector<LinOp*> &constraints,
                                std::vector<int> &constr_offsets) {
  /* Must specify an offset for each constraint */
  if (constraints.size() != constr_offsets.size()) {
    std::cerr << "Error: Invalid constraint offsets: ";
    std::cerr << "Offset must be the same length as constraints" << std::endl;
    exit(-1);
  }

  int offset_end = 0;
  /* Offsets must be monotonically increasing */
  for (unsigned i = 0; i < constr_offsets.size(); i++) {
    LinOp constr = *constraints[i];
    int offset_start = constr_offsets[i];
    offset_end = offset_start + constr.size[0] * constr.size[1];

    if (i + 1 < constr_offsets.size() && constr_offsets[i + 1] < offset_end) {
      std::cerr << "Error: Invalid constraint offsets: ";
      std::cerr << "Offsets are not monotonically increasing" << std::endl;
      exit(-1);
    }
  }
  return offset_end;
}

/* function: build_matrix
*
* Description: Given a list of linear operations, this function returns a data
* structure containing a sparse matrix representation of the cone program.
*
* Input: std::vector<LinOp *> constraints, our list of constraints represented
* as a linear operation tree
*
* Output: prob_data, a data structure which contains a sparse representation
* of the coefficient matrix, a dense representation of the constant vector,
* and maps containing our mapping from variables, and a map from the rows of our
* matrix to their corresponding constraint.
*
*/
ProblemData build_matrix(std::vector< LinOp* > constraints,
                         std::map<int, int> id_to_col) {
  ProblemData prob_data;
  int num_rows = get_total_constraint_length(constraints);
  prob_data.const_vec = std::vector<double> (num_rows, 0);
  prob_data.id_to_col = id_to_col;
  int vert_offset = 0;
  int horiz_offset  = 0;

  /* Build matrix one constraint at a time */
  for (unsigned i = 0; i < constraints.size(); i++) {
    LinOp constr = *constraints[i];
    process_constraint(constr, prob_data.V, prob_data.I, prob_data.J,
                       prob_data.const_vec, vert_offset,
                       prob_data.id_to_col, horiz_offset);
    prob_data.const_to_row[i] = vert_offset;
    vert_offset += constr.size[0] * constr.size[1];
  }
  prob_data.num_constraints = num_rows;
  return prob_data;
}

/*  See comment above for build_matrix. Requires specification of a vertical
    offset, VERT_OFFSET, for each constraint in the vector CONSTR_OFFSETS.

    Valid CONSTR_OFFSETS assume that a vertical offset is provided for each
    constraint and that the offsets are not overlapping. In particular,
    the vertical offset for constraint i + the size of constraint i must be
    less than the vertical offset for constraint i+1.
    */
ProblemData build_matrix(std::vector<LinOp*> constraints,
                         std::map<int, int> id_to_col,
                         std::vector<int> constr_offsets) {
  ProblemData prob_data;

  /* Function also verifies the offsets are valid */
  int num_rows = get_total_constraint_length(constraints, constr_offsets);
  prob_data.const_vec = std::vector<double> (num_rows, 0);
  prob_data.id_to_col = id_to_col;
  int horiz_offset  = 0;

  /* Build matrix one constraint at a time */
  for (unsigned i = 0; i < constraints.size(); i++) {
    LinOp constr = *constraints[i];
    int vert_offset = constr_offsets[i];
    process_constraint(constr, prob_data.V, prob_data.I, prob_data.J,
                       prob_data.const_vec, vert_offset,
                       prob_data.id_to_col, horiz_offset);
    prob_data.const_to_row[i] = vert_offset;
  }
  prob_data.num_constraints = num_rows;
  return prob_data;
}
