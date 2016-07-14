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

#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <vector>
#include <map>
#include "Utils.hpp"
#include <iostream>

/* Stores the result of calling BUILD_MATRIX on a collection of LinOp
 * trees. */
class ProblemData {
public:
  /* COO sparse matrix representation. V stores the data, I the row indices
   * and J the column indices. */
  std::vector<double> V;
  std::vector<int> I;
  std::vector<int> J;

  /* Dense matrix representation of the constant vector */
  std::vector<double> const_vec;

  /* Map of variable_id to column in the problemData matrix */
  std::map<int, int> id_to_col;

  /* Map of constant linOp's to row in the problemData matrix  */
  std::map<int, int> const_to_row;

  /* Number of constraints */
  int num_constraints;

  /* CSC representation */
  std::vector<double> vals;
  std::vector<int> row_idxs;
  std::vector<int> col_ptrs;

  /* convert COO representation to CSC */
  void toCSC(int num_variables) {
    int nnz = I.size();
    vals.assign(nnz, 0);
    row_idxs.assign(nnz, 0);
    col_ptrs.assign(num_variables + 1, 0);
    coo_tocsc(num_constraints, num_variables, nnz,
              &I[0], &J[0], &V[0],
              &col_ptrs[0], &row_idxs[0], &vals[0]);
  }

  /*******************************************
   * The functions below return problemData vectors as contiguous 1d
   * numpy arrays.
   *
   * Note the function prototypes must match CVXCanon.i exactly to
   * properly run and compile.
   *
   * Each function is wrapped using SWIG's numpy.i typemap, so can
   * be called in python using
   *
   *        problemData.getV(N)
   *
   * where N is the length of the vector V. The double *pointer VALUES
   * is generated by the wrapper, which allocates space for NUM_VALUES
   * elements. Thus, NUM_VALUES must be exactly the length of the array.
   ********************************************/

  /**
   * Returns the data vector V as a contiguous 1D numpy array.
   */
  void getV(double* values, int num_values) {
    for (int i = 0; i < num_values; i++) {
      values[i] = V[i];
    }
  }

  /**
   * Returns the row index vector I as a contiguous 1D numpy array.
   */
  void getI(double* values, int num_values) {
    for (int i = 0; i < num_values; i++) {
      values[i] = I[i];
    }
  }

  /**
   * Returns the column index vector J as a contiguous 1D numpy array.
   */
  void getJ(double* values, int num_values) {
    for (int i = 0; i < num_values; i++) {
      values[i] = J[i];
    }
  }

  /**
   * Returns the CONST_VEC as a contiguous 1D numpy array.
   */
  void getConstVec(double* values, int num_values) {
    for (int i = 0; i < num_values; i++) {
      values[i] = const_vec[i];
    }
  }
};

#endif