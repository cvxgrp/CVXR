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

// Some useful defines for Matricies/etc.

#ifndef UTILS_H
#define UTILS_H

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/lexical_cast.hpp>         // lexical cast
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

#include <Rcpp.h>
#include <RcppEigen.h>

#include <algorithm>
#include <set>
#define NULL_MATRIX Eigen::SparseMatrix<double>(0,0)

typedef Eigen::Matrix<int, Eigen::Dynamic, 1> Vector;
typedef Eigen::SparseMatrix<double> Matrix;
typedef std::map<int, Matrix> CoeffMap;
typedef Eigen::Triplet<double> Triplet;

// Taken from
// https://github.com/pv/scipy-work/blob/master/scipy/sparse/sparsetools/coo.h
/*
 * Compute B = A for COO matrix A, CSR matrix B
 *
 *
 * Input Arguments:
 *   I  n_row      - number of rows in A
 *   I  n_col      - number of columns in A
 *   I  nnz        - number of nonzeros in A
 *   I  Ai[nnz(A)] - row indices
 *   I  Aj[nnz(A)] - column indices
 *   T  Ax[nnz(A)] - nonzeros
 * Output Arguments:
 *   I Bp  - row pointer
 *   I Bj  - column indices
 *   T Bx  - nonzeros
 *
 * Note:
 *   Output arrays Bp, Bj, and Bx must be preallocated
 *
 * Note:
 *   Input:  row and column indices *are not* assumed to be ordered
 *
 *   Note: duplicate entries are carried over to the CSR represention
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 *
 */
template <class I, class T>
void coo_tocsr(const I n_row,
               const I n_col,
               const I nnz,
               const I Ai[],
               const I Aj[],
               const T Ax[],
               I Bp[],
               I Bj[],
               T Bx[])
{
  //compute number of non-zero entries per row of A
  std::fill(Bp, Bp + n_row, 0);

  for (I n = 0; n < nnz; n++) {
    Bp[Ai[n]]++;
  }

  //cumsum the nnz per row to get Bp[]
  for (I i = 0, cumsum = 0; i < n_row; i++) {
    I temp = Bp[i];
    Bp[i] = cumsum;
    cumsum += temp;
  }
  Bp[n_row] = nnz;

  //write Aj,Ax into Bj,Bx
  for (I n = 0; n < nnz; n++) {
    I row  = Ai[n];
    I dest = Bp[row];

    Bj[dest] = Aj[n];
    Bx[dest] = Ax[n];

    Bp[row]++;
  }

  for (I i = 0, last = 0; i <= n_row; i++) {
    I temp = Bp[i];
    Bp[i]  = last;
    last   = temp;
  }

  //now Bp,Bj,Bx form a CSR representation (with possible duplicates)
}

template<class I, class T>
void coo_tocsc(const I n_row,
               const I n_col,
               const I nnz,
               const I Ai[],
               const I Aj[],
               const T Ax[],
               I Bp[],
               I Bi[],
               T Bx[])
{ coo_tocsr<I, T>(n_col, n_row, nnz, Aj, Ai, Ax, Bp, Bi, Bx); }

#endif
