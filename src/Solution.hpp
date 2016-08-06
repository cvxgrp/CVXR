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

#ifndef SOLUTION_H
#define SOLUTION_H

#include <map>
#include "Utils.hpp"

/* Standardized solver status codes */
enum solverstatus {
  OPTIMAL,
  INFEASIBLE,
  UNBOUNDED,
  OPTIMAL_INACCURATE,
  INFEASIBLE_INACCURATE,
  UNBOUNDED_INACCURATE,
  SOLVER_ERROR
};
typedef solverstatus solverStatus;


class Solution {
public:
  /* solver STATUS */
  solverStatus status;

  double optimal_value;

  /* variable index to primal variable values */
  std::map<int, Eigen::MatrixXd> primal_values;

  /* constraint index to dual variable value.
     Note: Following CVXPY, only dual variable values for EQ and LEQ constraints are provided since
           cone constraints cannot be specified by the user
   */
  std::map<int, Eigen::MatrixXd> dual_values;
};

#endif
