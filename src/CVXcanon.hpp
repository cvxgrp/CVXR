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

#ifndef CVXCANON_H
#define CVXCANON_H

#include <vector>
#include "LinOp.hpp"
#include "Utils.hpp"
#include "Solution.hpp"

enum objectivesense {
  MINIMIZE,
  MAXIMIZE
};
typedef objectivesense Sense;

// Top Level Entry point
/**
 * Arguments:
 *  Sense: MAXIMIZE or MINIMIZE
 *
 *  Objective: LinOp expression tree representing the OBJECTIVE
 *
 *  Constraints: List of LinOp expression tree for each constraint
 *    Constraint expression trees are represented as LinOp nodes.
 *
 *    Constraint TYPE:  EQ,     // equality constraint
 *                      LEQ,    // non-negative orthant
 *                      SOC,    // second-order cone
 *                      EXP,    // exponential cone
 *                      SDP,    // semi-definite cone **** NOT CURRENTLY SUPPORTED
 *
 *    SIZE: constraint SIZE
 *
 *    DENSE_DATA: constraint ID
 *
 *    ARGUMENTS: type EQ, LEQ: ARGUMENTS[0] is the linOp expression tree
 *
 *               type SOC:    ARGUMENTS[0] is the expression tree for T
 *                            ARGUMENTS[1:N] are expression trees for X_ELEMS
 *
 *               type EXP:    ARGUMENTS[0] expression tree for X
 *                            ARGUMENTS[1] expression tree for Y
 *                            ARGUMENTS[2] expression tree for Z
 *
 *
 *  Solver options:
 *    feastol: the tolerance on the primal and dual residual
 *    abstol: the absolute tolerance on the duality gap
 *    reltol: the relative tolerance on the duality gap
 *    feastol_inacc: the tolerance on the primal and dual residual if reduced precisions
 *    abstol_inacc: the absolute tolerance on the duality gap if reduced precision
 *    reltol_inacc: the relative tolerance on the duality gap if reduced precision
 *    max_iters: the maximum numer of iterations
 *    verbose: signals to print on non zero value
 *
 *
 * Returns:
 *  SOLUTION objective for the solved problem (documented in Solution.hpp)
 */
Solution solve(Sense sense, LinOp* objective, std::vector< LinOp* > constraints,
               std::map<std::string, double> solver_options);
#endif