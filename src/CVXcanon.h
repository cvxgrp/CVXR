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

#ifndef CVXCANON_H
#define CVXCANON_H

#include <vector>
#include "LinOp.h"
#include "Utils.h"
#include "ProblemData.h"

#ifdef _R_INTERFACE_

// DPP tensor-based entry points.
// var_length = total number of variable columns (constant goes at column var_length).
// param_to_size: map from param_id to param_size (empty for non-DPP).
void build_matrix_2(std::vector< LinOp* > constraints,
		    int var_length,
		    std::map<int, int> id_to_col,
		    std::map<int, int> param_to_size,
		    Rcpp::XPtr<ProblemData> prob_data);

void build_matrix_3(std::vector<LinOp*> constraints,
		    int var_length,
		    std::map<int, int> id_to_col,
		    std::map<int, int> param_to_size,
		    std::vector<int> constr_offsets,
		    Rcpp::XPtr<ProblemData> prob_data);
#else

ProblemData build_matrix(std::vector< LinOp* > constraints, int var_length, std::map<int, int> id_to_col, std::map<int, int> param_to_size);
ProblemData build_matrix(std::vector< LinOp* > constraints, int var_length, std::map<int, int> id_to_col, std::map<int, int> param_to_size, std::vector<int> constr_offsets);

#endif

#endif
