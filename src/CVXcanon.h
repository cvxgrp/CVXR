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

// Top Level Entry point
void build_matrix_2(std::vector< LinOp* > constraints,
		    std::map<int, int> id_to_col,
		    Rcpp::XPtr<ProblemData> prob_data);

void build_matrix_3(std::vector<LinOp*> constraints,
		    std::map<int, int> id_to_col,
		    std::vector<int> constr_offsets,
		    Rcpp::XPtr<ProblemData> prob_data);
#else

// Top Level Entry point
ProblemData build_matrix(std::vector< LinOp* > constraints, std::map<int, int> id_to_col);
ProblemData build_matrix(std::vector< LinOp* > constraints, std::map<int, int> id_to_col, std::vector<int> constr_offsets);

#endif

#endif
