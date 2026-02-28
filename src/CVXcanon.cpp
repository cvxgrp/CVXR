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

// DPP tensor-based canonicalization.
// Mirrors CVXPY cvxcore/src/cvxcore.cpp.

#include "CVXcanon.h"
#include <iostream>
#include <map>
#include <utility>
#include "LinOp.h"
#include "LinOpOperations.h"
#include "ProblemData.h"

#ifdef _OPENMP
#include "omp.h"
#endif


/* function: add_matrix_to_vectors
*
* Adds a matrix to our sparse matrix triplet representation.
* Takes horizontal and vertical offsets for positioning within
* the larger matrix.
*/
void add_matrix_to_vectors(const Matrix &block, std::vector<double> &V,
                           std::vector<int>  &I, std::vector<int> &J,
                           int vert_offset, int horiz_offset){
	auto number_of_nonzeros = block.nonZeros();
	V.reserve(V.size() + number_of_nonzeros);
	I.reserve(I.size() + number_of_nonzeros);
	J.reserve(J.size() + number_of_nonzeros);

	for ( int k = 0; k < block.outerSize(); ++k ) {
		for ( Matrix::InnerIterator it(block, k); it; ++it ){
			V.push_back(it.value());
			I.push_back(it.row() + vert_offset);
			J.push_back(it.col() + horiz_offset);
		}
	}
}

/* function: process_constraint
*
* DPP tensor-based constraint processing.
* Mirrors CVXPY cvxcore.cpp process_constraint.
*
* Gets tensor coefficients via lin_to_tensor, then distributes
* the sparse triplets into the ProblemData's TensorV/I/J arrays.
* Constants (var_id == CONSTANT_ID) go at column var_length.
*/
void process_constraint(LinOp &lin, ProblemData &problemData,
                        int vert_offset, int var_length,
                        const std::map<int, int> &id_to_col) {
	Tensor coeffs = lin_to_tensor(lin);

	for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
		int param_id = it->first;
		/* Safety: skip param_ids not initialized in ProblemData.
		 * This happens when coeff_quad_form calls get_problem_matrix
		 * with empty param_to_size but the LinOp tree has PARAM nodes
		 * (because linOp_data is set for all data nodes). */
		if (problemData.TensorV.find(param_id) == problemData.TensorV.end()) {
			continue;
		}
		const DictMat& var_map = it->second;
		for (auto in_it = var_map.begin(); in_it != var_map.end(); ++in_it) {
			int var_id = in_it->first;
			const std::vector<Matrix>& blocks = in_it->second;
			for (unsigned i = 0; i < blocks.size(); ++i) {
				int horiz_offset;
				if (var_id == CONSTANT_ID) {
					horiz_offset = var_length;
				} else {
					horiz_offset = id_to_col.at(var_id);
				}

				/* Safety: skip if inner index exceeds initialized size */
				if (i >= problemData.TensorV[param_id].size()) {
					continue;
				}

				#ifdef _OPENMP
				#pragma omp critical
				#endif
				{
					add_matrix_to_vectors(blocks[i],
					    problemData.TensorV[param_id][i],
					    problemData.TensorI[param_id][i],
					    problemData.TensorJ[param_id][i],
					    vert_offset, horiz_offset);
				}
			}
		}
	}
}

/* Returns the number of rows assuming vertical stacking of constraints */
int get_total_constraint_length(std::vector< LinOp* > constraints){
	int result = 0;
	for (unsigned i = 0; i < constraints.size(); i++){
		result += constraints[i]->size[0] * constraints[i]->size[1];
	}
	return result;
}

/* Returns the number of rows using user-provided vertical offsets */
int get_total_constraint_length(std::vector<LinOp*> &constraints,
                                std::vector<int> &constr_offsets){
	if(constraints.size() != constr_offsets.size()){
#ifdef _R_INTERFACE_
		Rcpp::stop("Invalid constraint offsets: CONSTR_OFFSET must be the same length as CONSTRAINTS");
#else
		std::cerr << "Error: Invalid constraint offsets" << std::endl;
		exit(-1);
#endif
	}

	int offset_end = 0;
	for(unsigned i = 0; i < constr_offsets.size(); i++){
		LinOp *constr = constraints[i];
		int offset_start = constr_offsets[i];
		offset_end = offset_start + constr->size[0] * constr->size[1];

		if(i + 1 < constr_offsets.size() && constr_offsets[i + 1] < offset_end){
#ifdef _R_INTERFACE_
			Rcpp::stop("Invalid constraint offsets: offsets are not monotonically increasing");
#else
			std::cerr << "Error: Invalid constraint offsets" << std::endl;
			exit(-1);
#endif
		}
	}
	return offset_end;
}

/* Create ProblemData with tensor entries for each parameter.
 * Always creates CONSTANT_ID with 1 inner vector.
 * For each param_id, creates param_size inner vectors. */
ProblemData init_data_tensor(std::map<int, int> param_to_size) {
	ProblemData output;
	output.init_id(CONSTANT_ID, 1);
	typedef std::map<int, int>::iterator it_type;
	for (it_type it = param_to_size.begin(); it != param_to_size.end(); ++it) {
		int param_id = it->first;
		int param_size = it->second;
		output.init_id(param_id, param_size);
	}
	return output;
}


#ifdef _R_INTERFACE_

/* function: build_matrix_2
*
* DPP tensor-based build_matrix. Mirrors CVXPY build_matrix.
* var_length: total number of variable columns.
* id_to_col: pre-computed variable-id to column-offset mapping.
* param_to_size: parameter-id to parameter-size mapping (empty for non-DPP).
*/
void build_matrix_2(std::vector< LinOp* > constraints,
		    int var_length,
		    std::map<int, int> id_to_col,
		    std::map<int, int> param_to_size,
		    Rcpp::XPtr<ProblemData> prob_data) {

	ProblemData init = init_data_tensor(param_to_size);
	prob_data->TensorV = init.TensorV;
	prob_data->TensorI = init.TensorI;
	prob_data->TensorJ = init.TensorJ;

	int vert_offset = 0;
	for (unsigned i = 0; i < constraints.size(); i++){
		LinOp *constr = constraints[i];
		process_constraint(*constr, *prob_data, vert_offset,
		                   var_length, id_to_col);
		vert_offset += constr->size[0] * constr->size[1];
	}
}

/* build_matrix_3: with user-provided constraint offsets */
void build_matrix_3(std::vector<LinOp*> constraints,
		    int var_length,
		    std::map<int, int> id_to_col,
		    std::map<int, int> param_to_size,
		    std::vector<int> constr_offsets,
		    Rcpp::XPtr<ProblemData> prob_data){

	// Verify offsets
	get_total_constraint_length(constraints, constr_offsets);

	ProblemData init = init_data_tensor(param_to_size);
	prob_data->TensorV = init.TensorV;
	prob_data->TensorI = init.TensorI;
	prob_data->TensorJ = init.TensorJ;

	for (unsigned i = 0; i < constraints.size(); i++){
		LinOp *constr = constraints[i];
		int vert_offset = constr_offsets[i];
		process_constraint(*constr, *prob_data, vert_offset,
		                   var_length, id_to_col);
	}
}

#else

ProblemData build_matrix(std::vector< LinOp* > constraints,
                         int var_length,
                         std::map<int, int> id_to_col,
                         std::map<int, int> param_to_size) {
	ProblemData prob_data = init_data_tensor(param_to_size);
	int vert_offset = 0;
	for (unsigned i = 0; i < constraints.size(); i++){
		LinOp constr = *constraints[i];
		process_constraint(constr, prob_data, vert_offset,
		                   var_length, id_to_col);
		vert_offset += constr.size[0] * constr.size[1];
	}
	return prob_data;
}

ProblemData build_matrix(std::vector<LinOp*> constraints,
                         int var_length,
                         std::map<int, int> id_to_col,
                         std::map<int, int> param_to_size,
                         std::vector<int> constr_offsets){
	ProblemData prob_data = init_data_tensor(param_to_size);
	get_total_constraint_length(constraints, constr_offsets);
	for (unsigned i = 0; i < constraints.size(); i++){
		LinOp constr = *constraints[i];
		int vert_offset = constr_offsets[i];
		process_constraint(constr, prob_data, vert_offset,
		                   var_length, id_to_col);
	}
	return prob_data;
}

#endif
