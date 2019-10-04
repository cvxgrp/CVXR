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

#include "CVXcanon.h"
#include <iostream>
#include <map>
#include "LinOp.h"
#include "LinOpOperations.h"
#include "ProblemData.h"

void mul_by_const(Matrix &coeff_mat,
        std::map<int, Matrix > &rh_coeffs,
        std::map<int, Matrix > &result){

	typedef std::map<int, Matrix >::iterator it_type;
	for (it_type it = rh_coeffs.begin(); it != rh_coeffs.end(); ++it){
		int id = it->first;
		Matrix rh = it->second;
		/* Convert scalars (1x1 matrices) to primitive types */
		if (coeff_mat.rows() == 1 && coeff_mat.cols() == 1){
			double scalar = coeff_mat.coeffRef(0, 0);
			if(result.count(id) == 0)
				result[id] = scalar * rh;
			else
				result[id] += scalar * rh;
		} else if (rh.rows() == 1 && rh.cols() == 1) {
			double scalar = rh.coeffRef(0, 0);
			if(result.count(id) == 0)
				result[id] = coeff_mat * scalar;
			else
				result[id] = coeff_mat * scalar;
		} else{
			if(result.count(id) == 0 )
				result[id] = coeff_mat * rh;
			else
				result[id] += coeff_mat * rh;
		}
	}
}

std::map<int, Matrix > get_coefficient(LinOp &lin){
	std::map<int, Matrix > coeffs;
	if (lin.type == VARIABLE){
		std::map<int, Matrix> new_coeffs = get_variable_coeffs(lin);
		typedef std::map<int, Matrix >::iterator it_type;
		for(it_type it = new_coeffs.begin(); it != new_coeffs.end(); ++it){
			if(coeffs.count(it->first) == 0)
				coeffs[it->first] = it->second ;
			else
				coeffs[it->first] += it->second;
		}
	}
	else if (lin.has_constant_type()){
		/* ID will be CONSTANT_TYPE */
		std::map<int, Matrix> new_coeffs = get_const_coeffs(lin);
		typedef std::map<int, Matrix >::iterator it_type;
		for(it_type it = new_coeffs.begin(); it != new_coeffs.end(); ++it){
			if(coeffs.count(it->first) == 0)
				coeffs[it->first] = it->second;
			else
				coeffs[it->first] += it->second;
		}
	}
	else {
		/* Multiply the arguments of the function coefficient in order */
		std::vector<Matrix> coeff_mat = get_func_coeffs(lin); 
		for (unsigned i = 0; i < lin.args.size(); i++){
			Matrix coeff = coeff_mat[i];
			std::map<int, Matrix > rh_coeffs = get_coefficient(*lin.args[i]);
			std::map<int,  Matrix > new_coeffs;
			mul_by_const(coeff, rh_coeffs, new_coeffs);

			typedef std::map<int, Matrix>::iterator it_type;
			for (it_type it = new_coeffs.begin(); it != new_coeffs.end(); ++it){
				if(coeffs.count(it->first) == 0)
					coeffs[it->first] = it->second;
				else
					coeffs[it->first] += it->second;		
			}
		}
	}
	return coeffs;
}

int get_horiz_offset(int id, std::map<int, int> &offsets,
                     int &horiz_offset, LinOp &lin){
	if ( !offsets.count(id) ){
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
                           int &vert_offset, int &horiz_offset){
	for ( int k = 0; k < block.outerSize(); ++k ) {
		for ( Matrix::InnerIterator it(block, k); it; ++it ){
			V.push_back(it.value());

			/* Push back current row and column indices */
			I.push_back(it.row() + vert_offset);   	
			J.push_back(it.col() + horiz_offset);
		}
	}
}

void extend_constant_vec(std::vector<double> &const_vec, int &vert_offset,
                         Matrix &block){
	int rows = block.rows();
	for ( int k = 0; k < block.outerSize(); ++k ){
		for ( Matrix::InnerIterator it(block, k); it; ++it ){
			int idx = vert_offset + (it.col() * rows) + it.row();
			const_vec[idx] += it.value();
		}
	}
}

void process_constraint(LinOp & lin, std::vector<double> &V,
                        std::vector<int> &I, std::vector<int> &J,
                        std::vector<double> &constant_vec, int &vert_offset,
                        std::map<int, int> &id_to_col, int & horiz_offset){
	/* Get the coefficient for the current constraint */
	std::map<int, Matrix > coeffs = get_coefficient(lin);	

	typedef std::map<int, Matrix >::iterator it_type;
	for(it_type it = coeffs.begin(); it != coeffs.end(); ++it){
		int id = it->first;	// Horiz offset determined by the id
		Matrix block = it->second;
		if (id == CONSTANT_ID) { // Add to CONSTANT_VEC if linop is constant
			extend_constant_vec(constant_vec, vert_offset, block);	
		}
		else {
			int offset = get_horiz_offset(id, id_to_col, horiz_offset, lin);
			add_matrix_to_vectors(block, V, I, J, vert_offset, offset);
		}
	}
}

/* Returns the number of rows in the matrix assuming vertical stacking
	 of coefficient matrices */
int get_total_constraint_length(std::vector< LinOp* > constraints){
	int result = 0;
#ifdef _R_DEBUG_
	Rcpp::Rcout << "In get_total_constraint_length, size = " << constraints.size() << std::endl;
#endif
	for (unsigned i = 0; i < constraints.size(); i++){
#ifdef _R_DEBUG_
	  Rcpp::Rcout << "i, r, c, " << constraints[i]->size[0] << ", " << constraints[i]->size[1] << std::endl;
#endif
		result += constraints[i]->size[0] * constraints[i]->size[1];
	}
	return result;
}

/* Returns the number of rows in the matrix using the user provided vertical
	 offsets for each constraint. */
int get_total_constraint_length(std::vector<LinOp*> &constraints, 
									 							std::vector<int> &constr_offsets){
	/* Must specify an offset for each constraint */
	if(constraints.size() != constr_offsets.size()){
#ifdef _R_INTERFACE_
	        Rcpp::stop("Invalid constraint offsets: CONSTR_OFFSET must be the same length as CONSTRAINTS");
#else
  	        std::cerr << "Error: Invalid constraint offsets: ";
		std::cerr	<< "CONSTR_OFFSET must be the same length as CONSTRAINTS" << std::endl;
		exit(-1);
#endif
	}

	int offset_end = 0;
	/* Offsets must be monotonically increasing */
#ifdef _R_INTERFACE_
	for(unsigned i = 0; i < constr_offsets.size(); i++){
		LinOp *constr = constraints[i];
		int offset_start = constr_offsets[i];
		offset_end = offset_start + constr->size[0] * constr->size[1];

		if(i + 1 < constr_offsets.size() && constr_offsets[i + 1] < offset_end){
#ifdef _R_INTERFACE_
	                Rcpp::stop("Invalid constraint offsets: offsets are not monotonically increasing");
#else
			std::cerr << "Error: Invalid constraint offsets: ";
			std::cerr << "Offsets are not monotonically increasing" << std::endl;
			exit(-1);
#endif
		}
	}
#else
	for(unsigned i = 0; i < constr_offsets.size(); i++){
		LinOp constr = *constraints[i];
		int offset_start = constr_offsets[i];
		offset_end = offset_start + constr.size[0] * constr.size[1];

		if(i + 1 < constr_offsets.size() && constr_offsets[i + 1] < offset_end){
#ifdef _R_INTERFACE_
	                Rcpp::stop("Invalid constraint offsets: offsets are not monotonically increasing");
#else
			std::cerr << "Error: Invalid constraint offsets: ";
			std::cerr << "Offsets are not monotonically increasing" << std::endl;
			exit(-1);
#endif
		}
	}
#endif
	return offset_end;
}


#ifdef _R_INTERFACE_

/* function: build_matrix
*
* Description: Given a list of linear operations, this function returns a data
* structure containing a sparse matrix representation of the cone program.
*
* Input: std::vector<LinOp *> constraints, our list of constraints represented
* as a linear operation tree and an Rcpp XPtr to a newly created ProblemData instance,
* a data structure which contains a sparse representation
* of the coefficient matrix, a dense representation of the constant vector,
* and maps containing our mapping from variables, and a map from the rows of our
* matrix to their corresponding constraint.
* Output: None
*
*/
void build_matrix_2(std::vector< LinOp* > constraints,
		    std::map<int, int> id_to_col,
		    Rcpp::XPtr<ProblemData> prob_data) {
#ifdef _R_DEBUG_
  Rcpp::Rcout << "In Build_matrix_2" <<std::endl;    
#endif  

  int num_rows = get_total_constraint_length(constraints);
#ifdef _R_DEBUG_
  Rcpp::Rcout << "In Build_matrix_2 after length constraints " << num_rows <<std::endl;    
#endif  

  prob_data->const_vec = std::vector<double> (num_rows, 0);
#ifdef _R_DEBUG_
  Rcpp::Rcout << "Build_matrix_2 before id_to_col" <<std::endl;    
#endif  

  prob_data->id_to_col = id_to_col;
#ifdef _R_DEBUG_
  Rcpp::Rcout << "Build_matrix_2 after id_to_col" <<std::endl;    
#endif  

  int vert_offset = 0;
  int horiz_offset  = 0;

  /* Build matrix one constraint at a time */
  for (unsigned i = 0; i < constraints.size(); i++){
#ifdef _R_DEBUG_
    Rcpp::Rcout << "In Build_matrix_2 loop " << i << std::endl;    
#endif  
    LinOp *constr = constraints[i];
    process_constraint(*constr, prob_data->V, prob_data->I, prob_data->J,
		       prob_data->const_vec, vert_offset,
		       prob_data->id_to_col, horiz_offset);
    prob_data->const_to_row[i] = vert_offset;
    vert_offset += constr->size[0] * constr->size[1];
  }

}

/*  See comment above for build_matrix. Requires specification of a vertical
    offset, VERT_OFFSET, for each constraint in the vector CONSTR_OFFSETS. 
    
    Valid CONSTR_OFFSETS assume that a vertical offset is provided for each 
    constraint and that the offsets are not overlapping. In particular, 
    the vertical offset for constraint i + the size of constraint i must be
    less than the vertical offset for constraint i+1.
*/
void build_matrix_3(std::vector<LinOp*> constraints,
		    std::map<int, int> id_to_col,
		    std::vector<int> constr_offsets,
		    Rcpp::XPtr<ProblemData> prob_data){
  
  /* Function also verifies the offsets are valid */
  int num_rows = get_total_constraint_length(constraints, constr_offsets);
  prob_data->const_vec = std::vector<double> (num_rows, 0);
  prob_data->id_to_col = id_to_col;
  int horiz_offset  = 0;
  
  /* Build matrix one constraint at a time */
  for (unsigned i = 0; i < constraints.size(); i++){
    LinOp *constr = constraints[i];
#ifdef _R_DEBUG_
    Rcpp::Rcout << "Processing constraint " << i << std::endl;
#endif
    int vert_offset = constr_offsets[i];
    process_constraint(*constr, prob_data->V, prob_data->I, prob_data->J,
		       prob_data->const_vec, vert_offset,
		       prob_data->id_to_col, horiz_offset);
    prob_data->const_to_row[i] = vert_offset;
  }
}

#else

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
	for (unsigned i = 0; i < constraints.size(); i++){
		LinOp constr = *constraints[i];
		process_constraint(constr, prob_data.V, prob_data.I, prob_data.J,
		                   prob_data.const_vec, vert_offset,
		                   prob_data.id_to_col, horiz_offset);
		prob_data.const_to_row[i] = vert_offset;
		vert_offset += constr.size[0] * constr.size[1];
	}
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
                         std::vector<int> constr_offsets){
	ProblemData prob_data;

	/* Function also verifies the offsets are valid */
	int num_rows = get_total_constraint_length(constraints, constr_offsets);
	prob_data.const_vec = std::vector<double> (num_rows, 0);
	prob_data.id_to_col = id_to_col;
	int horiz_offset  = 0;

	/* Build matrix one constraint at a time */
	for (unsigned i = 0; i < constraints.size(); i++){
		LinOp constr = *constraints[i];
		int vert_offset = constr_offsets[i];
		process_constraint(constr, prob_data.V, prob_data.I, prob_data.J,
		                   prob_data.const_vec, vert_offset,
		                   prob_data.id_to_col, horiz_offset);
		prob_data.const_to_row[i] = vert_offset;
	}
	return prob_data;
}

#endif
