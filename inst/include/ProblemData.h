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

#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <vector>
#include <map>
#include <cassert>

/* Stores the result of calling BUILD_MATRIX on a collection of LinOp
 * trees. Uses tensor storage for DPP support.
 *
 * For non-DPP (no parameters): TensorV/I/J have a single outer key
 * (CONSTANT_ID = -1) with a single inner vector (index 0).
 * This degenerates to the same behavior as flat V/I/J.
 *
 * For DPP: TensorV/I/J have multiple outer keys (one per parameter +
 * CONSTANT_ID), each with a vector of inner vectors (one per parameter
 * element).
 */
class ProblemData {
public:
	/* COO sparse matrix representation stored per (param_id, element_idx).
	 * TensorV[param_id][element_idx] = vector of doubles (values)
	 * TensorI[param_id][element_idx] = vector of ints (row indices)
	 * TensorJ[param_id][element_idx] = vector of ints (col indices)
	 */
	std::map<int, std::vector<std::vector<double> > > TensorV;
	std::map<int, std::vector<std::vector<int> > > TensorI;
	std::map<int, std::vector<std::vector<int> > > TensorJ;

	// Pointers needed to extract V, I, J for a specific (param_id, vec_idx).
	int param_id;
	int vec_idx;

	// Initialize TensorV/I/J for the given parameter.
	void init_id(int new_param_id, int param_size) {
		assert(TensorV.count(new_param_id) == 0);
		std::vector<std::vector<double> > vecV(param_size);
		std::vector<std::vector<int> > vecI(param_size);
		std::vector<std::vector<int> > vecJ(param_size);
		TensorV[new_param_id] = vecV;
		TensorI[new_param_id] = vecI;
		TensorJ[new_param_id] = vecJ;
	}

	// Get length of V for current (param_id, vec_idx).
	int getLen() {
		return TensorV[param_id][vec_idx].size();
	}

	// Get the V data for current (param_id, vec_idx).
	std::vector<double>& getV() {
		return TensorV[param_id][vec_idx];
	}

	// Get the I data for current (param_id, vec_idx).
	std::vector<int>& getI() {
		return TensorI[param_id][vec_idx];
	}

	// Get the J data for current (param_id, vec_idx).
	std::vector<int>& getJ() {
		return TensorJ[param_id][vec_idx];
	}

	// Get all param_ids present in the tensor.
	std::vector<int> get_param_ids() {
		std::vector<int> ids;
		for (auto it = TensorV.begin(); it != TensorV.end(); ++it) {
			ids.push_back(it->first);
		}
		return ids;
	}

	// Get the number of inner vectors for a given param_id.
	int get_num_vecs(int pid) {
		if (TensorV.count(pid) == 0) return 0;
		return TensorV[pid].size();
	}

#ifdef _R_INTERFACE_
  ~ProblemData() {
#ifdef _R_DEBUG_
    Rcpp::Rcout << "ProblemData was destroyed!" <<std::endl;
#endif
  }
#endif
};

#endif
