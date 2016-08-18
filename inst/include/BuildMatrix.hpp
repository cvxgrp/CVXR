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

#ifndef BUILDMATRIX_H
#define BUILDMATRIX_H

#include <vector>
#include "LinOp.hpp"
#include "Utils.hpp"
#include "ProblemData.hpp"

ProblemData build_matrix(std::vector< LinOp* > constraints, std::map<int, int> id_to_col);
ProblemData build_matrix(std::vector< LinOp* > constraints, std::map<int, int> id_to_col, std::vector<int> constr_offsets);

#endif