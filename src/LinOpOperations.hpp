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

#ifndef LINOPOPERATIONS_H
#define LINOPOPERATIONS_H

#include <vector>
#include <map>
#include "Utils.hpp"
#include "LinOp.hpp"

std::map<int, Matrix> get_variable_coeffs(LinOp &lin);
std::map<int, Matrix> get_const_coeffs(LinOp &lin);
std::vector<Matrix> get_func_coeffs(LinOp &lin);
int get_id_data(LinOp &lin);

#endif