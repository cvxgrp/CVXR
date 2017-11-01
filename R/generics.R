#'
#' Sign Properties
#' 
#' Determine if an expression is positive, negative, or zero.
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A logical value.
#' @name sign-methods
NULL

#' @rdname sign-methods
#' @export
setGeneric("is_zero", function(object) { standardGeneric("is_zero") })

#' @rdname sign-methods
#' @export
setGeneric("is_positive", function(object) { standardGeneric("is_positive") })

#' @rdname sign-methods
#' @export
setGeneric("is_negative", function(object) { standardGeneric("is_negative") })

#'
#' Curvature of Expression
#' 
#' The curvature of an expression.
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A string indicating the curvature of the expression, either "CONSTANT", "AFFINE", "CONVEX, "CONCAVE", or "UNKNOWN".
#' @docType methods
#' @rdname curvature
#' @export
setGeneric("curvature", function(object) { standardGeneric("curvature") })

#'
#' Curvature Properties
#' 
#' Determine if an expression is constant, affine, convex, concave, quadratic, or piecewise linear (pwl).
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A logical value.
#' @name curvature-methods
NULL

#' @rdname curvature-methods
#' @export
setGeneric("is_constant", function(object) { standardGeneric("is_constant") })

#' @rdname curvature-methods
#' @export
setGeneric("is_affine", function(object) { standardGeneric("is_affine") })

#' @rdname curvature-methods
#' @export
setGeneric("is_convex", function(object) { standardGeneric("is_convex") })

#' @rdname curvature-methods
#' @export
setGeneric("is_concave", function(object) { standardGeneric("is_concave") })

#' @rdname curvature-methods
#' @export
setGeneric("is_quadratic", function(object) { standardGeneric("is_quadratic") })

#' @rdname curvature-methods
#' @export
setGeneric("is_pwl", function(object) { standardGeneric("is_pwl") })

#'
#' Is an Expression DCP?
#' 
#' Determine if an expression complies with the disciplined convex programming rules.
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A logical value indicating whether the expression is DCP compliant, i.e. no unknown curvatures.
#' @docType methods
#' @rdname is_dcp
#' @export
setGeneric("is_dcp", function(object) { standardGeneric("is_dcp") })

#'
#' Size of Expression
#'
#' The size of an expression.
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A vector with two elements \code{c(row, col)} representing the dimensions of the expression.
#' @docType methods
#' @rdname size
#' @export
setGeneric("size", function(object) { standardGeneric("size") })

#'
#' Size Properties
#'
#' Determine if an expression is a scalar, vector, or matrix.
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A logical value.
#' @name size-methods
NULL

#' @rdname size-methods
#' @export
setGeneric("is_scalar", function(object) { standardGeneric("is_scalar") })

#' @rdname size-methods
#' @export
setGeneric("is_vector", function(object) { standardGeneric("is_vector") })

#' @rdname size-methods
#' @export
setGeneric("is_matrix", function(object) { standardGeneric("is_matrix") })

# The value of the objective given the solver primal value.
setGeneric("primal_to_result", function(object, result) { standardGeneric("primal_to_result") })

#'
#' Get or Set Value
#'
#' Get or set the value of a variable, parameter, expression, or problem.
#' 
#' @param object A \linkS4class{Variable}, \linkS4class{Parameter}, \linkS4class{Expression}, or \linkS4class{Problem} object.
#' @param value A numeric scalar, vector, or matrix to assign to the object.
#' @return The numeric value of the variable, parameter, or expression. If any part of the mathematical object is unknown, return \code{NA}.
#' @name value-methods
NULL

#' @rdname value-methods
#' @export
setGeneric("value", function(object) { standardGeneric("value") })

#' @rdname value-methods
#' @export
setGeneric("value<-", function(object, value) { standardGeneric("value<-") })

# Internal method for saving a value.
setGeneric("save_value", function(object, value) { standardGeneric("save_value") })

# Information needed to reconstruct the expression aside from its arguments.
setGeneric("get_data", function(object) { standardGeneric("get_data") })

#'
#' Variable, Parameter, or Expression Name
#'
#' The string representation of a variable, parameter, or expression.
#' 
#' @param object A \linkS4class{Variable}, \linkS4class{Parameter}, or \linkS4class{Expression} object.
#' @return For \linkS4class{Variable} or \linkS4class{Parameter} objects, the value in the name slot. For \linkS4class{Expression} objects, a string indicating the nested atoms and their respective arguments.
#' @docType methods
#' @rdname name
#' @export
setGeneric("name", function(object) { standardGeneric("name") })

#'
#' Parts of an Expression
#'
#' List the variables, parameters, or constants in a canonical expression.
#'
#' @param object A \linkS4class{Canonical} expression.
#' @return A list of \linkS4class{Variable}, \linkS4class{Parameter}, or \linkS4class{Constant} objects.
#' @name expression-parts
NULL

#' @rdname expression-parts
#' @export
setGeneric("variables", function(object) { standardGeneric("variables") })

#' @rdname expression-parts
#' @export
setGeneric("parameters", function(object) { standardGeneric("parameters") })

#' @rdname expression-parts
#' @export
setGeneric("constants", function(object) { standardGeneric("constants") })

#'
#' Sub/Super-Gradient
#'
#' The (sub/super)-gradient of the expression with respect to each variable.
#' Matrix expressions are vectorized, so the gradient is a matrix. \code{NA} indicates variable values are unknown or outside the domain.
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A list mapping each variable to a sparse matrix.
#' @docType methods
#' @rdname grad
#' @export
setGeneric("grad", function(object) { standardGeneric("grad") })

#'
#' Domain
#'
#' A list of constraints describing the closure of the region where the expression is finite.
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A list of \linkS4class{Constraint} objects.
#' @docType methods
#' @rdname domain
#' @export
setGeneric("domain", function(object) { standardGeneric("domain") })

#
# Validate Value
#
# Check that the value satisfies a \linkS4class{Leaf}'s symbolic attributes.
# 
# @param object A \linkS4class{Leaf} object.
# @param val The assigned value.
# @return The value converted to proper matrix type.
setGeneric("validate_val", function(object, val) { standardGeneric("validate_val") })

#'
#' Canonicalize
#'
#' Computes the graph implementation of a canonical expression.
#' 
#' @param object A \linkS4class{Canonical} object.
#' @return A list of \code{list(affine expression, list(constraints))}.
#' @docType methods
#' @rdname canonicalize
#' @aliases canonical_form
#' @export
setGeneric("canonicalize", function(object) { standardGeneric("canonicalize") })

# The canonical form of the input.
setGeneric("canonical_form", function(object) { standardGeneric("canonical_form") })

#
# Gradient of an Atom
#
# The (sub/super) gradient of the atom with respect to each argument. Matrix expressions are vectorized, so the gradient is a matrix.
# 
# @param object An \linkS4class{Atom} object.
# @param values A list of numeric values for the arguments.
# @return A list of sparse matrices or \code{NA}.
setGeneric(".grad", function(object, values) { standardGeneric(".grad") })

#
# Domain of an Atom
#
# The constraints describing the domain of the atom.
# 
# @param object An \linkS4class{Atom} object.
# @return A list of \linkS4class{Constraint} objects.
setGeneric(".domain", function(object) { standardGeneric(".domain") })

#
# Gradient of an AxisAtom
#
# The (sub/super) gradient of the atom with respect to each argument. Matrix expressions are vectorized, so the gradient is a matrix. Takes the axis into account.
# 
# @param values A list of numeric values for the arguments.
# @return A list of sparse matrices or \code{NA}.
setGeneric(".axis_grad", function(object, values) { standardGeneric(".axis_grad") })

#
# Column Gradient of an Atom
#
# The (sub/super) gradient of the atom with respect to a column argument. Matrix expressions are vectorized, so the gradient is a matrix.
# @param value A numeric value for a column.
# @return A sparse matrix or \code{NA}.
setGeneric(".column_grad", function(object, value) { standardGeneric(".column_grad") })

# Positive definite inequalities
#' @rdname PSDConstraint
#' @export
setGeneric("%>>%", function(e1, e2) { standardGeneric("%>>%") })

#' @rdname PSDConstraint
#' @export
setGeneric("%<<%", function(e1, e2) { standardGeneric("%<<%") })

#'
#' Atom Size
#'
#' Determine the size of an atom based on its arguments.
#' 
#' @param object A \linkS4class{Atom} object.
#' @return A numeric vector \code{c(row, col)} indicating the size of the atom.
#' @rdname size_from_args
setGeneric("size_from_args", function(object) { standardGeneric("size_from_args") })

#'
#' Atom Sign
#'
#' Determine the sign of an atom based on its arguments.
#' 
#' @param object A \linkS4class{Atom} object.
#' @return A logical vector \code{c(is positive, is negative)} indicating the sign of the atom.
#' @rdname sign_from_args
setGeneric("sign_from_args", function(object) { standardGeneric("sign_from_args") })
setGeneric("validate_args", function(object) { standardGeneric("validate_args") })
setGeneric("get_data", function(object) { standardGeneric("get_data") })
setGeneric("to_numeric", function(object, values) { standardGeneric("to_numeric") })

#'
#' Curvature of an Atom
#'
#' Determine if an atom is convex, concave, or affine.
#' 
#' @param object A \linkS4class{Atom} object.
#' @return A logical value.
#' @name curvature-atom
NULL

#' @rdname curvature-atom
#' @export
setGeneric("is_atom_convex", function(object) { standardGeneric("is_atom_convex") })

#' @rdname curvature-atom
#' @export
setGeneric("is_atom_concave", function(object) { standardGeneric("is_atom_concave") })

#' @rdname curvature-atom
#' @export
setGeneric("is_atom_affine", function(object) { standardGeneric("is_atom_affine") })

#'
#' Curvature of Composition
#'
#' Determine whether a composition is increasing or decreasing in \code{idx}.
#' 
#' @param object A \linkS4class{Atom} object.
#' @param idx An index into the atom.
#' @return A logical value.
#' @name curvature-comp
NULL

#' @rdname curvature-comp
#' @export
setGeneric("is_incr", function(object, idx) { standardGeneric("is_incr") })

#' @rdname curvature-comp
#' @export
setGeneric("is_decr", function(object, idx) { standardGeneric("is_decr") })

#'
#' Graph Implementation
#' 
#' Reduces the atom to an affine expression and list of constraints.
#' 
#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @return A list of \code{list(LinOp for objective, list of constraints)}, where LinOp is a list representing the linear operator.
#' @docType methods
#' @rdname graph_implementation
setGeneric("graph_implementation", function(object, arg_objs, size, data) { standardGeneric("graph_implementation") })

# Constraint generic methods
setGeneric("id", function(object) { standardGeneric("id") })
setGeneric("residual", function(object) { standardGeneric("residual") })
setGeneric("violation", function(object) { standardGeneric("violation") })
setGeneric("num_cones", function(object) { standardGeneric("num_cones") })
setGeneric("cone_size", function(object) { standardGeneric("cone_size") })
setGeneric("dual_value", function(object) { standardGeneric("dual_value") })

# 
# Format Constraints
# 
# Formats constraints for the solver.
# 
# @param object A \linkS4class{Constraint} object.
# @param eq_constr A list of the equality constraints in the canonical problem.
# @param leq_constr A list of the inequality constraints in the canonical problem.
# @param dims A list with the dimensions of the conic constraints.
# @param solver A string representing the solver to be called.
# @return A list containing equality constraints, inequality constraints, and dimensions.
# @rdname format_constr
setGeneric("format_constr", function(object, eq_constr, leq_constr, dims, solver) { standardGeneric("format_constr") })
setGeneric("constr_type", function(object) { standardGeneric("constr_type") })
setGeneric("constr_id", function(object) { standardGeneric("constr_id") })

setGeneric("block_add", function(object, mat, block, vert_offset, horiz_offset, rows, cols, vert_step, horiz_step) { standardGeneric("block_add") })
setGeneric("place_x0", function(object, big_x, var_offsets) { standardGeneric("place_x0") })
setGeneric("place_Df", function(object, big_Df, Df, var_offsets, vert_offset) { standardGeneric("place_Df") })
setGeneric("place_H", function(object, big_H, H, var_offsets) { standardGeneric("place_H") })
setGeneric("extract_variables", function(object, x, var_offsets) { standardGeneric("extract_variables") })

# Problem generic methods
#'
#' Parts of a Problem
#'
#' Get the objective, constraints, status (from the last call to \code{solve}), size metrics, or solver statistics of a problem.
#'
#' @param object A \linkS4class{Problem} object.
#' @return The requested slot of the object.
#' @name problem-parts
NULL

#' @rdname problem-parts
#' @export
setGeneric("objective", function(object) { standardGeneric("objective") })

#' @rdname problem-parts
#' @export
setGeneric("constraints", function(object) { standardGeneric("constraints") })

#' @rdname problem-parts
#' @export
setGeneric("status", function(object) { standardGeneric("status") })
setGeneric("status<-", function(object, value) { standardGeneric("status<-") })

#' @rdname problem-parts
#' @export
setGeneric("size_metrics", function(object) { standardGeneric("size_metrics") })

#' @rdname problem-parts
#' @export
setGeneric("solver_stats", function(object) { standardGeneric("solver_stats") })
setGeneric("solver_stats<-", function(object, value) { standardGeneric("solver_stats<-") })

#'
#' Get Problem Data
#'
#' Get the problem data used in the call to the solver.
#' 
#' @param object A \linkS4class{Problem} object.
#' @param solver A string indicating the solver that the problem data is for.
#' @return A list of arguments for the solver.
#' @rdname get_problem_data
#' @export
setGeneric("get_problem_data", function(object, solver) { standardGeneric("get_problem_data") })

setGeneric("is_qp", function(object) { standardGeneric("is_qp") })
setGeneric("unpack_results", function(object, solver, results_dict) { standardGeneric("unpack_results") })
setGeneric(".handle_no_solution", function(object, status) { standardGeneric(".handle_no_solution") })
setGeneric("Problem.save_values", function(object, result_vec, objstore, offset_map) { standardGeneric("Problem.save_values") })
setGeneric(".save_dual_values", function(object, result_vec, constraints, constr_types) { standardGeneric(".save_dual_values") })
setGeneric(".update_problem_state", function(object, results_dict, sym_data, solver) { standardGeneric(".update_problem_state") })

# Problem data generic methods
setGeneric("get_objective", function(object) { standardGeneric("get_objective") })
setGeneric("get_eq_constr", function(object) { standardGeneric("get_eq_constr") })
setGeneric("get_ineq_constr", function(object) { standardGeneric("get_ineq_constr") })
setGeneric("get_nonlin_constr", function(object) { standardGeneric("get_nonlin_constr") })

# Solver generic methods
setGeneric("choose_solution", function(solver, results_dict) { standardGeneric("choose_solution") })
setGeneric("import_solver", function(solver) { standardGeneric("import_solver") })
setGeneric("nonlin_constr", function(solver) { standardGeneric("nonlin_constr") })
setGeneric("validate_solver", function(solver, constraints) { standardGeneric("validate_solver") })
setGeneric("validate_cache", function(solver, objective, constraints, cached_data) { standardGeneric("validate_cache") })
setGeneric("get_sym_data", function(solver, objective, constraints, cached_data) { standardGeneric("get_sym_data") })
setGeneric("get_matrix_data", function(solver, objective, constraints, cached_data) { standardGeneric("get_matrix_data") })
setGeneric("Solver.get_problem_data", function(solver, objective, constraints, cached_data) { standardGeneric("Solver.get_problem_data") })

setGeneric("matrix_intf", function(solver) { standardGeneric("matrix_intf") })
setGeneric("vec_intf", function(solver) { standardGeneric("vec_intf") })
setGeneric("split_constr", function(solver, constr_map) { standardGeneric("split_constr") })
setGeneric("Solver.solve", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) { standardGeneric("Solver.solve") })
setGeneric("format_results", function(solver, results_dict, data, cached_data) { standardGeneric("format_results") })

#'
#' Solver Capabilities
#' 
#' Determine if a solver is capable of solving a linear program (LP), second-order cone program (SOCP), semidefinite program (SDP), exponential cone program (EXP), or mixed-integer program (MIP).
#' 
#' @param solver A \linkS4class{Solver} object.
#' @return A logical value.
#' @name Solver-capable
NULL

#' @rdname Solver-capable
#' @export
setGeneric("lp_capable", function(solver) { standardGeneric("lp_capable") })

#' @rdname Solver-capable
#' @export
setGeneric("socp_capable", function(solver) { standardGeneric("socp_capable") })

#' @rdname Solver-capable
#' @export
setGeneric("sdp_capable", function(solver) { standardGeneric("sdp_capable") })

#' @rdname Solver-capable
#' @export
setGeneric("exp_capable", function(solver) { standardGeneric("exp_capable") })

#' @rdname Solver-capable
#' @export
setGeneric("mip_capable", function(solver) { standardGeneric("mip_capable") })

# Map of solver status code to CVXR status.
setGeneric("status_map", function(solver, status) { standardGeneric("status_map") })
