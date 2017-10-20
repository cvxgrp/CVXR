# DCP attribute generic methods
setGeneric("is_zero", function(object) { standardGeneric("is_zero") })
setGeneric("is_positive", function(object) { standardGeneric("is_positive") })
setGeneric("is_negative", function(object) { standardGeneric("is_negative") })
setGeneric("is_unknown", function(object) { standardGeneric("is_unknown") })

setGeneric("is_constant", function(object) { standardGeneric("is_constant") })
setGeneric("is_affine", function(object) { standardGeneric("is_affine") })
setGeneric("is_convex", function(object) { standardGeneric("is_convex") })
setGeneric("is_concave", function(object) { standardGeneric("is_concave") })
setGeneric("is_quadratic", function(object) { standardGeneric("is_quadratic") })
setGeneric("is_pwl", function(object) { standardGeneric("is_pwl") })
setGeneric("is_dcp", function(object) { standardGeneric("is_dcp") })
setGeneric("size", function(object) { standardGeneric("size") })
setGeneric("primal_to_result", function(object, result) { standardGeneric("primal_to_result") })

# Expression generic methods
setGeneric("value", function(object) { standardGeneric("value") })
setGeneric("value<-", function(object, value) { standardGeneric("value<-") })
setGeneric("save_value", function(object, value) { standardGeneric("save_value") })
setGeneric("get_data", function(object) { standardGeneric("get_data") })
setGeneric("curvature", function(object) { standardGeneric("curvature") })
setGeneric("is_scalar", function(object) { standardGeneric("is_scalar") })
setGeneric("is_vector", function(object) { standardGeneric("is_vector") })
setGeneric("is_matrix", function(object) { standardGeneric("is_matrix") })

setGeneric("name", function(object) { standardGeneric("name") })
setGeneric("variables", function(object) { standardGeneric("variables") })
setGeneric("parameters", function(object) { standardGeneric("parameters") })
setGeneric("constants", function(object) { standardGeneric("constants") })
setGeneric("grad", function(object) { standardGeneric("grad") })
setGeneric("domain", function(object) { standardGeneric("domain") })
setGeneric("validate_val", function(object, val) { standardGeneric("validate_val") })
setGeneric("canonical_form", function(object) { standardGeneric("canonical_form") })
setGeneric("canonicalize", function(object) { standardGeneric("canonicalize") })

setGeneric(".grad", function(object, values) { standardGeneric(".grad") })
setGeneric(".domain", function(object) { standardGeneric(".domain") })
setGeneric(".axis_grad", function(object, values) { standardGeneric(".axis_grad") })
setGeneric(".column_grad", function(object, value) { standardGeneric(".column_grad") })

# Positive definite inequalities
setGeneric("%>>%", function(e1, e2) { standardGeneric("%>>%") })
setGeneric("%<<%", function(e1, e2) { standardGeneric("%<<%") })

# Atom generic methods
setGeneric("validate_args", function(object) { standardGeneric("validate_args") })
setGeneric("size_from_args", function(object) { standardGeneric("size_from_args") })
setGeneric("sign_from_args", function(object) { standardGeneric("sign_from_args") })
setGeneric("get_data", function(object) { standardGeneric("get_data") })
setGeneric("to_numeric", function(object, values) { standardGeneric("to_numeric") })

setGeneric("is_atom_convex", function(object) { standardGeneric("is_atom_convex") })
setGeneric("is_atom_concave", function(object) { standardGeneric("is_atom_concave") })
setGeneric("is_atom_affine", function(object) { standardGeneric("is_atom_affine") })
setGeneric("is_incr", function(object, idx) { standardGeneric("is_incr") })
setGeneric("is_decr", function(object, idx) { standardGeneric("is_decr") })
setGeneric("graph_implementation", function(object, arg_objs, size, data) { standardGeneric("graph_implementation") })

# Constraint generic methods
setGeneric("id", function(object) { standardGeneric("id") })
setGeneric("residual", function(object) { standardGeneric("residual") })
setGeneric("violation", function(object) { standardGeneric("violation") })
setGeneric("num_cones", function(object) { standardGeneric("num_cones") })
setGeneric("cone_size", function(object) { standardGeneric("cone_size") })
setGeneric("dual_value", function(object) { standardGeneric("dual_value") })

#' 
#' Format Constraints
#' 
#' Formats constraints for the solver.
#' 
#' @param object A \linkS4class{Constraint} object.
#' @param eq_constr A list of the equality constraints in the canonical problem.
#' @param leq_constr A list of the inequality constraints in the canonical problem.
#' @param dims A list with the dimensions of the conic constraints.
#' @param solver A string representing the solver to be called.
#' @return A list containing equality constraints, inequality constraints, and dimensions.
#' @docType methods
#' @rdname format_constr
setGeneric("format_constr", function(object, eq_constr, leq_constr, dims, solver) { standardGeneric("format_constr") })
setGeneric("constr_type", function(object) { standardGeneric("constr_type") })
setGeneric("constr_id", function(object) { standardGeneric("constr_id") })

setGeneric("block_add", function(object, mat, block, vert_offset, horiz_offset, rows, cols, vert_step, horiz_step) { standardGeneric("block_add") })
setGeneric("place_x0", function(object, big_x, var_offsets) { standardGeneric("place_x0") })
setGeneric("place_Df", function(object, big_Df, Df, var_offsets, vert_offset) { standardGeneric("place_Df") })
setGeneric("place_H", function(object, big_H, H, var_offsets) { standardGeneric("place_H") })
setGeneric("extract_variables", function(object, x, var_offsets) { standardGeneric("extract_variables") })

# Problem generic methods
setGeneric("status", function(object) { standardGeneric("status") })
setGeneric("status<-", function(object, value) { standardGeneric("status<-") })
setGeneric("size_metrics", function(object) { standardGeneric("size_metrics") })
setGeneric("solver_stats", function(object) { standardGeneric("solver_stats") })
setGeneric("solver_stats<-", function(object, value) { standardGeneric("solver_stats<-") })
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
#' The types of convex problems that a solver is capable of solving.
#' 
#' @param solver A \linkS4class{Solver} object.
#' @return A logical value indicating the solver capability.
#' @docType methods
#' @rdname Solver-capable
NULL

#' @describeIn Solver-capable A logical value indicating whether the solver is capable of solving linear programs.
setGeneric("lp_capable", function(solver) { standardGeneric("lp_capable") })

#' @describeIn Solver-capable A logical value indicating whether the solver is capable of solving second-order cone programs.
setGeneric("socp_capable", function(solver) { standardGeneric("socp_capable") })

#' @describeIn Solver-capable A logical value indicating whether the solver is capable of solving semidefinite programs.
setGeneric("sdp_capable", function(solver) { standardGeneric("sdp_capable") })

#' @describeIn Solver-capable A logical value indicating whether the solver is capable of solving exponential cone programs.
setGeneric("exp_capable", function(solver) { standardGeneric("exp_capable") })

#' @describeIn Solver-capable A logical value indicating whether the solver is capable of solving mixed-integer programs.
setGeneric("mip_capable", function(solver) { standardGeneric("mip_capable") })

setGeneric("status_map", function(solver, status) { standardGeneric("status_map") })

