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
setGeneric("is_dcp", function(object) { standardGeneric("is_dcp") })

setGeneric("size", function(object) { standardGeneric("size") })
setGeneric("dcp_curvature", function(monotonicity, func_curvature, arg_sign, arg_curvature) { standardGeneric("dcp_curvature") })
setGeneric("DCPAttr.mul_elemwise", function(lh_exp, rh_exp) { standardGeneric("DCPAttr.mul_elemwise") })

# Expression generic methods
setGeneric("value", function(object) { standardGeneric("value") })
setGeneric("value<-", function(object, value) { standardGeneric("value<-") })
setGeneric("save_value", function(object) { standardGeneric("save_value") })
setGeneric("get_data", function(object) { standardGeneric("get_data") })
setGeneric("init_dcp_attr", function(object) { standardGeneric("init_dcp_attr") })
setGeneric("curvature", function(object) { standardGeneric("curvature") })
setGeneric("is_scalar", function(object) { standardGeneric("is_scalar") })
setGeneric("is_vector", function(object) { standardGeneric("is_vector") })
setGeneric("is_matrix", function(object) { standardGeneric("is_matrix") })

setGeneric("variables", function(object) { standardGeneric("variables") })
setGeneric("parameters", function(object) { standardGeneric("parameters") })
setGeneric("constants", function(object) { standardGeneric("constants") })
setGeneric("domain", function(object) { standardGeneric("domain") })
setGeneric("validate_val", function(object, val) { standardGeneric("validate_value") })
setGeneric("canonical_form", function(object) { standardGeneric("canonical_form") })
setGeneric("canonicalize", function(object) { standardGeneric("canonicalize") })

# Positive definite inequalities
setGeneric("%>>%", function(e1, e2) { standardGeneric("%>>%") })
setGeneric("%<<%", function(e1, e2) { standardGeneric("%<<%") })

# Atom generic methods
setGeneric("validate_args", function(object) { standardGeneric("validate_args") })
setGeneric("shape_from_args", function(object) { standardGeneric("shape_from_args") })
setGeneric("sign_from_args", function(object) { standardGeneric("sign_from_args") })
setGeneric("func_curvature", function(object) { standardGeneric("func_curvature") })
setGeneric("monotonicity", function(object) { standardGeneric("monotonicity") })
setGeneric("get_data", function(object) { standardGeneric("get_data") })
setGeneric("name", function(object) { standardGeneric("name") })
setGeneric("to_numeric", function(object, values) { standardGeneric("to_numeric") })

setGeneric("Atom.dcp_curvature", function(curvature, args, monotonicities) { standardGeneric("Atom.dcp_curvature") })
setGeneric("graph_implementation", function(object, arg_objs, size, data) { standardGeneric("graph_implementation") })
setGeneric("sum_squares", function(expr) { standardGeneric("sum_squares") })

# Constraint generic methods
setGeneric("id", function(object) { standardGeneric("id") })
setGeneric("violation", function(object) { standardGeneric("violation") })
setGeneric("num_cones", function(object) { standardGeneric("num_cones") })
setGeneric("cone_size", function(object) { standardGeneric("cone_size") })
setGeneric("dual_variable", function(object) { standardGeneric("dual_variable") })
setGeneric("format_constr", function(object, eq_constr, leq_constr, dims, solver) { standardGeneric("format_constr") })
setGeneric("constr_type", function(object) { standardGeneric("constr_type") })

# Problem generic methods
setGeneric(".reset_cache", function(object) { standardGeneric(".reset_cache") })
setGeneric("cvxr_solve", function(object, solver, ignore_dcp, warm_start, verbose, parallel, ...) { standardGeneric("cvxr_solve") })

# Problem data generic methods
setGeneric("reset_param_data", function(object) { standardGeneric("reset_param_data") })
setGeneric(".dummy_constr", function(object) { standardGeneric(".dummy_constr") })
setGeneric("get_data", function(object) { standardGeneric("get_data") })
setGeneric("get_objective", function(object) { standardGeneric("get_objective") })
setGeneric("get_eq_constr", function(object) { standardGeneric("get_eq_constr") })
setGeneric("get_ineq_constr", function(object) { standardGeneric("get_ineq_constr") })
setGeneric(".init_matrix_cache", function(object, constraints, x_length) { standardGeneric(".init_matrix_cache") })
setGeneric(".lin_matrix", function(object, mat_cache, caching) { standardGeneric(".lin_matrix") })
setGeneric(".cache_to_matrix", function(object, mat_cache) { standardGeneric(".cache_to_matrix") })

# Solver generic methods
setGeneric("validate_solver", function(solver, constraints) { standardGeneric("validate_solver") })
setGeneric("validate_cache", function(solver, objective, constraints, cached_data) { standardGeneric("validate_cache") })
setGeneric("get_sym_data", function(solver, objective, constraints, cached_data) { standardGeneric("get_sym_data") })
setGeneric("get_matrix_data", function(solver, objective, constraints, cached_data) { standardGeneric("get_matrix_data") })
setGeneric("get_problem_data", function(solver, objective, constraints, cached_data) { standardGeneric("get_problem_data") })

setGeneric("matrix_intf", function(solver) { standardGeneric("matrix_intf") })
setGeneric("vec_intf", function(solver) { standardGeneric("vec_intf") })
setGeneric("split_constr", function(solver, constr_map) { standardGeneric("split_constr") })
setGeneric("cvxr_solve_int", function(solver, objective, constraints, cached_data, warm_start, verbose, solver_opts) { standardGeneric("cvxr_solve_int") })
setGeneric("format_results", function(solver, results_dict, data, cached_data) { standardGeneric("format_results") })

