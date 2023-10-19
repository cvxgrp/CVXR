######################################################################
#                      Solver utility functions.
######################################################################

#'
#' Stacks the values of the given variables.
#'
#' @param variables a list of variables.
#' @param default value to use when variable value is NA
#' @param byrow logical value indicating whether to fill the matrix by rows. Defaults to FALSE.
stack_vals <- function(variables, default, byrow = FALSE) {
  value <- list()
  for(variable in variables) {
    val <- value(variable)
    if(is.na(val))   # Unknown values.
      mat <- matrix(rep(default, size(variable)), ncol = 1)
    else
      mat <- matrix(val, ncol = 1, byrow = byrow)
    value <- c(value, list(mat))
  }
  # return(do.call("c", value))
  return(do.call("rbind", value))
}

expcone_permutor <- function(n_cones, exp_cone_order) {
  if(n_cones == 0)
    return(c())
  order <- rep(exp_cone_order, n_cones)   # e.g., c(1,0,2, 1,0,2, 1, 0,2, ...)
  offsets <- 3*do.call("c", lapply(seq(n_cones), function(i) { rep(i-1, 3) }))   # e.g., c(0,0,0, 3,3,3, 6,6,6, ...)
  perm <- order + offsets
  return(perm)
}

# #'
# #' Organize the constraints into a dictionary keyed by constraint names.
# #'
# #' @param constraints a list of constraints.
# #' @return A list of constraint types where constr_map[[cone_type]] maps to a list.
# group_constraints <- function(constraints) {
#  constr_map <- list()
#  for(c in constraints)
#    constr_map[[class(c)]] <- c(constr_map[[class(c)]], c)
#  return(constr_map)
# }

#'
#' Gets a specified value of a dual variable.
#'
#' @param result_vec A vector containing the dual variable values.
#' @param offset An offset to get correct index of dual values.
#' @param constraint A list of the constraints in the problem.
#' @return A list of a dual variable value and its offset.
extract_dual_value <- function(result_vec, offset, constraint) {
  value <- result_vec[seq(offset + 1, length.out = size(constraint))]
  if(size(constraint) == 1)
    value <- as.numeric(value)
  offset <- offset + size(constraint)
  return(list(value, offset))
}

#'
#' Gets the values of the dual variables.
#'
#' @param result_vec A vector containing the dual variable values.
#' @param parse_func Function handle for the parser.
#' @param constraints A list of the constraints in the problem.
#' @return A map of constraint ID to dual variable value.
get_dual_values <- function(result_vec, parse_func, constraints) {
  dual_vars <- list()
  offset <- 0
  for(constr in constraints) {
    # TODO: Reshape based on dual variable size.
    parsed <- parse_func(result_vec, offset, constr)
    dual_vars[[as.character(id(constr))]] <- parsed[[1]]
    offset <- parsed[[2]]
  }
  return(dual_vars)
}

#'
#' The ReductionSolver class.
#'
#' Generic interface for a solver that uses reduction semantics.
#'
#' @slot DIMS The key that maps to ConeDims in the data returned by perform(). There are separate ConeDims classes for cone programs vs. QPs. See matrix stuffing functions for details.
#' @slot MIP_CAPABLE Can the solver handle mixed-integer programs?
#' @rdname ReductionSolver-class
setClass("ReductionSolver", representation(DIMS = "character", VAR_ID = "character", DUAL_VAR_ID = "character", EQ_CONSTR = "character", NEQ_CONSTR = "character", MIP_CAPABLE = "logical"),   # Keys for inverse data. Internal use only!
                            prototype(DIMS = "dims", VAR_ID = "var_id", DUAL_VAR_ID = "dual_var_id", EQ_CONSTR = "eq_constr", NEQ_CONSTR = "other_constr", MIP_CAPABLE = FALSE), contains = "Reduction")

# Solver capabilities.
#' @param solver,object,x A \linkS4class{ReductionSolver} object.
#' @describeIn ReductionSolver Can the solver handle mixed-integer programs?
setMethod("mip_capable", "ReductionSolver", function(solver) { solver@MIP_CAPABLE })

#' @describeIn ReductionSolver Returns the name of the solver
setMethod("name", "ReductionSolver", function(x) { stop("Unimplemented") })

#' @describeIn ReductionSolver Imports the solver
setMethod("import_solver", "ReductionSolver", function(solver) { stop("Unimplemented") })

#' @describeIn ReductionSolver Is the solver installed?
setMethod("is_installed", "ReductionSolver", function(solver) {
  tryCatch({
    import_solver(solver)
  }, error = function(e) {
    solver_str <- ifelse(is(solver, "character"), solver, name(solver))
    warning("Encountered unexpected error importing solver ", solver_str)
  })
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ReductionSolver Solve a problem represented by data returned from perform.
setMethod("solve_via_data", "ReductionSolver", function(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts, solver_cache) {
  ## if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
  stop("Unimplemented")
})

#' @param problem A \linkS4class{Problem} object.
#' @describeIn ReductionSolver Solve the problem and return a Solution object.
setMethod("reduction_solve", "ReductionSolver", function(object, problem, warm_start, verbose, solver_opts) {
  res <- perform(object, problem)
  object <- res[[1]]
  data <- res[[2]]
  inv_data <- res[[3]]
  solution <- solve_via_data(object, data, warm_start, verbose, solver_opts)
  return(invert(object, solution, inv_data))
})

#'
#' The ConstantSolver class.
#'
ConstantSolver <- setClass("ConstantSolver", prototype(MIP_CAPABLE = TRUE), contains = "ReductionSolver")

#' @param problem A \linkS4class{Problem} object.
#' @describeIn ConstantSolver Is the solver capable of solving the problem?
setMethod("accepts", signature(object = "ConstantSolver", problem = "Problem"), function(object, problem) {
  return(length(variables(problem)) == 0)
})

#' @describeIn ConstantSolver Returns a list of the ConstantSolver, Problem, and an empty list.
setMethod("perform", signature(object = "ConstantSolver", problem = "Problem"), function(object, problem) {
  return(list(object, problem, list()))
})

#' @param solution A \linkS4class{Solution} object to invert.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn ConstantSolver Returns the solution.
setMethod("invert", signature(object = "ConstantSolver", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) {
  return(solution)
})

#' @describeIn ConstantSolver Returns the name of the solver.
setMethod("name", "ConstantSolver", function(x) { return("CONSTANT_SOLVER") })

#' @describeIn ConstantSolver Imports the solver.
setMethod("import_solver", "ConstantSolver", function(solver) { TRUE })

#' @describeIn ConstantSolver Is the solver installed?
setMethod("is_installed", "ConstantSolver", function(solver) { TRUE })


#' @param data Data for the solver.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ConstantSolver Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "ConstantSolver", function(object, data, warm_start, verbose, solver_opts, solver_cache = NULL) {
  ## if (missing(solver_cache)) solver_cache  <- new.env(parent=emptyenv())
  return(reduction_solve(object, data, warm_start, verbose, solver_opts))
})

#' @param problem A \linkS4class{Problem} object.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param solver_opts A list of Solver specific options
#' @describeIn ConstantSolver Solve the problem and return a \linkS4class{Solution} object.
setMethod("reduction_solve", "ConstantSolver", function(object, problem, warm_start, verbose, solver_opts) {
  if(all(sapply(problem@constraints, function(c) { !is.na(value(c)) })))
    return(Solution(OPTIMAL, value(problem@objective), list(), list(), list()))
  else
    return(Solution(INFEASIBLE, NA_real_, list(), list(), list()))
})

##################################################
#                 SOLVING CHAIN
##################################################

SolvingChain.is_lp <- function(object) {
  # Is problem a linear program?
  for(c in object@constraints) {
    if(!(is(c, "EqConstraint") || is(c, "ZeroConstraint")) || is_pwl(c@args[[1]]))
      return(FALSE)
  }
  
  for(var in variables(object)) {
    if(is_psd(var) || is_nsd(var))
      return(FALSE)
  }
  
  return(is_dcp(object) && is_pwl(object@objective@args[[1]]))
}

SolvingChain.solve_as_qp(problem, candidates) {
  conic_not_qp_slvs <- c()
  for(s in candidates$conic_solvers) {
    if(!(s %in% candidates$qp_solvers))
      conic_not_qp_slvs <- c(conic_not_qp_slvs, s)
  }
  
  if(SolvingChain.is_lp(problem) && length(conic_not_qp_slvs) > 0) {
    # OSQP can take many iterations for LPs; use a conic solver instead.
    # GUROBI and CPLEX QP/LP interfaces are more efficient -> Use them instead of conic if applicable.
    return(FALSE)
  }
  return(length(candidates$qp_solvers) > 0 && Qp2SymbolicQp.accepts(problem))
}

#'
#' Builds a chain that rewrites a problem into an intermediate representation suitable for numeric reductions.
#'
#' @param problem The problem for which to build a chain.
#' @param candidates List of candidate solvers divided into qp_solvers and conic_solvers.
#' @param gp If TRUE, the problem is parsed as a Disciplined Geometric Program instead of a Disciplined Convex Program
#' @return A list of Reduction objects that can be used to convert the problem to an intermediate form.
SolvingChain.reductions_for_problem_class(problem, candidates, gp = FALSE, solver_opts = NULL) {
  reductions <- list()
  
  # TODO: Handle boolean constraints.
  if(Complex2Real.accepts(problem))
    reductions <- c(reductions, list(Complex2Real()))
  if(gp)
    reductions <- c(reductions, list(Dgp2Dcp()))
  
  if(!gp && !is_dcp(problem)) {
    append <- build_non_disciplined_error_msg(problem, "DCP")
    if(is_dgp(problem))
      append <- paste(append, "However, the problem does follow DGP rules. Consider calling solve() with gp = TRUE", sep = "\n")
    else if(is_dqcp(problem))
      append <- paste(append, "However, the problem does follow DQCP rules. Consider calling solve() with qcp = TRUE", sep = "\n")
    stop(paste("Problem does not follow DCP rules. Specifically:", append, sep = "\n"))
  } else if(gp && !is_dgp(problem)) {
    append <- build_non_disciplined_error_msg(problem, "DGP")
    if(is_dcp(problem))
      append <- paste(append, "However, the problem does follow DCP rules. Consider calling solve() with gp = FALSE", sep = "\n")
    else if(is_dqcp(problem))
      append <- paste(append, "However, the problem does follow DQCP rules. Consider calling solve() with qcp = TRUE", sep = "\n")
    stop("Problem does not follow DGP rules. ", append)
  }
  
  # Dcp2Cone and Qp2SymbolicQp require problems to minimize their objectives.
  if(inherits(problem@objective, "Maximize"))
    reductions <- c(reductions, list(FlipObjective()))
  
  if(!is.null(solver_opts) && "use_quad_obj" %in% names(solver_opts))
    use_quad <- solver_opts$use_quad_obj
  else
    use_quad <- TRUE
  
  if(SolvingChain.solve_as_qp(problem, candidates) && use_quad)
    reductions <- c(reductions, list(CvxAttr2Constr(), Qp2SymbolicQp()))
  else {
    # Canonicalize into a conic problem.
    if(length(candidates$conic_solvers) == 0)
      stop("Problem could not be reduced to a QP, and no conic solvers exist among candidate solvers")
  }
  
  constr_types <- sapply(problem@constraints, class)
  if("FiniteSet" %in% constr_types)
    reductions <- c(reductions, list(Valinvec2mixedint()))
  
  return(reductions)
}

#'
#' Build a reduction chain from a problem to an installed solver.
#' 
#' Note that if the supplied problem has zero variables, the solver parameter 
#' will be ignored.
#'
#' @param problem The problem for which to build a chain.
#' @param candidates A list of candidate solvers.
#' @param gp If TRUE, the problem is parsed as a Disciplined Geometric Program instead of a Disciplined Convex Program.
#' @param enforce_dpp When TRUE, a DPPError will be throw when trying to parse a non-DPP problem (instead of just a warning). Defaults to FALSE.
#' @param ignore_dpp When TRUE, DPP problems will be treated as non-DPP, which may speed up compilation. Defaults to FALSE.
#' @param canon_backend Specifies which backend to use for canonicalization, which can affect compilation time. Defaults to NA, i.e., selecting the default backend ("CPP").
#' @param solver_opts List of additional arguments to pass to the solver. 
#' @return A \linkS4class{SolvingChain} that can be used to solve the problem.
construct_solving_chain <- function(problem, candidates, gp = FALSE, enforce_dpp = FALSE, ignore_dpp = FALSE, canon_backend = NA_character_, solver_opts = NULL) {
  if(length(variables(problem)) == 0)
    return(SolvingChain(reductions = list(ConstantSolver())))
  reductions <- SolvingChain.reductions_for_problem_class(problem, candidates, gp, solver_opts)
  
  # Process DPP status of the problem.
  dpp_context <- ifelse(gp, "dgp", "dcp")
  dpp_error_msg <- paste("You are solving a parameterized problem that is not DPP. ",
                         "Because the problem is not DPP, subsequent solves will not be "
                         "faster than the first one. For more information, see the "
                         "documentation on Discplined Parametrized Programming, at\n"
                         "\thttps://www.cvxpy.org/tutorial/advanced/index.html#"
                         "disciplined-parametrized-programming", sep = "")
  
  if(ignore_dpp || !is_dpp(problem, dpp_context)) {
    # No warning for ignore_dpp.
    if(ignore_dpp)
      reductions <- c(list(EvalParams()), reductions)
    else if(!enforce_dpp) {
      warning(dpp_error_msg)
      reductions <- c(list(EvalParams()), reductions)
    } else
      stop(dpp_error_msg)
  } else if(any(sapply(parameters(problem), is_complex)))
    reductions <- c(list(EvalParams()), reductions)
  else {   # Compilation with DPP.
    n_parameters <- sum(parameters(problem), function(param) { prod(dim(param)) })
    if(n_parameters >= PARAM_THRESHOLD)
      warning("Your problem has too many parameters for efficient DPP compilation. We suggest setting ignore_dpp = TRUE")
  }
  
  # Conclude with matrix stuffing; choose one of the following paths:
  #   1) QpMatrixStuffing -> [a QpSolver]
  #   2) ConeMatrixStuffing -> [a ConicSolver]
  if(!is.null(solver_opts) && "use_quad_obj" %in% names(solver_opts))
    use_quad <- solver_opts$use_quad_obj
  else
    use_quad <- TRUE
  
  if(SolvingChain.solve_as_qp(problem, candidates) && use_quad) {
    # Canonicalize as a QP.
    solver <- candidates$qp_solvers[[1]]
    solver_instance <- SOLVER_MAP_QP[[solver]]
    reductions <- c(reductions, list(QpMatrixStuffing(canon_backend = canon_backend), solver_instance))
    return(SolvingChain(reductions = reductions))
  }
  
  # Canonicalize as a cone program.
  if(length(candidates$conic_solvers) == 0)
    stop(paste("Problem could not be reduced to a QP, and no conic solvers exist among candidate solvers (",
               paste(unlist(candidates), collapse = ","), ")", sep = ""))
  
  # We use constr_types to infer an incomplete list of cones that the solver will need after canonicalization.
  constr_types <- sapply(problem@constraints, class)
  constr_types <- unique(constr_types)
  ex_cos <- constr_types[constr_types %in% EXOTIC_CONES]   # The way we populate ex_cos will need to change if and when we have atoms that require exotic cones.
  approx_cos <- constr_types[constr_types %in% APPROX_CONES]
  for(co in ex_cos) {
    sim_cos <- EXOTIC_CONES[[co]]   # Get the set of required simple cones.
    if(!(app_cos %in% constr_types))
      constr_types <- c(constr_types, app_cos)
    constr_types <- constr_types[which(constr_types != co)]
  }
  
  # We now go over individual elementary cones support by CVXR (SOC, ExpCone, 
  # NonNegConstraint, ZeroConstraint, PSDConstraint, PowCone3D) and check if
  # they've appeared in constr_types or if the problem has an atom requiring 
  # that cone.
  cones <- c()
  atoms <- atoms(problem)
  if("SOC" %in% constr_types || any(atoms %in% SOC_ATOMS))
    cones <- c(cones, "SOC")
  if("ExpCone" %in% constr_types || any(atoms %in% EXP_ATOMS))
    cones <- c(cones, "ExpCone")
  if(any(sapply(c("IneqConstraint", "NonPosConstraint", "NonNegConstraint"), function(t) { t %in% constr_types } )) || any(atoms %in% NONPOS_ATOMS))
    cones <- c(cones, "NonNegConstraint")
  if("EqConstraint" %in% constr_types || "ZeroConstraint" %in% constr_types)
    cones <- c(cones, "ZeroConstraint")
  if("PSDConstraint" %in% constr_types || any(atoms %in% PSD_ATOMS) || any(variables(problem), function(v) { is_psd(v) || is_nsd(v) }))
    cones <- c(cones, "PSDConstraint")
  if("PowCone3D" %in% constr_types) {
    # If we add in atoms that specifically use the 3D power cone (rather than 
    # the ND power cone), then we'll need to check for those atoms here as well.
    cones <- c(cones, "PowCone3D")
  }

  # Here, we make use of the observation that canonicalization only
  # increases the number of constraints in our problem.
  has_constr <- length(cones) > 0 || length(problem@constraints) > 0
  
  for(solver %in% candidates$conic_solvers) {
    solver_instance <- SOLVER_MAP_CONIC[[solver]]
    
    # Cones supported for MI problems may differ from non MI.
    if(is_mixed_integer(problem))
      supported_constraints <- mi_supported_constraints(solver_instance)
    else
      supported_constraints <- supported_constraints(solver_instance)
    
    if(all(cones %in% supported_constraints) && (has_constr || !requires_constr(solver_instance))) {
      if(length(ex_cost) > 0)
        reductions <- c(reductions, list(Exotic2Common()))
      if("RelEntrConeQuad" %in% approx_cos || "OpRelEntrConeQuad" %in% approx_cos)
        reductions <- c(reductions, list(QuadApprox()))
      
      # Should the objective be canonicalized to a quadratic?
      if(!is.null(solver_opts) && "use_quad_obj" %in% names(solver_opts))
        use_quad_obj <- solver_opts$use_quad_obj
      else
        use_quad_obj <- TRUE
      quad_obj <- use_quad_obj && supports_quad_obj(solver_instance) && has_quadratic_term(problem@objective@expr)
      reductions <- c(reductions, list(Dcp2Cone(quad_obj = quad_obj), CvxAttr2Constr(), ConeMatrixStuffing(quad_obj = quad_obj, canon_backend = canon_backend), solver_instance))
      return(SolvingChain(reductions = reductions))
    }
  }
  
  stop(paste("Either candidate conic solvers (",
       paste(candidates$conic_solvers, sep = " ", collapse = ","), ") do not support the cones output by the problem (",
       paste(cones, sep = " ", collapse = ","), "), or there are not enough constraints in the problem.", sep = ""))
}

setClassUnion("ReductionSolverORNULL", c("ReductionSolver", "NULL"))

#'
#' The SolvingChain class.
#'
#' This class represents a reduction chain that ends with a solver.
#'
#' @rdname SolvingChain-class
.SolvingChain <- setClass("SolvingChain", representation(solver = "ReductionSolverORNULL"), prototype(solver = NULL), contains = "Chain")
SolvingChain <- function(problem = NULL, reductions = list()) { .SolvingChain(problem = problem, reductions = reductions) }

setMethod("initialize", "SolvingChain", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  if(length(.Object@reductions) == 0)
    stop("Solving chains must terminate with a ReductionSolver")

  last <- .Object@reductions[[length(.Object@reductions)]]
  if(!is(last, "ReductionSolver"))
    stop("Solving chains must terminate with a ReductionSolver.")
  .Object@solver <- last
  return(.Object)
})

# Create and return a new SolvingChain by concatenating chain with this instance.
#' @param object A \linkS4class{SolvingChain} object.
#' @param chain A \linkS4class{Chain} to prepend.
#' @describeIn SolvingChain Create and return a new SolvingChain by concatenating chain with this instance.
setMethod("prepend", signature(object = "SolvingChain", chain = "Chain"), function(object, chain) {
  SolvingChain(reductions = c(chain@reductions, object@reductions))
})

#' @param problem The problem to solve.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param solver_opts A list of Solver specific options
#' @describeIn SolvingChain Applies each reduction in the chain to the problem, solves it,
#' and then inverts the chain to return a solution of the supplied problem.
setMethod("reduction_solve", signature(object = "SolvingChain", problem = "Problem"), function(object, problem, warm_start, verbose, solver_opts) {
  tmp <- perform(object, problem)
  object <- tmp[[1]]
  data <- tmp[[2]]
  inverse_data <- tmp[[3]]
  solution <- solve_via_data(object@solver, data, warm_start, verbose, solver_opts)
  return(invert(object, solution, inverse_data))
})

#' @param data Data for the solver.
#' @describeIn SolvingChain Solves the problem using the data output by the an apply invocation.
setMethod("reduction_solve_via_data", "SolvingChain", function(object, problem, data, warm_start, verbose, solver_opts = list()) {
  return(solve_via_data(object@solver, data, warm_start, verbose, solver_opts, problem@.solver_cache))
})

#'
#' Construct Intermediate Chain
#'
#' Builds a chain that rewrites a problem into an intermediate representation suitable for numeric reductions.
#'
#' @param problem The problem for which to build a chain.
#' @param candidates A list of candidate solvers.
#' @param gp A logical value indicating whether the problem is a geometric program.
#' @return A \linkS4class{Chain} object that can be used to convert the problem to an intermediate form.
setMethod("construct_intermediate_chain", signature(problem = "Problem", candidates = "list"), function(problem, candidates, gp = FALSE) {
  reductions <- list()
  if(length(variables(problem)) == 0)
    return(Chain(reductions = reductions))

  # TODO: Handle boolean constraints.
  if(Complex2Real.accepts(problem))
    reductions <- c(reductions, Complex2Real())
  if(gp)
    reductions <- c(reductions, Dgp2Dcp())
  
  if(!gp && !is_dcp(problem)) {
    append <- build_non_disciplined_error_msg(problem, "DCP")
    if(is_dgp(problem))
      append <- paste(append, "However, the problem does follow DGP rules. Consider calling solve() with gp = TRUE", sep = "\n")
    else if(is_dqcp(problem))
      append <- paste(append, "However, the problem does follow DQCP rules. Consider calling solve() with qcp = TRUE", sep = "\n")
    stop(paste("Problem does not follow DCP rules. Specifically:", append, sep = "\n"))
  } else if(gp && !is_dgp(problem)) {
    append <- build_non_disciplined_error_msg(problem, "DGP")
    if(is_dcp(problem))
      append <- paste(append, "However, the problem does follow DCP rules. Consider calling solve with gp = FALSE", sep = "\n")
    else if(is_dqcp(problem))
      append <- paste(append, "However, the problem does follow DQCP rules. Consider calling solve() with qcp = TRUE", sep = "\n")
    stop(paste("Problem does not follow DGP rules.", append))
  }

  # Dcp2Cone and Qp2SymbolicQp require problems to minimize their objectives.
  if(inherits(problem@objective, "Maximize"))
    reductions <- c(reductions, FlipObjective())

  # First, attempt to canonicalize the problem to a linearly constrained QP.
  if(length(candidates$qp_solvers) > 0 && Qp2SymbolicQp.accepts(problem)) {
    reductions <- c(reductions, list(CvxAttr2Constr(), Qp2SymbolicQp()))
    return(Chain(reductions = reductions))
  }

  # Canonicalize it to conic problem.
  if(length(candidates$conic_solvers) == 0)
    stop("Problem could not be reduced to a QP, and no conic solvers exist among candidate solvers")
  reductions <- c(reductions, list(Dcp2Cone(), CvxAttr2Constr()))
  return(Chain(reductions = reductions))
})

################################################
#                 BISECTION
################################################
Bisection.lower_problem <- function(problem) {
  # Evaluates lazy constraints
  Problem(Minimize(0), c(problem@constraints, lapply(problem@.lazy_constraints, function(con) { con() })))
}

Bisection.infeasible <- function(problem, result) {
  # is.null(problem) || problem@status %in% c(INFEASIBLE, INFEASIBLE_INACCURATE)
  is.null(problem) || result$status %in% c(INFEASIBLE, INFEASIBLE_INACCURATE)
}

# TODO: Eliminate the need for reconstructing the problem (via Bisection.lower_problem). 
# Right now, some constraints are lazy, namely the callable constraints, because for 
# certain values of the parameters they become invalidated. Either support lazy 
# constraints that do not need recompilation, or (if possible) find rewritings that 'just work'.
Bisection.find_bisection_interval <- function(problem, t, solver = NULL, low = NA_real_, high = NA_real_, max_iters = 100) {
  # Finds an interval for bisection
  if(is.na(low))
    low <- ifelse(is_nonneg(t), 0, -1)
  if(is.na(high))
    high <- ifelse(is_nonpos(t), 0, 1)
  
  infeasible_low <- is_nonneg(t)
  feasible_high <- is_nonpos(t)
  for(i in seq(max_iters)) {
    if(!feasible_high) {
      value(t) <- high
      lowered <- Bisection.lower_problem(problem)
      result <- solve(lowered, solver = solver)
      if(Bisection.infeasible(lowered, result)) {
        low <- high
        high <- 2*high
        next
      } else if(result$status %in% SOLUTION_PRESENT)
        feasible_high <- TRUE
      else
        stop("Solver failed with status ", result$status)
    }
    
    if(!infeasible_low) {
      value(t) <- low
      lowered <- Bisection.lower_problem(problem)
      result <- solve(lowered, solver = solver)
      if(Bisection.infeasible(lowered, result))
        infeasible_low <- TRUE
      else if(result$status %in% SOLUTION_PRESENT) {
        high <- low
        low <- 2*low
        next
      } else
        stop("Solver failed with status ", result$status)
    }
    
    if(infeasible_low && feasible_high)
      return(c(low, high))
  }
  
  stop("Unable to find suitable interval for bisection; your problem may be unbounded.")
}

Bisection.bisect_int <- function(problem, solver, t, low, high, tighten_lower, tighten_higher, eps = 1e-6, verbose = FALSE, max_iters = 100) {
  # Bisect problem on the parameter t
  verbose_freq <- 5
  soln <- NA_real_
  
  for(i in seq(max_iters)) {
    if(low > high)
      stop("low must be less than or equal to high")
    
    if(!is.na(soln) && (high - low) <= eps) {
      # the previous iteration might have been infeasible, but
      # the tighten* functions might have narrowed the interval
      # to the optimal value in the previous iteration (hence the
      # soln is not None check)
      return(c(soln, low, high))
    }
    
    query_pt <- (low + high) / 2.0
    if(verbose && i %% verbose_freq == 0) {
      print(paste("(iteration ", i, ") lower bound: ", low, sep = ""))
      print(paste("(iteration ", i, ") upper bound: ", high, sep = ""))
      print(paste("(iteration ", i, ") query point: ", query_pt, sep = ""))
    }
    value(t) <- query_pt
    lowered <- Bisection.lowered_problem(problem)
    result <- solve(lowered, solver = solver)
    
    if(Bisection.infeasible(lowered, result)) {
      if(verbose && i %% verbose_freq == 0)
        print(paste("(iteration ", i, ") query was infeasible"))
      low <- tighten_lower(query_pt)
    } else if(result$status %in% SOLUTION_PRESENT) {
      if(verbose && i %% verbose_freq == 0)
        print(paste("(iteration ", i, ") query was feasible with solution ", result$solution, sep = ""))
      soln <- result$solution
      high <- tighten_higher(query_pt)
    } else {
      if(verbose)
        print("Aborting; the solver failed...")
      stop("Solver failed with status ", result$status)
    }
  }
}

#'
#' Bisection on a one-parameter family of DCP problems.
#' 
#' Bisects on a one-parameter family of DCP problems emitted by Dqcp2Dcp.
#' 
#' @param problem A \linkS4class{Problem} emitted by Dqcp2Dcp.
#' @param solver The solver to use for bisection.
#' @param low (Optional) Lower bound for bisection.
#' @param high (Optional) Upper bound for bisection.
#' @param eps Terminate bisection when width of interval is strictly less than eps.
#' @param verbose A logical value indicating whether to print verbose output related to bisection.
#' @param max_iters The maximum number of iterations to run the bisection.
#' @return A \linkS4class{Solution} object.
setMethod("bisect", "Problem", function(problem, solver = NULL, low = NA_real_, high = NA_real_, eps = 1e-6, verbose = FALSE, max_iters = 100, max_iters_interval_search = 100) {
  if(!.hasSlot(problem, ".bisection_data"))
    stop("bisect only accepts problems emitted by Dqcp2Dcp (with .bisection_data slot)")
  tmp <- problem@.bisection_data
  feas_problem <- tmp[[1]]
  t <- tmp[[2]]
  tighten_lower <- tmp[[3]]
  tighten_higher <- tmp[[4]]
  
  if(verbose)
    cat("\n******************************************************",
        "**************************\n",
        "Preparing to bisect problem\n\n", as.character(Bisection.lower_problem(problem)), "\n", sep = "")
  
  lowered_feas <- Bisection.lower_problem(feas_problem)
  result <- solve(lowered_feas, solver = solver)
  if(Bisection.infeasible(lowered_feas, result)) {
    if(verbose)
      print("Problem is infeasible")
    return(failure_solution(INFEASIBLE))
  }
  
  if(is.na(low) || is.na(high)) {
    if(verbose)
      print("Finding interval for bisection...")
    tmp <- Bisection.find_bisection_interval(problem, t, solver, low, high, max_iters_interval_search)
    low <- tmp[[1]]
    high <- tmp[[2]]
  }
  
  if(verbose) {
    print(paste("initial lower bound:", low))
    print(paste("initial upper bound:", high))
  }
  
  tmp <- Bisection.bisect_int(problem, solver, t, low, high, tighten_lower, tighten_higher, eps, verbose, max_iters)
  soln <- tmp[[1]]
  low <- tmp[[2]]
  high <- tmp[[3]]
  soln@opt_val <- (low + high) / 2.0
  if(verbose)
    cat("Bisection completed, with lower bound ", low, " and upper bound ", high,
          "\n******************************************",
          "**************************************\n", sep = "")
  return(soln)
})

################################################
#             MATRIX COMPRESSION
################################################
# TODO: This uses matrix pointers (indptr) in Python, so it should be implemented in C++.

################################################
#                 KKT SOLVER
################################################
# A custom KKT solver for CVXOPT that can handle redundant constraints.
# Uses regularization and iterative refinement.

REG_EPS <- 1e-9   # Regularization constant.

setup_ldl_factor <- function(c, G, h, dims, A, b) {
  # The meanings of arguments in this function are identical to those of the
  # function cvxopt.solvers.conelp. Refer to CVXOPT documentation
  #
  # https://cvxopt.org/userguide/coneprog.html#linear-cone-programs
  #
  # for more information.
  #
  # Note: CVXOPT allows G and A to be passed as dense matrix objects. However,
  # this function will only ever be called with spmatrix objects. If creating
  # a custom kktsolver of your own, you need to conform to this sparse matrix
  # assumption.
  
  factor <- kkt_ldl(G, dims, A)
  return(factor)
}

kkt_ldl <- function(G, dims, A) {
  # Returns a function handle "factor", which conforms to the CVXOPT
  # custom KKT solver specifications:
  #
  #     https://cvxopt.org/userguide/coneprog.html#exploiting-structure.
  #
  # For convenience, we provide a short outline for how this function works.
  #
  # First, we allocate workspace for use in "factor". The factor function is
  # called with data (H, W). Once called, the factor function computes an LDL
  # factorization of the 3 x 3 system:
  #
  #     [ H           A'   G'*W^{-1}  ]
  #     [ A           0    0          ].
  #     [ W^{-T}*G    0   -I          ]
  #
  # Once that LDL factorization is computed, "factor" constructs another
  # inner function, called "solve". The solve function uses the newly
  # constructed LDL factorization to compute solutions to linear systems of
  # the form
  #
  #     [ H     A'   G'    ]   [ ux ]   [ bx ]
  #     [ A     0    0     ] * [ uy ] = [ by ].
  #     [ G     0   -W'*W  ]   [ uz ]   [ bz ]
  #
  # The factor function concludes by returning a reference to the solve function.
  #
  # Note: In the 3 x 3 system, H is n x n, A is p x n, and G is N x n, where
  # N = dims[['l']] + sum(dims[['q']]) + sum(sapply(dims[['s']], function(k) { k^2 })). 
  # For cone programs, H is the zero matrix.
  
  p <- nrow(A)
  n <- ncol(A)
  ldK <- n + p + dims$l + sum(dims$q) + sum(sapply(dims$s, function(k) { as.integer(k*(k+1)/2) } ))
  
  stop("CVXOPT KKT solver is unimplemented. Need to write most of the code in C++.")
  
  # TODO: These should all be CVXOPT matrix objects.
  K <- matrix(0, nrow = ldK, ncol = ldK)
  ipiv <- matrix(0, nrow = ldK, ncol = 1)
  u <- matrix(0, nrow = ldK, ncol = 1)
  g <- matrix(0, nrow = nrow(G), ncol = 1)
  
  factor <- function(W, H = NA_real_) {
    K <- 0*K
    if(!is.na(H))
      K[1:n, 1:n] <- H
    K[(n+1):(n+p+1), 1:n] <- A
    
    for(k in seq(n)) {
      g <- G[,k]
      # TODO: These should be calls to functions in CVXOPT (cvxopt.misc).
      # scale(g, W, trans='T', inverse='I')
      # pack(g, K, dims, 0, offsety=k*ldK + n + p)
    }
    K[seq((ldK+1)*(p+n) + 1, nrow(K), by = ldK+1)] <- -1.0
    
    # Add positive regularization in 1x1 block and negative in 2x2 block.
    K[seq(1, (ldK+1)*n, by = ldK+1)] <- K[seq(1, (ldK+1)*n, by = ldK+1)] + REG_EPS
    K[seq((ldK+1)*n + 1, nrow(K), by = ldK+1)] <- K[seq((ldK+1)*n + 1, nrow(K), by = ldK+1)] - REG_EPS
    # TODO: lapack.sytrf(K, ipiv)
    
    solve <- function(x, y, z) {
      # Solve
      #
      #     [ H          A'   G'*W^{-1}  ]   [ ux   ]   [ bx        ]
      #     [ A          0    0          ] * [ uy   [ = [ by        ]
      #     [ W^{-T}*G   0   -I          ]   [ W*uz ]   [ W^{-T}*bz ]
      #
      # and return ux, uy, W*uz.
      #
      # On entry, x, y, z contain bx, by, bz.  On exit, they contain
      # the solution ux, uy, W*uz.
      
      # TODO: This function needs to be implemented using CVXOPT/BLAS/LAPACK functions.
      # blas.copy(x, u)
      # blas.copy(y, u, offsety=n)
      # scale(z, W, trans='T', inverse='I')
      # pack(z, u, dims, 0, offsety=n + p)
      # lapack.sytrs(K, ipiv, u)
      # blas.copy(u, x, n=n)
      # blas.copy(u, y, offsetx=n, n=p)
      # unpack(u, z, dims, 0, offsetx=n + p)
    }
    
    return(solve)
  }
  
  return(factor)
}

