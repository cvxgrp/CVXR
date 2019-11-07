#'
#' The Solver class.
#'
#' This virtual class represents the generic interface for a solver.
#'
#' @name Solver-class
#' @aliases Solver
#' @rdname Solver-class
Solver <- setClass("Solver", contains = "VIRTUAL")

#'
#' Choose a Solver
#'
#' Determines the appropriate solver.
#'
#' @param constraints A list of canonicalized constraints.
#' @return A \linkS4class{Solver} object.
#' @rdname Solver-choose_solver
Solver.choose_solver <- function(constraints) {
  constr_map <- SymData.filter_constraints(constraints)
  # If no constraints, use ECOS.
  if(length(constraints) == 0)
    return(ECOS())
  # If mixed integer constraints, use ECOS_BB.
  else if(length(constr_map[[BOOL_MAP]]) > 0 || length(constr_map[[INT_MAP]]) > 0)
    return(ECOS_BB())
  # If SDP, defaults to CVXOPT.
  else if(length(constr_map[[SDP_MAP]]) > 0) {
    # if(requireNamespace("cccp")) {
    #  return(CVXOPT())
    # } else
      return(SCS())
  }
  # Otherwise use ECOS
  else
    return(ECOS())
}

setMethod("is_installed", "Solver", function(solver) {
  tryCatch({
    import_solver(solver)
  }, error = function(e) {
    FALSE
  })
})

#' @param solver A \linkS4class{Solver} object.
#' @describeIn Solver Verify the solver can solve the problem.
setMethod("validate_solver", "Solver", function(solver, constraints) {
  # Check the solver is installed.
  if(!import_solver(solver))
    stop("The solver ", name(solver), " is not installed.")

  # Check the solver can solve the problem.
  constr_map <- SymData.filter_constraints(constraints)

  if( (length(constr_map[[BOOL_MAP]]) > 0 || length(constr_map[[INT_MAP]]) > 0) && !mip_capable(solver))
    Solver._reject_problem(solver, "it cannot solve mixed-integer problems")
  else if(length(constr_map[[SDP_MAP]]) > 0 && !sdp_capable(solver))
    Solver._reject_problem(solver, "it cannot solve semidefinite problems")
  else if(length(constr_map[[EXP_MAP]]) > 0 && !exp_capable(solver))
    Solver._reject_problem(solver, "it cannot solve exponential cone problems")
  else if(length(constr_map[[SOC_MAP]]) > 0 && !socp_capable(solver))
    Solver._reject_problem(solver, "it cannot solve second-order cone problems")
  else if(length(constraints) == 0 && name(solver) %in% c("SCS", "GLPK"))
    Solver._reject_problem(solver, "it cannot solve unconstrained problems")
})

#
# Reject Problem
#
# Raise an error indicating that the solver cannot solve a problem.
#
# @param solver A \linkS4class{Solver} object.
# @param reason A short description of the reason the problem cannot be solved by this solver.
# @rdname Solver-reject_problem
Solver._reject_problem <- function(solver, reason) {
  message <- paste("The solver", name(solver), "cannot solve the problem because", reason)
  stop(message)
}

setMethod("validate_cache", "Solver", function(solver, objective, constraints, cached_data) {
  prob_data <- cached_data[[name(solver)]]
  if(!is.null(prob_data@sym_data) && (!isTRUE(all.equal(objective, prob_data@sym_data@objective)) ||
                                      !isTRUE(all.equal(constraints, prob_data@sym_data@constraints)))) {
    prob_data@sym_data <- NULL
    prob_data@matrix_data <- NULL
  }
  cached_data[[name(solver)]] <- prob_data
  cached_data
})

setMethod("get_sym_data", "Solver", function(solver, objective, constraints, cached_data) {
  cached_data <- validate_cache(solver, objective, constraints, cached_data)
  prob_data <- cached_data[[name(solver)]]
  if(is.null(prob_data@sym_data))
    prob_data@sym_data <- SymData(objective, constraints, solver)
  cached_data[[name(solver)]] <- prob_data
  cached_data
})

setMethod("get_matrix_data", "Solver", function(solver, objective, constraints, cached_data) {
  cached_data <- get_sym_data(solver, objective, constraints, cached_data)
  sym_data <- cached_data[[name(solver)]]@sym_data
  prob_data <- cached_data[[name(solver)]]
  if(is.null(prob_data@matrix_data))
    prob_data@matrix_data <- MatrixData(sym_data, solver)   # TODO: Update this constructor for nonlinear constraints
  cached_data[[name(solver)]] <- prob_data
  cached_data
})

setMethod("Solver.get_problem_data", "Solver", function(solver, objective, constraints, cached_data) {
  cached_data <- get_sym_data(solver, objective, constraints, cached_data)
  sym_data <- cached_data[[name(solver)]]@sym_data
  cached_data <- get_matrix_data(solver, objective, constraints, cached_data)
  matrix_data <- cached_data[[name(solver)]]@matrix_data

  data <- list()
  obj <- get_objective(matrix_data)
  eq <- get_eq_constr(matrix_data)
  ineq <- get_ineq_constr(matrix_data)

  data[[C_KEY]] <- obj[[1]]
  data[[OFFSET]] <- obj[[2]]
  data[[A_KEY]] <- eq[[1]]
  data[[B_KEY]] <- eq[[2]]
  data[[G_KEY]] <- ineq[[1]]
  data[[H_KEY]] <- ineq[[2]]
  data[[DIMS]] <- sym_data@.dims

  conv_idx <- Solver._noncvx_id_to_idx(data[[DIMS]], sym_data@.var_offsets, sym_data@.var_sizes)
  data[[DIMS]] <- conv_idx$dims
  data[[BOOL_IDX]] <- unique(conv_idx$bool_idx)
  data[[INT_IDX]] <- unique(conv_idx$int_idx)
  data
})

#' @describeIn Solver A logical value indicating whether nonlinear constraints are needed.
setMethod("nonlin_constr", "Solver", function(solver) { FALSE })

#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @describeIn Solver Call the solver on the canonicalized problem.
setMethod("Solver.solve", "Solver", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) { stop("Unimplemented") })

#' @param results_dict A list containing the solver output.
#' @param data A list containing information about the problem.
#' @describeIn Solver Convert raw solver output into standard list of results.
setMethod("format_results", "Solver", function(solver, results_dict, data, cached_data) { stop("Unimplemented") })

# Is the problem a mixed-integer program?
Solver.is_mip <- function(data) {
  length(data[[BOOL_IDX]]) > 0 || length(data[[INT_IDX]]) > 0
}

#
# Non-Convex ID to Index
#
# Converts the non-convex constraint variable IDs in dims into indices.
#
# @param dims The dimensions of the cones.
# @param var_offsets A list mapping variabld ID to horizontal offset.
# @param var_sizes A list mapping variable ID to  variable dimensions.
# @return A list of indices for the boolean and integer variables.
# @rdname Solver-noncvx_id_to_idx
Solver._noncvx_id_to_idx <- function(dims, var_offsets, var_sizes) {
  if(BOOL_IDS %in% names(dims)) {
    bool_idx <- lapply(dims[[BOOL_IDS]], function(var_id) {
      var_id <- as.character(var_id)
      offset <- var_offsets[[var_id]]
      size <- var_sizes[[var_id]]
      offset + seq(1, size[1]*size[2], by = 1)
    })
    bool_idx <- unlist(bool_idx)   # Should I call unique on this? Seeing dupes in dims[[BOOL_IDS]]
    dims[[BOOL_IDS]] <- NULL
  } else
    bool_idx <- integer(0)

  if(INT_IDS %in% names(dims)) {
    int_idx <- lapply(dims[[INT_IDS]], function(var_id) {
      var_id <- as.character(var_id)
      offset <- var_offsets[[var_id]]
      size <- var_sizes[[var_id]]
      offset + seq(1, size[1]*size[2], by = 1)
    })
    int_idx <- unlist(int_idx)   # Should I call unique on this? Seeing dupes in dims[[INT_IDS]]
    dims[[INT_IDS]] <- NULL
  } else
    int_idx <- integer(0)

  list(dims = dims, bool_idx = bool_idx, int_idx = int_idx)
}

#'
#' The ECOS class.
#'
#' This class is an interface for the ECOS solver.
#'
#' @references A. Domahidi, E. Chu, and S. Boyd. "ECOS: An SOCP solver for Embedded Systems." \emph{Proceedings of the European Control Conference}, pp. 3071-3076, 2013.
#' @seealso \code{\link[ECOSolveR]{ECOS_csolve}} and the \href{https://github.com/embotech/ecos}{ECOS Official Repository}.
#' @name ECOS-class
#' @aliases ECOS
#' @rdname ECOS-class
setClass("ECOS", contains = "Solver")

#' @rdname ECOS-class
#' @examples
#' ecos <- ECOS()
#' lp_capable(ecos)
#' sdp_capable(ecos)
#' socp_capable(ecos)
#' exp_capable(ecos)
#' mip_capable(ecos)
#' @export
ECOS <- function() {
    new("ECOS")
    ##ECOS$new()
}

# ECOS capabilities
#' @param object,solver An \linkS4class{ECOS} object.
#' @describeIn ECOS ECOS can handle linear programs.
setMethod("lp_capable", "ECOS", function(solver) { TRUE })

#' @describeIn ECOS ECOS can handle second-order cone programs.
setMethod("socp_capable", "ECOS", function(solver) { TRUE })

#' @describeIn ECOS ECOS cannot handle semidefinite programs.
setMethod("sdp_capable", "ECOS", function(solver) { FALSE })

#' @describeIn ECOS ECOS can handle exponential cone programs.
setMethod("exp_capable", "ECOS", function(solver) { TRUE })

#' @describeIn ECOS ECOS cannot handle mixed-integer programs.
setMethod("mip_capable", "ECOS", function(solver) { FALSE })

#'
#' ECOS Status Map
#'
#' Map of ECOS status to CVXR status.
#'
#' @param solver A \linkS4class{ECOS} object.
#' @param status An exit code returned by ECOS:
#' \describe{
#'    \item{ECOS_OPTIMAL (0)}{Problem solved to optimality.}
#'    \item{ECOS_PINF (1)}{Found certificate of primal infeasibility.}
#'    \item{ECOS_DINF (2)}{Found certificate of dual infeasibility.}
#'    \item{ECOS_INACC_OFFSET (10)}{Offset exitflag at inaccurate results.}
#'    \item{ECOS_MAXIT (-1)}{Maximum number of iterations reached.}
#'    \item{ECOS_NUMERICS (-2)}{Search direction unreliable.}
#'    \item{ECOS_OUTCONE (-3)}{\eqn{s} or \eqn{z} got outside the cone, numerics?}
#'    \item{ECOS_SIGINT (-4)}{Solver interrupted by a signal/ctrl-c.}
#'    \item{ECOS_FATAL (-7)}{Unknown problem in solver.}
#' }
#' @return A string indicating the status, either "optimal", "infeasible", "unbounded", "optimal_inaccurate", "infeasible_inaccurate", "unbounded_inaccurate", or "solver_error".
#' @docType methods
#' @rdname ECOS-status_map
setMethod("status_map", "ECOS", function(solver, status) {
    if(status == 0) {
        OPTIMAL
    } else if(status == 1) {
        INFEASIBLE
    } else if(status == 2) {
        UNBOUNDED
    } else if(status == 10) {
        OPTIMAL_INACCURATE
    } else if(status == 11) {
        INFEASIBLE_INACCURATE
    } else if(status == 12) {
        UNBOUNDED_INACCURATE
    } else if(status %in% c(-1, -2, -3, -4, -7)) {
        SOLVER_ERROR
    } else stop("ECOS status unrecognized: ", status)
})

#' @describeIn ECOS The name of the solver.
#' @export
setMethod("name", "ECOS", function(object) { ECOS_NAME })

#' @describeIn ECOS Imports the ECOSolveR library.
setMethod("import_solver", "ECOS", function(solver) { TRUE }) ## ECOSolveR is in imports!

setMethod("split_constr", "ECOS", function(solver, constr_map) {
  list(eq_constr = constr_map[[EQ_MAP]], ineq_constr = constr_map[[LEQ_MAP]], nonlin_constr = list())
})

#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @importFrom ECOSolveR ECOS_csolve ecos.control
#' @describeIn ECOS Call the solver on the canonicalized problem.
setMethod("Solver.solve", "ECOS", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)
  data[[DIMS]]['e'] <- data[[DIMS]][[EXP_DIM]]

  # TODO: Naras please fix this by making ECOSolveR handle type conversion, e.g. logical -> integer
  ## Naras Answer: Fixed in ECOSolveR version >= 0.4
  ## if(prod(dim(data[[G_KEY]])) == 0) data[[G_KEY]] <- NULL
  ## if(prod(dim(data[[A_KEY]])) == 0) data[[A_KEY]] <- NULL
  ## data[[DIMS]] <- lapply(data[[DIMS]], function(dim) { as.integer(dim) })
  solver_opts <- ECOSolveR::ecos.control()
  solver_opts$VERBOSE <- as.integer(verbose)
  other_opts <- list(...)
  solver_opts[names(other_opts)] <- other_opts

  results_dict <- ECOSolveR::ECOS_csolve(c = data[[C_KEY]], G = data[[G_KEY]], h = data[[H_KEY]], dims = data[[DIMS]], A = data[[A_KEY]], b = data[[B_KEY]], control = solver_opts)
  format_results(solver, results_dict, data, cached_data)
})

#' @param results_dict A list containing the solver output.
#' @param data A list containing information about the problem.
#' @describeIn ECOS Convert raw solver output into standard list of results.
setMethod("format_results", "ECOS", function(solver, results_dict, data, cached_data) {
  new_results <- list()
  status <- status_map(solver, results_dict$retcodes[["exitFlag"]])
  new_results[[STATUS]] <- status

  # Timing data
  new_results[[SOLVE_TIME]] <- results_dict$timing[["tsolve"]]
  new_results[[SETUP_TIME]] <- results_dict$timing[["tsetup"]]
  new_results[[NUM_ITERS]] <- results_dict$retcodes[["iter"]]

  if(new_results[[STATUS]] %in% SOLUTION_PRESENT) {
    primal_val <- results_dict$summary[["pcost"]]
    new_results[[VALUE]] <- primal_val + data[[OFFSET]]
    new_results[[PRIMAL]] <- results_dict$x
    new_results[[EQ_DUAL]] <- results_dict$y
    new_results[[INEQ_DUAL]] <- results_dict$z
  }
  new_results
})

#'
#' The ECOS_BB class.
#'
#' This class is an interface for the ECOS BB (branch-and-bound) solver.
#'
#' @references A. Domahidi, E. Chu, and S. Boyd. "ECOS: An SOCP solver for Embedded Systems." \emph{Proceedings of the European Control Conference}, pp. 3071-3076, 2013.
#' @seealso \code{\link[ECOSolveR]{ECOS_csolve}} and the \href{https://github.com/embotech/ecos}{ECOS Official Repository}.
#' @name ECOS_BB-class
#' @aliases ECOS_BB
#' @rdname ECOS_BB-class
setClass("ECOS_BB", contains = "ECOS")

#' @rdname ECOS_BB-class
#' @examples
#' ecos_bb <- ECOS_BB()
#' lp_capable(ecos_bb)
#' sdp_capable(ecos_bb)
#' socp_capable(ecos_bb)
#' exp_capable(ecos_bb)
#' mip_capable(ecos_bb)
#' @export
ECOS_BB <- function() {
  new("ECOS_BB")
  ##ECOS_BB$new()
}

# ECOS_BB capabilities
#' @param object,solver A \linkS4class{ECOS_BB} object.
#' @describeIn ECOS_BB ECOS_BB can handle linear programs.
setMethod("lp_capable", "ECOS_BB", function(solver) { TRUE })

#' @describeIn ECOS_BB ECOS_BB can handle second-order cone programs.
setMethod("socp_capable", "ECOS_BB", function(solver) { TRUE })

#' @describeIn ECOS_BB ECOS_BB cannot handle semidefinite programs.
setMethod("sdp_capable", "ECOS_BB", function(solver) { FALSE })

#' @describeIn ECOS_BB ECOS_BB cannot handle exponential cone programs.
setMethod("exp_capable", "ECOS_BB", function(solver) { FALSE })

#' @describeIn ECOS_BB ECOS_BB can handle mixed-integer programs.
setMethod("mip_capable", "ECOS_BB", function(solver) { TRUE })

# EXITCODES from ECOS_BB
# MI_OPTIMAL_SOLN (ECOS_OPTIMAL)                               ECOS_BB found optimal solution
# MI_INFEASIBLE (ECOS_PINF)                                    ECOS_BB proved problem is infeasible
# MI_UNBOUNDED (ECOS_DINF)                                     ECOS_BB proved problem is unbounded
# MI_MAXITER_FEASIBLE_SOLN (ECOS_OPTIMAL + ECOS_INACC_OFFSET)  ECOS_BB hit maximum iterations, but a feasible solution was found and the best seen feasible solution was returned
# MI_MAXITER_NO_SOLN (ECOS_PINF + ECOS_INACC_OFFSET)           ECOS_BB hit maximum iterations without finding a feasible solution
# MI_MAXITER_UNBOUNDED (ECOS_DINF + ECOS_INACC_OFFSET)         ECOS_BB hit maximum interations without finding a feasible solution that was unbounded

#' @describeIn ECOS_BB The name of the solver.
#' @export
setMethod("name", "ECOS_BB", function(object) { ECOS_BB_NAME })

#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @importFrom ECOSolveR ECOS_csolve ecos.control
#' @describeIn ECOS_BB Call the solver on the canonicalized problem.
setMethod("Solver.solve", "ECOS_BB", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)

  # TODO: Naras please fix this by making ECOSolveR handle type conversion, e.g. logical -> integer
  # TODO: Naras please fix this by making ECOSolveR handle type conversion, e.g. logical -> integer
  ## Naras Answer: Fixed in ECOSolveR version >= 0.4
  ## if(prod(dim(data[[G_KEY]])) == 0) data[[G_KEY]] <- NULL
  ## if(prod(dim(data[[A_KEY]])) == 0) data[[A_KEY]] <- NULL
  ## data[[DIMS]] <- lapply(data[[DIMS]], function(dim) { as.integer(dim) })
  storage.mode(data[[BOOL_IDX]]) <- "integer"
  storage.mode(data[[INT_IDX]]) <- "integer"

  solver_opts <- ECOSolveR::ecos.control()
  solver_opts$VERBOSE <- as.integer(verbose)
  other_opts <- list(...)
  solver_opts[names(other_opts)] <- other_opts

  results_dict <- ECOSolveR::ECOS_csolve(c = data[[C_KEY]], G = data[[G_KEY]], h = data[[H_KEY]], dims = data[[DIMS]], A = data[[A_KEY]], b = data[[B_KEY]],
                                         bool_vars = data[[BOOL_IDX]], int_vars = data[[INT_IDX]], control = solver_opts)
  format_results(solver, results_dict, data, cached_data)
})

#'
#' The SCS class.
#'
#' This class is an interface for the SCS solver.
#'
#' @references B. O'Donoghue, E. Chu, N. Parikh, and S. Boyd. "Conic Optimization via Operator Splitting and Homogeneous Self-Dual Embedding." \emph{Journal of Optimization Theory and Applications}, pp. 1-27, 2016. \url{https://doi.org/10.1007/s10957-016-0892-3}.
#' @seealso \code{\link[scs]{scs}} and the \href{https://github.com/cvxgrp/scs}{SCS Github}.
#' @name SCS-class
#' @aliases SCS
#' @rdname SCS-class
setClass("SCS", contains = "ECOS")

#' @rdname SCS-class
#' @examples
#' scs <- SCS()
#' lp_capable(scs)
#' sdp_capable(scs)
#' socp_capable(scs)
#' exp_capable(scs)
#' mip_capable(scs)
#' @export

#' @export
SCS <- function() {
    new("SCS")
    ##SCS$new()
}

# SCS capabilities
#' @param object,solver A \linkS4class{SCS} object.
#' @describeIn SCS SCS can handle linear programs.
setMethod("lp_capable", "SCS", function(solver) { TRUE })

#' @describeIn SCS SCS can handle second-order cone programs.
setMethod("socp_capable", "SCS", function(solver) { TRUE })

#' @describeIn SCS SCS can handle semidefinite programs.
setMethod("sdp_capable", "SCS", function(solver) { TRUE })

#' @describeIn SCS SCS can handle exponential cone programs.
setMethod("exp_capable", "SCS", function(solver) { TRUE })

#' @describeIn SCS SCS cannot handle mixed-integer programs.
setMethod("mip_capable", "SCS", function(solver) { FALSE })

#'
#' SCS Status Map
#'
#' Map of SCS status to CVXR status.
#'
#' @param solver A \linkS4class{SCS} object.
#' @param status An exit code returned by SCS.
#' @return A string indicating the status, either "optimal", "infeasible", "unbounded", "optimal_inaccurate", "infeasible_inaccurate", "unbounded_inaccurate", or "solver_error".
#' @docType methods
#' @rdname SCS-status_map
setMethod("status_map", "SCS", function(solver, status) {
  if(status == "Solved") OPTIMAL
  else if(status == "Solved/Inaccurate") OPTIMAL_INACCURATE
  else if(status == "Unbounded") UNBOUNDED
  else if(status == "Unbounded/Inaccurate") UNBOUNDED_INACCURATE
  else if(status == "Infeasible") INFEASIBLE
  else if(status == "Infeasible/Inaccurate") INFEASIBLE_INACCURATE
  else if(status %in% c("Failure", "Indeterminate")) SOLVER_ERROR
  else stop("SCS status unrecognized: ", status)
})

#' @describeIn SCS The name of the solver.
#' @export
setMethod("name", "SCS", function(object) { SCS_NAME })

#' @describeIn SCS Imports the scs library.
setMethod("import_solver", "SCS", function(solver) { TRUE }) ## Sincs scs is in imports!

setMethod("split_constr", "SCS", function(solver, constr_map) {
  list(eq_constr = c(constr_map[[EQ_MAP]], constr_map[[LEQ_MAP]]), ineq_constr = list(), nonlin_constr = list())
})

#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @importFrom scs scs
#' @describeIn SCS Call the solver on the canonicalized problem.
setMethod("Solver.solve", "SCS", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)

  # Set the options to be VERBOSE plus any user-specific options
  solver_opts <- list(...)
  solver_opts$verbose <- verbose
  ## Fix for acceleration_loopback parameter set to 20 by default in SCS
  if (is.null(solver_opts$acceleration_lookback)) solver_opts$acceleration_lookback  <- 10L

  scs_args <- list(A = data[[A_KEY]], b = data[[B_KEY]], obj = data[[C_KEY]], cone = data[[DIMS]])

  # If warm starting, add old primal and dual variables
  solver_cache <- cached_data[name(solver)]
  if(warm_start && !is.na(solver_cache$prev_result)) {
    stop("Warm start currently unimplemented")
    scs_args$x <- solver_cache$prev_result$x
    scs_args$y <- solver_cache$prev_result$y
    scs_args$s <- solver_cache$prev_result$s
  }

  results_dict <- do.call(scs::scs, c(scs_args, list(control = solver_opts)))
  # if(length(solver_opts) == 0)
  #  results_dict <- scs::scs(A = data[[A_KEY]], b = data[[B_KEY]], obj = data[[C_KEY]], cone = data[[DIMS]])
  # else
  #  results_dict <- scs::scs(A = data[[A_KEY]], b = data[[B_KEY]], obj = data[[C_KEY]], cone = data[[DIMS]], control = solver_opts)
  format_results(solver, results_dict, data, cached_data)
})

#' @param results_dict A list containing the solver output.
#' @param data A list containing information about the problem.
#' @describeIn SCS Convert raw solver output into standard list of results.
setMethod("format_results", "SCS", function(solver, results_dict, data, cached_data) {
  solver_cache <- cached_data[name(solver)]
  dims <- data[[DIMS]]
  new_results <- list()
  status <- status_map(solver, results_dict$info$status)
  new_results[[STATUS]] <- status

  # Timing and iteration data
  new_results[[SOLVE_TIME]] <- results_dict$info$solveTime/1000
  new_results[[SETUP_TIME]] <- results_dict$info$setupTime/1000
  new_results[[NUM_ITERS]] <- results_dict$info$iter

  if(new_results[[STATUS]] %in% SOLUTION_PRESENT) {
    # Save previous result for possible future warm start
    solver_cache$prev_result <- list(x = results_dict$x, y = results_dict$y, s = results_dict$s)
    primal_val <- results_dict$info$pobj
    new_results[[VALUE]] <- primal_val + data[[OFFSET]]
    new_results[[PRIMAL]] <- results_dict$x
    if(dims[[EQ_DIM]] > 0)
      new_results[[EQ_DUAL]] <- results_dict$y[1:dims[[EQ_DIM]]]
    else
      new_results[[EQ_DUAL]] <- numeric(0)

    y <- results_dict$y[(dims[[EQ_DIM]]+1):length(results_dict$y)]
    if(is.null(dims[[SDP_DIM]])) {
      old_sdp_sizes <- 0
      new_sdp_sizes <- 0
    } else {
      old_sdp_sizes <- sum(floor(dims[[SDP_DIM]] * (dims[[SDP_DIM]] + 1)/2))
      new_sdp_sizes <- sum(dims[[SDP_DIM]]^2)
    }
    if(is.null(dims[[SOC_DIM]]))
      y_offset <- dims[[LEQ_DIM]]
    else
      y_offset <- dims[[LEQ_DIM]] + sum(dims[[SOC_DIM]])
    y_true <- rep(0, length(y) + (new_sdp_sizes - old_sdp_sizes))
    y_true_offset <- y_offset
    if(y_true_offset > 0)
      y_true[1:y_true_offset] <- y[1:y_offset]

    # Expand SDP duals from lower triangular to full matrix, scaling off diagonal entries by 1/sqrt(2)
    for(n in dims[[SDP_DIM]]) {
      if(n > 0) {
        tri <- y[(y_offset + 1):(y_offset + floor(n*(n+1)/2))]
        y_true[(y_true_offset + 1):(y_true_offset + n^2)] <- SCS.tri_to_full(tri, n)
        y_true_offset <- y_true_offset + n^2
        y_offset <- y_offset + floor(n*(n+1)/2)
      }
    }

    if(length(y_true) > y_true_offset)
      y_true[(y_true_offset + 1):length(y_true)] <- y[(y_offset + 1):length(y)]
    new_results[[INEQ_DUAL]] <- y_true
  } else {
    # No result to save
    solver_cache$prev_result <- NULL
  }
  return(new_results)
})

# Expands floor(n*(n+1)/2) lower triangular entries to a full matrix with off-diagonal entries scaled by 1/sqrt(2).
SCS.tri_to_full <- function(lower_tri, n) {
  # Expands floor(n*(n+1)/2) lower triangular to full matrix, with off-diagonal entries scaled by 1/sqrt(2)
  full <- matrix(0, nrow = n, ncol = n)
  for(col in 1:n) {
    for(row in col:n) {
      idx <- row - col + floor(n*(n+1)/2) - floor((n-col+1)*(n-col+2)/2) + 1
      if(row != col) {
        full[row, col] <- lower_tri[idx]/sqrt(2)
        full[col, row] <- lower_tri[idx]/sqrt(2)
      } else
        full[row, col] <- lower_tri[idx]
    }
  }
  return(matrix(full, nrow = n^2))
}

# # TODO: This is a Python solver, which is partially ported to R via the cccp library.
# setClass("CVXOPT", contains = "Solver")
# CVXOPT <- function() {
#   stop("Unimplemented solver")
#   new("CVXOPT")
#   ##CVXOPT$new()
# }

## #'
## #' The LS class.
## #'
## #' This class represents a linearly constrained least squares solver using R's \code{base::solve} function.
## #' LS is incapable of solving any general cone program and must be invoked through a special path.
## #'
## #' @name LS-class
## #' @aliases LS
## #' @rdname LS-class
## setClass("LS", contains = "Solver")

## #' @rdname LS-class
## #' @export
## LS <- function() {
##   stop("Unimplemented solver")
##   new("LS")
## }

# # LS is incapable of solving any general cone program and must be invoked through a special path
# #' @param object,solver A \linkS4class{LS} object.
# #' @describeIn LS Returns \code{FALSE} since LS must be invoked through a special path.
# setMethod("lp_capable", "LS", function(solver) { FALSE })
#
# #' @describeIn LS Returns \code{FALSE} since LS must be invoked through a special path.
# setMethod("socp_capable", "LS", function(solver) { FALSE })
#
# #' @describeIn LS Returns \code{FALSE} since LS must be invoked through a special path.
# setMethod("sdp_capable", "LS", function(solver) { FALSE })
#
# #' @describeIn LS Returns \code{FALSE} since LS must be invoked through a special path.
# setMethod("exp_capable", "LS", function(solver) { FALSE })
#
# #' @describeIn LS Returns \code{FALSE} since LS must be invoked through a special path.
# setMethod("mip_capable", "LS", function(solver) { FALSE })
#
# #' @describeIn LS The name of the solver.
# #' @export
# setMethod("name", "LS", function(object) { LS_NAME })
#
# #' @describeIn LS Imports the Matrix library.
# setMethod("import_solver", "LS", function(solver) { requireNamespace("Matrix") })
#
# setMethod("split_constr", "LS", function(solver, constr_map) {
#   list(eq_constr = constr_map[[EQ_MAP]], ineq_constr = constr_map[[LEQ_MAP]], nonlin_constr = list())
# })
#
# #' @param objective A list representing the canonicalized objective.
# #' @param constraints A list of canonicalized constraints.
# #' @param cached_data A list mapping solver name to cached problem data.
# #' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
# #' @param verbose A logical value indicating whether to print solver output.
# #' @param ... Additional arguments to the solver.
# #' @describeIn LS Call the solver on the canonicalized problem.
# setMethod("Solver.solve", "LS", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
#   sym_data <- get_sym_data(solver, objective, constraints)
#   id_map <- sym_data@var_offsets
#   N <- sym_data@x_length
#   extractor <- QuadCoeffExtractor(id_map, N)   # TODO: QuadCoeffExtractor is unimplemented. See cvxpy/utilities/quadratic.py
#
#   # Extract the coefficients
#   coeffs <- get_coeffs(extractor, objective@args[[1]])
#   P <- coeffs$Ps[[1]]
#   q <- as.numeric(coeffs$Q)
#   r <- coeffs$R[[1]]
#
#   # Forming the KKT system
#   if(length(constraints) > 0) {
#     Cs <- lapply(constraints, function(c) {
#       coeffs <- get_coeffs(extractor, c@.expr)
#       if(length(coeffs) > 1)
#         coeffs[2:length(coeffs)]
#       else
#         c()
#     })
#     As <- do.call("rbind", lapply(Cs, function(C) { C[[1]] }))
#     bs <- as.numeric(sapply(Cs, function(C) { C[[2]] }))
#     lhs <- rbind(cbind(2*P, t(As), cbind(As, NA)))    # TODO: Fix this
#     rhs <- c(-q, -bs)
#   } else {   # Avoid calling rbind with empty list
#     lhs <- 2*P
#     rhs <- -q
#   }
#
#   # Actually solve the KKT system
#   tryCatch({
#       sol <- base::solve(lhs, rhs)
#       if(N > 0)
#         x <- sol[1:N]
#       else
#         x <- c()
#
#       if(length(sol) > N)
#         nu <- sol[(N+1):length(sol)]
#       else
#         nu <- c()
#
#       p_star <- t(x) %*% (P %*% x + q) + r
#     }, warning = function(w) {
#       x <- NA
#       nu <- NA
#       p_star <- NA
#     })
#
#   results_dict <- list()
#   results_dict[[PRIMAL]] <- x
#   results_dict[[EQ_DUAL]] <- nu
#   results_dict[[VALUE]] <- primal_to_result(objective, p_star)
#
#   format_results(solver, results_dict, NA, cached_data)
# })
#
# #' @param results_dict A list containing the solver output.
# #' @param data A list containing information about the problem.
# #' @describeIn LS Convert raw solver output into standard list of results.
# setMethod("format_results", "LS", function(solver, results_dict, data, cached_data) {
#   new_results <- results_dict
#   if(is.na(results_dict[[PRIMAL]]))
#     new_results[[STATUS]] <- INFEASIBLE
#   else
#     new_results[[STATUS]] <- OPTIMAL
#   new_results
# })


