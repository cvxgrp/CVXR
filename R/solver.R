#'
#' The Solver class.
#'
#' This virtual class represents the generic interface for a solver.
#' 
#' @name Solver-class
#' @aliases Solver
#' @rdname Solver-class
Solver <- setClass("Solver", contains = "VIRTUAL")

#
# Choose a Solver
#
# Determines the appropriate solver.
# 
# @param constraints A list of canonicalized constraints.
# @return A \linkS4class{Solver} object.
# @rdname Solver-choose_solver
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
    if("cvxopt" %in% installed.packages()) {
      require(cvxopt)
      return(CVXOPT())
    } else
      return(SCS())
  }
  # Otherwise use ECOS
  else
    return(ECOS())
}

#' @describeIn Solver Raises an exception if the solver cannot solver the problem.
#' @param solver A \linkS4class{Solver} object.
#' @param constraints A list of canonicalized constraints
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

#' @describeIn Solver Clears the cache if the objective or constraints changed.
#' @param solver A \linkS4class{Solver} object.
#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @return The updated \code{cached_data}.
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

#' @describeIn Solver Returns the symbolic data for the problem.
#' @param solver A \linkS4class{Solver} object.
#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @return A \linkS4class{SymData} object holding the symbolic data for the problem.
setMethod("get_sym_data", "Solver", function(solver, objective, constraints, cached_data) {
  cached_data <- validate_cache(solver, objective, constraints, cached_data)
  prob_data <- cached_data[[name(solver)]]
  if(is.null(prob_data@sym_data))
    prob_data@sym_data <- SymData(objective, constraints, solver)
  cached_data[[name(solver)]] <- prob_data
  cached_data
})

#' @describeIn Solver Returns the numeric data for the problem.
#' @param solver A \linkS4class{Solver} object.
#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @return A \linkS4class{SymData} object holding the symbolic data for the problem.
setMethod("get_matrix_data", "Solver", function(solver, objective, constraints, cached_data) {
  cached_data <- get_sym_data(solver, objective, constraints, cached_data)
  sym_data <- cached_data[[name(solver)]]@sym_data
  prob_data <- cached_data[[name(solver)]]
  if(is.null(prob_data@matrix_data))
    prob_data@matrix_data <- MatrixData(sym_data, solver)   # TODO: Update this constructor for nonlinear constraints
  cached_data[[name(solver)]] <- prob_data
  cached_data
})

#' @describeIn Solver Returns the argument for the call to the solver.
#' @param solver A \linkS4class{Solver} object.
#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @return A list of the arguments needed for the solver.
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
  data[[BOOL_IDX]] <- conv_idx$bool_idx
  data[[INT_IDX]] <- conv_idx$int_idx
  data
})

#' @describeIn Solver A logical value indicating whether nonlinear constraints are needed.
setMethod("nonlin_constr", "Solver", function(solver) { FALSE })

#'
#' Call to Solver
#' 
#' Returns the result of the call to the solver.
#' 
#' @param solver A \linkS4class{Solver} object.
#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @return A list containing the status, optimal value, primal variable, and dual variables for the equality and inequality constraints.
#' @docType methods
#' @rdname Solver-solve
setMethod("Solver.solve", "Solver", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) { stop("Unimplemented") })

#' 
#' Format Solver Results
#' 
#' Converts the solver output into standard form.
#' 
#' @param solver A \linkS4class{Solver} object.
#' @param results_dict A list containing the solver output.
#' @param data A list containing information about the problem.
#' @param cached_data A list mapping solver name to cached problem data.
#' @return A list containing the solver output in standard form.
#' @docType methods
#' @rdname format_results
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
#' @references A. Domahidi, E. Chu, and S. Boyd. "ECOS: An SOCP solver for Embedded Systems." \emph{Proceedings of the European Control Conference}, pp. 3071-3076, 2013. \url{http://web.stanford.edu/~boyd/papers/ecos.html}.
#' @seealso \code{\link[ECOSolveR]{ECOS_csolve}} and the \href{https://www.embotech.com/ECOS}{ECOS Official Site}.
#' @name ECOS-class
#' @rdname ECOS-class
#' @export
ECOS <- setClass("ECOS", contains = "Solver")

#' @name ECOS
#' @rdname ECOS-class
#' @export
ECOS <- function() {
    new("ECOS")
    ##ECOS$new()
}

# ECOS capabilities
#' @rdname Solver-capable
setMethod("lp_capable", "ECOS", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("socp_capable", "ECOS", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("sdp_capable", "ECOS", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("exp_capable", "ECOS", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("mip_capable", "ECOS", function(solver) { FALSE })

#' 
#' ECOS Status Map
#' 
#' Map of ECOS status to CVXR status.
#'
#' @param solver A \linkS4class{ECOS} solver object.
#' @param status An exit code returned by ECOS:
#' \itemize{
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
setMethod("name", "ECOS", function(object) { ECOS_NAME })

#' @describeIn ECOS Imports the ECOSolveR library.
setMethod("import_solver", "ECOS", function(solver) { requireNamespace("ECOSolveR") })

#' @describeIn ECOS The interface for matrices passed to the solver.
setMethod("matrix_intf", "ECOS", function(solver) { DEFAULT_SPARSE_INTF })

#' @describeIn ECOS The interface for vectors passed to the solver.
setMethod("vec_intf", "ECOS", function(solver) { DEFAULT_INTF })

#' @describeIn ECOS Extracts the equality, inequality, and nonlinear constraints.
#' @param solver A \linkS4class{ECOS} object.
#' @param constr_map A list of canonicalized constraints.
#' @return A list of equality, inequality, and nonlinear constraints.
setMethod("split_constr", "ECOS", function(solver, constr_map) {
  list(eq_constr = constr_map[[EQ_MAP]], ineq_constr = constr_map[[LEQ_MAP]], nonlin_constr = list())
})

#' @rdname Solver-solve
setMethod("Solver.solve", "ECOS", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)
  data[[DIMS]]['e'] <- data[[DIMS]][[EXP_DIM]]

  # TODO: Naras please fix this by making ECOSolveR handle type conversion, e.g. logical -> integer
  if(prod(dim(data[[G_KEY]])) == 0) data[[G_KEY]] <- NULL
  if(prod(dim(data[[A_KEY]])) == 0) data[[A_KEY]] <- NULL
  data[[DIMS]] <- lapply(data[[DIMS]], function(dim) { as.integer(dim) })
  solver_opts <- ECOSolveR::ecos.control()
  solver_opts$VERBOSE <- as.integer(verbose)
  other_opts <- list(...)
  solver_opts[names(other_opts)] <- other_opts

  results_dict <- ECOSolveR::ECOS_csolve(c = data[[C_KEY]], G = data[[G_KEY]], h = data[[H_KEY]], dims = data[[DIMS]], A = data[[A_KEY]], b = data[[B_KEY]], control = solver_opts)
  format_results(solver, results_dict, data, cached_data)
})

#' @rdname format_results
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
#' @references A. Domahidi, E. Chu, and S. Boyd. "ECOS: An SOCP solver for Embedded Systems." \emph{Proceedings of the European Control Conference}, pp. 3071-3076, 2013. \url{http://web.stanford.edu/~boyd/papers/ecos.html}.
#' @seealso \code{\link[ECOSolveR]{ECOS_csolve}} and the \href{https://www.embotech.com/ECOS}{ECOS Official Site}.
#' @name ECOS_BB-class
#' @rdname ECOS_BB-class
#' @export
setClass("ECOS_BB", contains = "ECOS")

#' @name ECOS_BB
#' @rdname ECOS_BB-class
#' @export
ECOS_BB <- function() {
  new("ECOS_BB")
  ##ECOS_BB$new()
}

# ECOS_BB capabilities
#' @rdname Solver-capable
setMethod("lp_capable", "ECOS_BB", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("socp_capable", "ECOS_BB", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("sdp_capable", "ECOS_BB", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("exp_capable", "ECOS_BB", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("mip_capable", "ECOS_BB", function(solver) { TRUE })

# EXITCODES from ECOS_BB
# MI_OPTIMAL_SOLN (ECOS_OPTIMAL)                               ECOS_BB found optimal solution
# MI_INFEASIBLE (ECOS_PINF)                                    ECOS_BB proved problem is infeasible
# MI_UNBOUNDED (ECOS_DINF)                                     ECOS_BB proved problem is unbounded
# MI_MAXITER_FEASIBLE_SOLN (ECOS_OPTIMAL + ECOS_INACC_OFFSET)  ECOS_BB hit maximum iterations, but a feasible solution was found and the best seen feasible solution was returned
# MI_MAXITER_NO_SOLN (ECOS_PINF + ECOS_INACC_OFFSET)           ECOS_BB hit maximum iterations without finding a feasible solution
# MI_MAXITER_UNBOUNDED (ECOS_DINF + ECOS_INACC_OFFSET)         ECOS_BB hit maximum interations without finding a feasible solution that was unbounded

#' @describeIn ECOS_BB The name of the solver.
setMethod("name", "ECOS_BB", function(object) { ECOS_BB_NAME })

#' @rdname Solver-solve
setMethod("Solver.solve", "ECOS_BB", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)

  # TODO: Naras please fix this by making ECOSolveR handle type conversion, e.g. logical -> integer
  if(prod(dim(data[[G_KEY]])) == 0) data[[G_KEY]] <- NULL
  if(prod(dim(data[[A_KEY]])) == 0) data[[A_KEY]] <- NULL
  data[[DIMS]] <- lapply(data[[DIMS]], function(dim) { as.integer(dim) })
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
#' @references B. O'Donoghue, E. Chu, N. Parikh, and S. Boyd. "Conic Optimization via Operator Splitting and Homogeneous Self-Dual Embedding." \emph{Journal of Optimization Theory and Applications}, pp. 1-27, 2016. \url{https://web.stanford.edu/~boyd/papers/scs.html}.
#' @seealso \code{\link[scs]{scs}} and the \href{https://github.com/cvxgrp/scs}{SCS Github}.
#' @name SCS-class
#' @rdname SCS-class
#' @export
setClass("SCS", contains = "ECOS")

#' @name SCS
#' @rdname SCS-class
#' @export
SCS <- function() {
    new("SCS")
    ##SCS$new()
}

# SCS capabilities
#' @rdname Solver-capable
setMethod("lp_capable", "SCS", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("socp_capable", "SCS", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("sdp_capable", "SCS", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("exp_capable", "SCS", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("mip_capable", "SCS", function(solver) { FALSE })

#' 
#' SCS Status Map
#' 
#' Map of SCS status to CVXR status.
#'
#' @param solver A \linkS4class{SCS} solver object.
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
setMethod("name", "SCS", function(object) { SCS_NAME })

#' @describeIn SCS Imports the scs library.
setMethod("import_solver", "SCS", function(solver) { requireNamespace("scs") })

#' @describeIn SCS Extracts the equality, inequality, and nonlinear constraints.
#' @param solver A \linkS4class{SCS} object.
#' @param constr_map A list of canonicalized constraints.
#' @return A list of equality, inequality, and nonlinear constraints.
setMethod("split_constr", "SCS", function(solver, constr_map) {
  list(eq_constr = c(constr_map[[EQ_MAP]], constr_map[[LEQ_MAP]]), ineq_constr = list(), nonlin_constr = list())
})

#' @rdname Solver-solve
setMethod("Solver.solve", "SCS", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)

  # Set the options to be VERBOSE plus any user-specific options
  solver_opts <- list(...)
  solver_opts$verbose <- verbose
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

#' @rdname format_results
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

# TODO: This is a Python solver, which is partially ported to R via the cccp library.
setClass("CVXOPT", contains = "Solver")
CVXOPT <- function() {
  stop("Unimplemented solver")
  new("CVXOPT")
  ##CVXOPT$new()
}

#'
#' The LS class.
#' 
#' This class represents a linearly constrained least squares solver using R's \code{base::solve} function.
#' LS is capable of solving any general cone program and must be invoked through a special path.
#'
#' @name LS-class
#' @rdname LS-class
#' @export
setClass("LS", contains = "Solver")

#' @name LS
#' @rdname LS-class
#' @export
LS <- function() {
  stop("Unimplemented solver")
  new("LS")
}

# LS is incapable of solving any general cone program and must be invoked through a special path
#' @rdname Solver-capable
setMethod("lp_capable", "LS", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("socp_capable", "LS", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("sdp_capable", "LS", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("exp_capable", "LS", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("mip_capable", "LS", function(solver) { FALSE })

#' @describeIn LS The name of the solver.
setMethod("name", "LS", function(object) { LS_NAME })

#' @describeIn LS Imports the Matrix library.
setMethod("import_solver", "LS", function(solver) { requireNamespace("Matrix") })

#' @describeIn LS Extracts the equality, inequality, and nonlinear constraints.
#' @param solver A \linkS4class{LS} object.
#' @param constr_map A list of canonicalized constraints.
#' @return A list of equality, inequality, and nonlinear constraints.
setMethod("split_constr", "LS", function(solver, constr_map) {
  list(eq_constr = constr_map[[EQ_MAP]], ineq_constr = constr_map[[LEQ_MAP]], nonlin_constr = list())
})

#' @rdname Solver-solve
setMethod("Solver.solve", "LS", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  sym_data <- get_sym_data(solver, objective, constraints)
  id_map <- sym_data@var_offsets
  N <- sym_data@x_length
  extractor <- QuadCoeffExtractor(id_map, N)   # TODO: QuadCoeffExtractor is unimplemented. See cvxpy/utilities/quadratic.py
  
  # Extract the coefficients
  coeffs <- get_coeffs(extractor, objective@args[[1]])
  P <- coeffs$Ps[[1]]
  q <- as.numeric(coeffs$Q)
  r <- coeffs$R[[1]]
  
  # Forming the KKT system
  if(length(constraints) > 0) {
    Cs <- lapply(constraints, function(c) { 
      coeffs <- get_coeffs(extractor, c@.expr)
      if(length(coeffs) > 1)
        coeffs[2:length(coeffs)]
      else
        c()
    })
    As <- do.call("rbind", lapply(Cs, function(C) { C[[1]] }))
    bs <- as.numeric(sapply(Cs, function(C) { C[[2]] }))
    lhs <- rbind(cbind(2*P, t(As), cbind(As, NA)))    # TODO: Fix this
    rhs <- c(-q, -bs)
  } else {   # Avoid calling rbind with empty list
    lhs <- 2*P
    rhs <- -q
  }
  
  # Actually solve the KKT system
  tryCatch({
      sol <- solve(lhs, rhs)
      if(N > 0)
        x <- sol[1:N]
      else
        x <- c()
      
      if(length(sol) > N)
        nu <- sol[(N+1):length(sol)]
      else
        nu <- c()
      
      p_star <- t(x) %*% (P %*% x + q) + r
    }, warning = function(w) {
      x <- NA
      nu <- NA
      p_star <- NA
    })
  
  results_dict <- list()
  results_dict[[PRIMAL]] <- x
  results_dict[[EQ_DUAL]] <- nu
  results_dict[[VALUE]] <- primal_to_result(objective, p_star)
  
  format_results(solver, results_dict, NA, cached_data)
})

#' @rdname format_results
setMethod("format_results", "LS", function(solver, results_dict, data, cached_data) {
  new_results <- results_dict
  if(is.na(results_dict[[PRIMAL]]))
    new_results[[STATUS]] <- INFEASIBLE
  else
    new_results[[STATUS]] <- OPTIMAL
  new_results
})

#'
#' The MOSEK class.
#' 
#' This class is an interface for the commercial MOSEK solver.
#'
#' @references E. Andersen and K. Andersen. "The MOSEK Interior Point Optimizer for Linear Programming: an Implementation of the Homogeneous Algorithm." \emph{High Performance Optimization}, vol. 33, pp. 197-232, 2000.
#' @seealso \code{\link{Rmosek}{mosek}} and the \href{https://www.mosek.com/products/mosek/}{MOSEK Official Site}.
#' @name MOSEK-class
#' @rdname MOSEK-class
#' @export
setClass("MOSEK", contains = "Solver")

#' @name MOSEK
#' @rdname MOSEK-class
#' @export
MOSEK <- function() {
  stop("Unimplemented")
  new("MOSEK") 
}

#' @rdname Solver-capable
setMethod("lp_capable", "MOSEK", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("socp_capable", "MOSEK", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("sdp_capable", "MOSEK", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("exp_capable", "MOSEK", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("mip_capable", "MOSEK", function(solver) { FALSE })

#' 
#' MOSEK Status Map
#' 
#' Map of MOSEK status to CVXR status.
#'
#' @param solver A \linkS4class{MOSEK} solver object.
#' @param status An exit code returned by MOSEK. See the \href{http://docs.mosek.com/8.0/dotnetfusion/solution_status.html}{MOSEK documentation} for details.
#' @return A string indicating the status, either "optimal", "infeasible", "unbounded", "optimal_inaccurate", "infeasible_inaccurate", "unbounded_inaccurate", or "solver_error".
#' @docType methods
#' @rdname MOSEK-status_map
setMethod("status_map", "MOSEK", function(solver, status) {
  if(status == "OPTIMAL")
    OPTIMAL
  else if(status == "PRIMAL_INFEASIBLE_CER")
    INFEASIBLE
  else if(status == "DUAL_INFEASIBLE_CER")
    UNBOUNDED
  else if(status == "NEAR_OPTIMAL")
    OPTIMAL_INACCURATE
  else if(status == "NEAR_PRIMAL_INFEASIBLE_CER")
    INFEASIBLE_INACCURATE
  else if(status == "NEAR_DUAL_INFEASIBLE_CER")
    UNBOUNDED_INACCURATE
  else
    SOLVER_ERROR
})

#' @describeIn MOSEK The name of the solver.
setMethod("name", "MOSEK", function(object) { MOSEK_NAME })

#' @describeIn MOSEK Imports the Rmosek library.
setMethod("import_solver", "MOSEK", function(solver) { requireNamespace("Rmosek") })

#' @describeIn MOSEK Extracts the equality, inequality, and nonlinear constraints.
#' @param solver A \linkS4class{MOSEK} object.
#' @param constr_map A list of canonicalized constraints.
#' @return A list of equality, inequality, and nonlinear constraints.
setMethod("split_constr", "MOSEK", function(solver, constr_map) {
  list(eq_constr = constr_map[[EQ_MAP]], ineq_constr = constr_map[[LEQ_MAP]], nonlin_constr = list())
})

#' @rdname Solver-solve
setMethod("Solver.solve", "MOSEK", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)
  
  A <- data[[A_KEY]]
  b <- data[[B_KEY]]
  G <- data[[G_KEY]]
  h <- data[[H_KEY]]
  c <- data[[C_KEY]]
  dims <- data[[DIMS]]
  problem <- list(sense = "minimize")
  
  # Size of problem
  numvar <- length(c) + sum(dims[[SOC_DIM]])
  numcon <- length(b) + dims[[LEQ_DIM]] + sum(dims[[SOC_DIM]]) + sum(dims[[SDP_DIM]]^2)
  
  # TODO: Fix crash on empty problem
  
  # Objective
  problem$c <- c      # Objective coefficients
  problem$c0 <- 0     # Objective constant
  problem$bx <- rbind(blx = rep(-Inf, numvar), bux = rep(Inf, numvar))   # Lower and upper variable bounds
  
  # SDP variables
  if(sum(dims[[SDP_DIM]]) > 0)
    problem$bardim <- dims[[SDP_DIM]]   # Semidefinite variable dimensions
  
  # Linear equality and linear inequality constraints
  if(nrow(A) > 0 && nrow(G) > 0)
    constraints_matrix <- rbind(A, G)
  else if(nrow(A) > 0)
    constraints_matrix <- A
  else
    constraints_matrix <- G
  if(dims[[LEQ_DIM]] == length(h))
    problem$bc <- rbind(blc = c(b, rep(-Inf, dims[[LEQ_DIM]])), buc = c(b, h))
  else
    problem$bc <- rbind(blc = c(b, rep(-Inf, dims[[LEQ_DIM]]), h[(1 + dims[[LEQ_DIM]]):length(h)]), buc = c(b, h))
  # h_leq <- h[1:dims[[LEQ_DIM]]]
  # sdp_total_dims <- sum(dims[[SDP_DIM]]^2)
  # soc_sdp_dims <- sum(dims[[SOC_DIM]]) + sdp_total_dims
  # h_soc_sdp <- h[(1 + dims[[LEQ_DIM]]):(1 + dims[[LEQ_DIM]] + soc_sdp_dims)]
  # problem$bc <- rbind(blc = c(b, rep(-Inf, dims[[LEQ_DIM]]), h_soc_sdp), buc = c(b, h_leq, h_soc_sdp))
  
  # Cone constraints
  num_cones <- length(dims[[SOC_DIM]])
  if(num_cones > 0) {
    cur_var_idx <- length(c)
    cur_con_idx <- length(b) + dims[[LEQ_DIM]]
    cones <- matrix(list(), nrow = 2, ncol = num_cones)
    
    for(k in 1:num_cones) {
      size_cone <- dims[[SOC_DIM]][k]
      
      # Add an identity for each cone
      id_mat <- sparseMatrix(i = cur_con_idx + 1:size_cone, j = cur_var_idx + 1:size_cone, x = 1)
      constraints_matrix <- rbind(constraints_matrix, id_mat)
      
      # Add a cone constraint
      cones[,k] <- list("QUAD", seq(cur_var_idx + 1, cur_var_idx + size_cone))
      cur_var_idx <- cur_var_idx + size_cone
      cur_con_idx <- cur_con_idx + size_cone
    }
    rownames(cones) <- c("type", "sub")
    problem$cones <- cones
  }
  problem$A <- constraints_matrix

  # SDP constraints
  num_sdp <- length(dims[[SDP_DIM]])
  if(num_sdp > 0) {
    for(k in 1:num_sdp) {
      size_matrix <- dims[[SDP_DIM]][k]
      for(i in 1:length(size_matrix)) {
        for(j in 1:length(size_matrix)) {
          if(i == j)
            coeff <- 1
          else
            coeff <- 0.5
          # TODO: Finish this
          cur_con_idx <- cur_con_idx + 1
        }
      }
    }
  }
      
  requireNamespace("Rmosek")  
  results_dict <- rmosek(problem, opts = list(getinfo = TRUE, soldetail = TRUE, verbose = verbose, ...))
  format_results(solver, results_dict, data, cached_data)
})

#' @describeIn MOSEK Chooses between the basic and interior point solution.
#' @param solver A \linkS4class{MOSEK} object.
#' @param results_dict A list of the results returned by the solver.
#' @return A list containing the preferred solution (\code{solist}) and status of the preferred solution (\code{solsta}).
setMethod("choose_solution", "MOSEK", function(solver, results_dict) {
  requireNamespace("Rmosek")
  rank <- function(status) {
    # Rank solutions: optimal > near_optimal > anything else > None
    if(status == "OPTIMAL")
      return(3)
    else if(status == "NEAR_OPTIMAL")
      return(2)
    else if(!is.null(status))
      return(1)
    else
      return(0)
  }
  
  # As long as interior solution is not worse, take it (for backward compatibility)
  solsta_bas <- results_dict$bas$solsta
  solsta_itr <- results_dict$itr$solsta
  if(rank(solsta_itr) >= rank(solsta_bas))
    list(solist = results_dict$sol$itr, solsta = solsta_itr)
  else
    list(solist = results_dict$sol$bas, solsta = solsta_bas)
})

#' @rdname format_results
setMethod("format_results", "MOSEK", function(solver, results_dict, data, cached_data) {
  requireNamespace("Rmosek")
  sol <- choose_solution(solver, results_dict)
  
  new_results <- list()
  new_results[[STATUS]] <- status_map(sol$solsta)
  new_results[[SOLVE_TIME]] <- results_dict$dinfo$OPTIMIZER_TIME
  new_results[[SETUP_TIME]] <- results_dict$dinfo$PRESOLVE_TIME
  new_results[[NUM_ITERS]] <- results_dict$iinfo$INTPNT_ITER
  
  if(new_results[[STATUS]] %in% SOLUTION_PRESENT) {
    # Get primal variable values
    new_results[[PRIMAL]] <- sol$solist$xx

    # Get objective value
    new_results[[VALUE]] <- sol$solist$pobjval + data[[OFFSET]]
    
    # TODO: Check if signs are inverted in MOSEK
    y <- sol$solist$slc - sol$solist$suc
    new_results[[EQ_DUAL]] <- y[1:length(data[[B_KEY]])]
    if(length(data[[B_KEY]]) < length(y))
      new_results[[INEQ_DUAL]] <- y[(length(data[[B_KEY]]) + 1):(length(data[[B_KEY]]) + data[[DIMS]][[LEQ_DIM]])]
  }
  new_results
})

#'
#' The GUROBI class.
#'
#' This class is an interface for the commercial GUROBI solver.
#' 
#' @references \emph{Gurobi optimizer reference manual version 5.0,} Gurobi Optimization, Inc., Houston, Texas, July 2012.
#' @seealso \href{http://www.gurobi.com/documentation/7.5/refman/r_api_overview.html}{GUROBI Official Site}.
#' @name GUROBI-class
#' @rdname GUROBI-class
#' @export
setClass("GUROBI", contains = "Solver")

#' @name GUROBI
#' @rdname GUROBI-class
#' @export
GUROBI <- function() {
  stop("Unimplemented")
  new("GUROBI")
}

#' @rdname Solver-capable
setMethod("lp_capable", "GUROBI", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("socp_capable", "GUROBI", function(solver) { TRUE })

#' @rdname Solver-capable
setMethod("sdp_capable", "GUROBI", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("exp_capable", "GUROBI", function(solver) { FALSE })

#' @rdname Solver-capable
setMethod("mip_capable", "GUROBI", function(solver) { TRUE })

#' 
#' GUROBI Status Map
#' 
#' Map of GUROBI status to CVXR status.
#'
#' @param solver A \linkS4class{GUROBI} solver object.
#' @param status An exit code returned by GUROBI. See the \href{http://www.gurobi.com/documentation/7.5/refman/optimization_status_codes.html}{GUROBI documentation} for details.
#' @return A string indicating the status, either "optimal", "infeasible", "unbounded", "optimal_inaccurate", "infeasible_inaccurate", "unbounded_inaccurate", or "solver_error".
#' @docType methods
#' @rdname GUROBI-status_map
setMethod("status_map", "GUROBI", function(solver, status) {
  if(status == 2)
    OPTIMAL
  else if(status == 3)
    INFEASIBLE
  else if(status == 5)
    UNBOUNDED
  else if(status == 9)
    OPTIMAL_INACCURATE
  else if(status %in% c(4, 6, 7, 8, 10, 11, 12, 13))
    SOLVER_ERROR
  else
    stop("GUROBI status unrecognized: ", status)
})

#' @describeIn GUROBI The name of the solver.
setMethod("name", "GUROBI", function(object) { GUROBI_NAME })

#' @describeIn GUROBI Imports the gurobi library.
setMethod("import_solver", "GUROBI", function(solver) { requireNamespace("gurobi") })

#' @describeIn GUROBI Extracts the equality, inequality, and nonlinear constraints.
#' @param solver A \linkS4class{GUROBI} object.
#' @param constr_map A list of canonicalized constraints.
#' @return A list of equality, inequality, and nonlinear constraints.
setMethod("split_constr", "GUROBI", function(solver, constr_map) {
  list(eq_constr = c(constr_map[[EQ_MAP]], constr_map[[LEQ_MAP]]), ineq_constr = list(), nonlin_constr = list())
})

#' @rdname format_results
setMethod("format_results", "GUROBI", function(solver, results_dict, data, cached_data) {
  dims <- data[[DIMS]]
  if(results_dict$status != SOLVER_ERROR) {
    solver_cache <- cached_data[[name(solver)]]
    solver_cache@prev_result <- list(vbasis = results_dict$vbasis, cbasis = results_dict$cbasis, c = data[[C_KEY]], A = data[[A_KEY]], b = data[[B_KEY]])
  }
  
  new_results <- list()
  new_results[[STATUS]] <- results_dict$status
  new_results[[SOLVE_TIME]] <- results_dict$runtime
  if(new_results[[STATUS]] %in% SOLUTION_PRESENT) {
    primal_val <- results_dict$objval
    new_results[[VALUE]] <- primal_val + data[[OFFSET]]
    new_results[[PRIMAL]] <- results_dict$x
    if(!Solver.is_mip(data)) {
      duals <- results_dict$pi
      if(dims[[EQ_DIM]] > 0)
        new_results[[EQ_DUAL]] <- duals[1:dims[[EQ_DIM]]]
      if(length(duals) > dims[[EQ_DIM]])
        new_results[[INEQ_DUAL]] <- duals[(dims[[EQ_DIM]] + 1):length(duals)]
    }
  }
})

####################
#                  #
# Solver utilities #
#                  #
####################
# solver_intf <- list(ECOS(), ECOS_BB(), CVXOPT(), GLPK(), GLPK_MI(), CBC(), SCS(), GUROBI(), Elemental(), MOSEK(), LS())
solver_intf <- list(ECOS(), ECOS_BB(), SCS())
SOLVERS <- solver_intf
names(SOLVERS) <- sapply(solver_intf, function(solver) { name(solver) })

installed_solvers <- function() {
  installed <- list()
  for(i in 1:length(SOLVERS)) {
    if(is_installed(SOLVERS[i]))
      installed <- c(installed, names(SOLVERS)[i])
  }
  installed
}
