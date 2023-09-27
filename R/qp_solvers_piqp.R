
#'
#' QP interface for the PIQP solver.
#'
#' @name PIQP-class
#' @aliases PIQP
#' @rdname PIQP-class
#' @export
setClass("PIQP", contains = "QpSolver")

#' @rdname PIQP-class
#' @export
PIQP <- function() { new("PIQP") }

#' @param solver,object,x A \linkS4class{PIQP} object.
#' @param status A status code returned by the solver.
#' @param default A status string to return if no status code match is found. If \code{default = NA}, this method will return an error when there is no match.
#' @describeIn PIQP Converts status returned by the PIQP solver to its respective CVXR status.
setMethod("status_map", "PIQP", function(solver, status, default = NA_character_) {
  PIQP_STATUS_MAP <- list(
    "1" = OPTIMAL,
    "-1" = USER_LIMIT, # Maxiter reached
    "-2" = INFEASIBLE, # primal infeasible
    "-3" = UNBOUNDED,  # dual infeasible
    "-8" = SOLVER_ERROR, # numerical error during solving
    "-10" = SOLVER_ERROR # invalid settings provided to solver
  )

  status_string <- PIQP_STATUS_MAP[[as.character(status)]]
  if (is.null(status_string)) {
    if (is.na(default))
      stop("PIQP status unrecognized: ", status)
    else
      default
  } else {
    status_string
  }
})

#' @describeIn PIQP Returns the name of the solver.
setMethod("name", "PIQP", function(x) { PIQP_NAME })

#' @describeIn PIQP Imports the solver.
setMethod("import_solver", "PIQP", function(solver) { requireNamespace("piqp", quietly = TRUE) })

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn PIQP Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "PIQP", solution = "list", inverse_data = "InverseData"), function(object, solution, inverse_data) {

  attr <- list(SOLVE_TIME = solution$info$run_time,
               EXTRA_STATS =  solution)

  # Map PIQP statuses back to CVXR statuses.
  status <- status_map(object, solution$info$status_val, default = SOLVER_ERROR)

  if (status %in% SOLUTION_PRESENT) {
    opt_val <- solution$info$obj_val + inverse_data[[OFFSET]]
    primal_vars <- list(object@VAR_ID = as.matrix(solution$x))
    dual_vars <- list(object@DUAL_VAR_ID = solution$y)
    attr[[NUM_ITERS]] <- solution$info$iter
    sol <- Solution(status, opt_val, primal_vars, dual_vars, attr)
  } else
    sol <- failure_solution(status)
  return(sol)
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @importFrom piqp piqp
#' @describeIn PIQP Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "PIQP", function(object, data, warm_start, verbose, solver_opts, solver_cache = NULL) {
  if(is.null(solver_cache)) solver_cache  <- new.env(parent = emptyenv())

  P <- data[[P_KEY]]
  A <- data[[A_KEY]]
  F <- data[[F_KEY]]
  q <- data[[Q_KEY]]
  b <- data[[B_KEY]]
  g <- data[[G_KEY]]

  ## Ensure our solver defaults and reconcile user overrides/preferences
  solver_opts <- reconcile_solver_options(solver_opts, SOLVER_DEFAULT_PARAM$PIQP)

  # Caching is not implemented it seems, but piqp does allow model updates. TODO later
  # Initialize and solve problem.
  solver <- piqp::piqp(P = P, c = q, A = A, b = b, G = F, h = g,
                       settings = c(list(verbose = verbose), solver_opts))
  results <- solver$solve()

  if(!is.null(solver_cache)) {
    solver_cache[[name(object)]] <- list(solver, data, results)
  }

  return(list(results, solver_cache))
})
