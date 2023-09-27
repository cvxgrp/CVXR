
#'
#' QP interface for the OSQP solver.
#'
#' @name OSQP-class
#' @aliases OSQP
#' @rdname OSQP-class
#' @export
setClass("OSQP", contains = "QpSolver")

#' @rdname OSQP-class
#' @export
OSQP <- function() { new("OSQP") }

#' @param solver,object,x A \linkS4class{OSQP} object.
#' @param status A status code returned by the solver.
#' @param default A status string to return if no status code match is found. If \code{default = NA}, this method will return an error when there is no match.
#' @describeIn OSQP Converts status returned by the OSQP solver to its respective CVXR status.
setMethod("status_map", "OSQP", function(solver, status, default = NA_character_) {
  OSQP_STATUS_MAP <- list(
    "1" = OPTIMAL,
    "2" = OPTIMAL_INACCURATE,
    "-2" = SOLVER_ERROR,           # Maxiter reached
    "-3" = INFEASIBLE,
    "3" = INFEASIBLE_INACCURATE,
    "-4" = UNBOUNDED,
    "4" = UNBOUNDED_INACCURATE,
    "-6" = USER_LIMIT,
    "-5" = SOLVER_ERROR,           # Interrupted by user
    "-10" = SOLVER_ERROR
  )

  status_string <- OSQP_STATUS_MAP[[as.character(status)]]
  if (is.null(status_string)) {
    if (is.na(default))
      stop("OSQP status unrecognized: ", status)
    else
      default
  } else {
    status_string
  }
})

#' @describeIn OSQP Returns the name of the solver.
setMethod("name", "OSQP", function(x) { OSQP_NAME })

#' @describeIn OSQP Imports the solver.
## setMethod("import_solver", "OSQP", function(solver) { requireNamespace("osqp", quietly = TRUE) })
## Since OSQP is a requirement, this is always TRUE
setMethod("import_solver", "OSQP", function(solver) { TRUE })

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn OSQP Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "OSQP", solution = "list", inverse_data = "InverseData"), function(object, solution, inverse_data) {

  attr <- list(SOLVE_TIME = solution$info$run_time,
               EXTRA_STATS =  solution)

  # Map OSQP statuses back to CVXR statuses.
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
#' @importFrom osqp osqp
#' @describeIn OSQP Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "OSQP", function(object, data, warm_start, verbose, solver_opts, solver_cache = NULL) {
  if(is.null(solver_cache))
    solver_cache  <- new.env(parent = emptyenv())

  P <- data[[P_KEY]]
  q <- data[[Q_KEY]]
  A <- Matrix(do.call(rbind, list(data[[A_KEY]], data[[F_KEY]])), sparse = TRUE)
  data$Ax <- A
  uA <- c(data[[B_KEY]], data[[G_KEY]])
  data$u <- uA
  lA <- c(data[[B_KEY]], rep(-Inf, length(data[[G_KEY]])))
  data$l <- lA

  # Overwrite defaults of eps_abs = eps_rel = 1e-3, max_iter = 4000.
  if(is.null(solver_opts$eps_abs))
      solver_opts$eps_abs <- 1e-5
  if(is.null(solver_opts$eps_rel))
    solver_opts$eps_rel <- 1e-5
  if(is.null(solver_opts$max_iter))
    solver_opts$max_iter <- 10000

  # Use cached data.
  if(warm_start && name(object) %in% names(solver_cache)) {  ## x %in% NULL will return FALSE always, so ok
    cache <- solver_cache[[name(object)]]
    solver <- cache[[1L]]
    old_data <- cache[[2L]]
    results <- cache[[3L]]

    new_args <- list()
    for(key in c("q", "l", "u")) {
      if(any(data[[key]] != old_data[[key]]))
        new_args[[key]] <- data[[key]]
    }

    factorizing <- FALSE
    if(any(dim(P) != old_data[[P_KEY]]) || any(P != old_data[[P_KEY]])) {
      P_triu <- Matrix::triu(P)
      new_args$Px <- P_triu
      factorizing <- TRUE
    }
    if(any(dim(A) != old_data$Ax) || any(A != old_data$Ax)) {
      new_args$Ax <- A
      factorizing <- TRUE
    }

    if(length(new_args) > 0)
      do.call(solver$Update, new_args)

    # Map OSQP statuses back to CVXR statuses.
    status <- status_map(object, results$info$status_val, default = SOLVER_ERROR)
    if(status == OPTIMAL)
      solver$WarmStart(results$x, results$y)

    # Polish if factorizing.
    if(is.null(solver_opts$polish))
      solver_opts$polish <- TRUE
    solver$UpdateSettings(c(list(verbose = verbose), solver_opts))
  } else {
    # Initialize and solve problem.
    if(is.null(solver_opts$polish))
      solver_opts$polish <- TRUE
    solver <- osqp::osqp(P, q, A, lA, uA, c(list(verbose = verbose), solver_opts))
  }

  results <- solver$Solve()

  if(!is.null(solver_cache)) {
    solver_cache[[name(object)]] <- list(solver, data, results)
  }
  return(list(results, solver_cache))
})
