## CVXPY SOURCE: cvxpy/reductions/solvers/conic_solvers/glpk_mi_conif.py

#' An interface for the GLPK MI solver.
#'
#' @name GLPK_MI-class
#' @aliases GLPK_MI
#' @rdname GLPK_MI-class
#' @export
setClass("GLPK_MI", prototype = list(MIP_CAPABLE = TRUE,
                                     SUPPORTED_CONSTRAINTS = supported_constraints(ConicSolver()),
                                     MI_SUPPORTED_CONSTRAINTS = supported_constraints(ConicSolver())),
         contains = "GLPK")

#' @rdname GLPK_MI-class
#' @export
GLPK_MI <- function() { new("GLPK_MI") }

# Map of GLPK_MI status to CVXR status.
#' @param solver,object,x A \linkS4class{GLPK_MI} object.
#' @param status A status code returned by the solver.
#' @describeIn GLPK_MI Converts status returned by the GLPK_MI solver to its respective CVXPY status.
setMethod("status_map", "GLPK_MI", function(solver, status) {
  if(status == 5)
    OPTIMAL
  else if(status == 2)
    SOLUTION_PRESENT #Not sure if feasible is the same thing as this
  else if(status == 3 | status == 4)
    INFEASIBLE
  else if(status == 1 | status == 6)
    UNBOUNDED
  else
    stop("GLPK_MI status unrecognized: ", status)
})

#' @describeIn GLPK_MI Returns the name of the solver.
setMethod("name", "GLPK_MI", function(x) { GLPK_MI_NAME })

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn GLPK_MI Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "GLPK_MI", function(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
  if(verbose)
    solver_opts$verbose <- verbose
  solver_opts$canonicalize_status <- FALSE
  if(!is.null(feastol)){
    warning("A value has been set for feastol, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }
  if(!is.null(reltol)){
    warning("A value has been set for reltol, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }
  if(!is.null(abstol)){
    warning("A value has been set for abstol, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }
  if(!is.null(num_iter)){
    warning("A value has been set for num_iter, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }

  # Construct problem data.
  c <- data[[C_KEY]]
  dims <- data[[ConicSolver()@dims]]
  nvar <- length(c)
  A <- data[[A_KEY]]
  b <- data[[B_KEY]]
  if(nrow(A) == 0)
    A <- Matrix(0, nrow = 0, ncol = length(c))

  G <- data[[G_KEY]]
  h <- data[[H_KEY]]
  if(nrow(G) == 0)
    G <- Matrix(0, nrow = 0, ncol = length(c))

  mat <- rbind(A, G)
  rhs <- c(b, h)
  bounds <- list(lower = list(ind = seq_along(c), val = rep(-Inf, nvar)))
  types <- rep("C", nvar)
  bools <- data[[BOOL_IDX]]
  ints <- data[[INT_IDX]]
  if (length(bools) > 0) {
    types[bools] <- "B"
  }
  if (length(ints) > 0) {
    types[ints] <- "I"
  }

  results_dict <- Rglpk::Rglpk_solve_LP(obj = c,
                                        mat = slam::as.simple_triplet_matrix(mat),
                                        dir = c(rep("==", dims@zero),
                                                rep("<=", dims@nonpos)),
                                        rhs = rhs,
                                        bounds = bounds,
                                        types = types,
                                        control = solver_opts,
                                        max = FALSE)

  # Convert results to solution format.
  solution <- list()

  solution[[STATUS]] <- status_map(object, results_dict$status)
  if(solution[[STATUS]] %in% SOLUTION_PRESENT) {
    ## Get primal variable values
    solution[[PRIMAL]] <- results_dict$solution
    ## Get objective value
    solution[[VALUE]] <- results_dict$optimum
    # solution[[EQ_DUAL]] <- results_dict$auxiliary[[1]]   # TODO: How do we get the dual variables?
    # solution[[INEQ_DUAL]] <- results_dict$auxiliar[[2]]
    solution[[EQ_DUAL]] <- list()
    solution[[INEQ_DUAL]] <- list()
  }
  solution
})
