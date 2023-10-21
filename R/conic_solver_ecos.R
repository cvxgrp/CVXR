#'
#' Utility method for formatting a ConeDims instance into a dictionary
#' that can be supplied to ECOS.
#'
#' @param cone_dims A \linkS4class{ConeDims} instance.
#' @return A dictionary of cone dimensions as a list
#' @export
ECOS.dims_to_solver_dict <- function(cone_dims) {
  list(l = as.integer(cone_dims@nonpos),
       q = lapply(cone_dims@soc, as.integer),
       e = as.integer(cone_dims@exp))
}

#'
#' An interface for the ECOS solver
#'
#' @name ECOS-class
#' @aliases ECOS
#' @rdname ECOS-class
#' @export
setClass("ECOS", contains = "ConicSolver")


setMethod("initialize", "ECOS",
          function(.Object, ...) {
            ##.Object <- callNextMethod(.Object, ...)
            ## Ensure EXP_CONE_ORDER is set
            .Object@EXP_CONE_ORDER <- c(0L, 1L, 2L)
            .Object@SUPPORTED_CONSTRAINTS <- c(.Object@SUPPORTED_CONSTRAINTS, "SOC", "ExpCone")
            .Object
          })

## The function below not needed usually!
## #' @rdname ECOS-class
## #' @export
## ECOS <- function() { new("ECOS") }

# EXITCODES from ECOS
# ECOS_OPTIMAL  (0)   Problem solved to optimality
# ECOS_PINF     (1)   Found certificate of primal infeasibility
# ECOS_DINF     (2)   Found certificate of dual infeasibility
# ECOS_INACC_OFFSET (10)  Offset exitflag at inaccurate results
# ECOS_MAXIT    (-1)  Maximum number of iterations reached
# ECOS_NUMERICS (-2)  Search direction unreliable
# ECOS_OUTCONE  (-3)  s or z got outside the cone, numerics?
# ECOS_SIGINT   (-4)  solver interrupted by a signal/ctrl-c
# ECOS_FATAL    (-7)  Unknown problem in solver

# Map of ECOS status to CVXR status.
#' @param solver,object,x A \linkS4class{ECOS} object.
#' @param status A status code returned by the solver.
#' @describeIn ECOS Converts status returned by the ECOS solver to its respective CVXPY status.
setMethod("status_map", "ECOS", function(solver, status) {
  if(status == 0)
    return(OPTIMAL)
  else if(status == 1)
    return(INFEASIBLE)
  else if(status == 2)
    return(UNBOUNDED)
  else if(status == 10)
    return(OPTIMAL_INACCURATE)
  else if(status == 11)
    return(INFEASIBLE_INACCURATE)
  else if(status == 12)
    return(UNBOUNDED_INACCURATE)
  else if(status %in% c(-1, -2, -3, -4, -7))
    return(SOLVER_ERROR)
  else
    stop("ECOS status unrecognized: ", status)
})

#' @describeIn ECOS Imports the solver
setMethod("import_solver", "ECOS", function(solver) { requireNamespace("ECOSolveR", quietly = TRUE) })

#' @describeIn ECOS Returns the name of the solver
setMethod("name", "ECOS", function(solver) { ECOS_NAME })

#' @param problem A \linkS4class{Problem} object.
#' @describeIn ECOS Returns a new problem and data for inverting the new solution.
setMethod("perform", signature(object = "ECOS", problem = "Problem"), function(object, problem) {
  ## CHECK: SHOULD problem not be of class ParamConeProg ? ASSUMING so, below because accepts method for
  ## ConicSolver checks for that.

  data <- list()
  inv_data <- list(object@VAR_ID = id(problem@x))
  # Format constraints
  #
  # ECOS requires constraints to be specified in the following order:
  # 1. zero cone
  # 2. non-negative orthant
  # 3. soc
  # 4. exponential

  if (! problem@formatted) {
    problem <- format_constr(object, problem, object@EXP_CONE_ORDER)
  }

  data[[PARAM_PROB]] <- problem
  data[[object@dims]] <- inv_data[[object@dims]] <- problem@cone_dims

  constr_map <- problem@constr_map
  inv_data[[EQ_CONSTR]] <- constr_map$ZeroConstraint
  inv_data[[NEQ_CONSTR]] <- c(constr_map$NonNegConstraint, constr_map$SOC, constr_map$ExpCone)

  len_eq_ind <- seq_len(problem@cone_dims$zero)

  output_list <- apply_parameters(problem)
  data[[C_KEY]] <- output_list$c
  inv_data[[OFFSET_KEY]] <- output_list$d

  mat <- output_list$A[len_eq_ind, ]
  if (nrow(mat) == 0L) {
    data[[A_KEY]] <- NULL
  } else {
    data[[A_KEY]] <- -mat
  }

  vec <- as.numeric(output_list$b[len_eq_ind, ])
  if (length(vec) == 0L) {
    data[[B_KEY]] <- NULL
  } else {
    data[[B_KEY]] <- vec
  }

  mat <- output_list$A[-len_eq_seq, ]
  if (nrow(mat) == 0L) {
    data[[G_KEY]] <- NULL
  } else {
    data[[G_KEY]] <- -mat
  }

  vec <- as.numeric(output_list$b[-len_eq_ind, ])
  if (length(vec) == 0L) {
    data[[H_KEY]] <- NULL
  } else {
    data[[H_KEY]] <- vec
  }

  ## CHECK why do you need to return object? Leave it out as it's not been touched!!
  ## return(list(object, data = data, inv_data = inv_data))
  list(data = data, inv_data = inv_data)
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose An integer number indicating level of solver verbosity.
#' @param feastol The feasible tolerance on the primal and dual residual.
#' @param reltol The relative tolerance on the duality gap.
#' @param abstol The absolute tolerance on the duality gap.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ReductionSolver Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "ECOS", function(object,
                                             data,
                                             warm_start,
                                             verbose,
                                             feastol,
                                             reltol,
                                             abstol,
                                             num_iter,
                                             solver_opts,
                                             solver_cache = new.env(parent = emptyenv()))

  if (is.null(feastol)) {
      feastol <- SOLVER_DEFAULT_PARAM$ECOS$feastol
  }
  if (is.null(reltol)) {
      reltol <- SOLVER_DEFAULT_PARAM$ECOS$reltol
  }
  if (is.null(abstol)) {
      abstol <- SOLVER_DEFAULT_PARAM$ECOS$abstol
  }
  if (is.null(num_iter)) {
      num_iter <- SOLVER_DEFAULT_PARAM$ECOS$maxit
  }
  ecos_opts <- ECOSolveR::ecos.control(maxit = as.integer(num_iter),
                                       feastol = feastol,
                                       reltol = reltol,
                                       abstol = abstol,
                                       verbose = as.integer(verbose))
  ecos_opts[names(solver_opts)] <- solver_opts

  ECOSolveR::ECOS_csolve(c = data[[C_KEY]],
                         G = data[[G_KEY]],
                         h = data[[H_KEY]],
                         dims = ECOS.dims_to_solver_dict(data[[object@dims]]),
                         A = data[[A_KEY]],
                         b = data[[B_KEY]],
                         control = ecos_opts)
})

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn ECOS Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "ECOS", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
  status <- status_map(object, solution$retcodes[["exitFlag"]])

  # Timing data.
  attr <- list()
  attr[[SOLVE_TIME]] <- solution$timing[["tsolve"]]
  attr[[SETUP_TIME]] <- solution$timing[["tsetup"]]
  attr[[NUM_ITERS]] <- solution$retcodes[["iter"]]

  if(status %in% SOLUTION_PRESENT) {
    primal_val <- solution$summary[["pcost"]]
    opt_val <- primal_val + inverse_data[[OFFSET]]
    primal_vars <- list()
    var_id <- inverse_data[[object@var_id]]
    primal_vars[[as.character(var_id)]] <- as.matrix(solution$x)
    eq_dual <- get_dual_values(solution$y, extract_dual_value, inverse_data[[object@eq_constr]])
    leq_dual <- get_dual_values(solution$z, extract_dual_value, inverse_data[[object@neq_constr]])
    eq_dual <- utils::modifyList(eq_dual, leq_dual)
    dual_vars <- eq_dual

    return(Solution(status, opt_val, primal_vars, dual_vars, attr))
  } else
    return(failure_solution(status))
})

