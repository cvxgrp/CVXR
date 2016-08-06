#'
#' The Minimize class.
#'
#' This class represents a minimization problem.
#'
#' @slot expr The expression to minimize.
#' @aliases Minimize
#' @export
Minimize <- setClass("Minimize", representation(expr = "ConstValORExpr"))

setMethod("initialize", "Minimize", function(.Object, expr) {
    .Object@expr <- as.Constant(expr)
    if(!all(size(.Object@expr) == c(1,1)))
      stop("The objective must resolve to a scalar")
    return(.Object)
})

setMethod("get_data", "Minimize", function(object) { list() })
setMethod("canonical_form", "Minimize", function(object) { canonicalize(object) })
setMethod("canonicalize", "Minimize", function(object) { canonical_form(object@expr) })
setMethod("variables", "Minimize", function(object) { variables(object@expr) })
setMethod("parameters", "Minimize", function(object) { parameters(object@expr) })
setMethod("is_dcp", "Minimize", function(object) { is_convex(object@expr) })

#'
#' The Maximize class.
#'
#' This class represents a maximization problem.
#'
#' @slot expr The expression to maximize.
#' @aliases Maximize
#' @export
Maximize <- setClass("Maximize", contains = "Minimize")

setMethod("-", signature(e1 = "Minimize", e2 = "missing"), function(e1, e2) { Maximize(expr = - e1@expr) })
setMethod("+", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { Minimize(e1@expr + e2@expr) })
setMethod("+", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { stop("Problem does not follow DCP rules") })
setMethod("-", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { e1 + -e2 })
setMethod("-", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { e1 + -e2 })
setMethod("*", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) {
  if(e2 >= 0) Minimize(expr = e1@expr * e2) else Maximize(expr = e1@expr * e2)
})
setMethod("*", signature(e1 = "Maximize", e2 = "numeric"), function(e1, e2) {
  if(e2 < 0) Minimize(expr = e1@expr * e2) else Maximize(expr = e1@expr * e2)
})
setMethod("/", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) { e1 * (1.0/e2) })

setMethod("-", signature(e1 = "Maximize", e2 = "missing"), function(e1, e2) { Minimize(-e1@expr) })
setMethod("+", signature(e1 = "Maximize", e2 = "Maximize"), function(e1, e2) { Maximize(expr = e1@expr + e2@expr) })
setMethod("+", signature(e1 = "Maximize", e2 = "Minimize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

setMethod("canonicalize", "Maximize", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constraints <- canon[[2]]
  list(neg_expr(obj), constraints)
})

setMethod("is_dcp", "Maximize", function(object) { is_concave(object@expr) })

#'
#' The Problem class.
#'
#' This class represents the convex optimization problem.
#'
#' @slot objective The expression to minimize or maximize.
#' @slot constraints (Optional) A list of constraints on the problem variables.
#' @aliases Problem
#' @export
.Problem <- setClass("Problem", representation(objective = "Minimize", constraints = "list", value = "numeric", status = "character", .cached_data = "list", .separable_problems = "list"),
                    prototype(constraints = list(), value = NA_real_, status = NA_character_, .cached_data = list(), .separable_problems = list()),
                    validity = function(object) {
                      if(!(class(object@objective) %in% c("Minimize", "Maximize")))
                        stop("[Problem: objective] objective must be of class Minimize or Maximize")
                      if(!is.na(object@value))
                        stop("[Problem: value] value should not be set by user")
                      if(!is.na(object@status))
                        stop("[Problem: status] status should not be set by user")
                      if(length(object@.cached_data) > 0)
                        stop("[Problem: .cached_data] .cached_data is an internal slot and should not be set by user")
                      if(length(object@.separable_problems) > 0)
                        stop("[Problem: .separable_problems] .separable_problems is an internal slot and should not be set by user")
                      return(TRUE)
                    })

Problem <- function(objective, constraints = list()) {
  .Problem(objective = objective, constraints = constraints)
}

setMethod("initialize", "Problem", function(.Object, ..., objective, constraints = list(), value = NA_real_, status = NA_character_, .cached_data = list(), .separable_problems = list()) {
  .Object@objective <- objective
  .Object@constraints <- constraints
  .Object@value <- value
  .Object@status <- status

  # Cached processed data for each solver.
  .Object@.cached_data <- list()
  .Object <- .reset_cache(.Object)

  # List of separable (sub)problems
  .Object@.separable_problems <- .separable_problems
  .Object
})

CachedProblem <- function(objective, constraints) { list(objective = objective, constraints = constraints) }
SolveResult <- function(opt_value, status, primal_values, dual_values) { list(opt_value = opt_value, status = status, primal_values = primal_values, dual_values = dual_values) }

setMethod(".reset_cache", "Problem", function(object) {
  for(solver_name in SOLVERS)
    object@.cached_data[[solver_name]] <- ProblemData()
  object@.cached_data[[PARALLEL]] <- CachedProblem(NA, NULL)
  object
})

setMethod("is_dcp", "Problem", function(object) {
  is_dcp_list <- lapply(c(object@objective, object@constraints), is_dcp)
  all(unlist(is_dcp_list))
})

setMethod("+", signature(e1 = "Problem", e2 = "missing"), function(e1, e2) { Problem(objective = e1@objective, constraints = e1@constraints) })
setMethod("-", signature(e1 = "Problem", e2 = "missing"), function(e1, e2) { Problem(objective = -e1@objective, constraints = e1@constraints) })
setMethod("+", signature(e1 = "Problem", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e1@objective + e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})
setMethod("-", signature(e1 = "Problem", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e1@objective - e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})
setMethod("*", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) {
  Problem(objective = e1@objective * e2, constraints = e1@constraints)
})
setMethod("*", signature(e1 = "numeric", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e2 * e1@objective, constraints = e1@constraints)
})
setMethod("/", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) {
  Problem(objective = e1@objective * (1.0/e2), constraints = e1@constraints)
})

setMethod("canonicalize", "Problem", function(object) {
  obj_canon <- canonical_form(object@objective)
  canon_constr <- sapply(object@constraints, function(x) { canonical_form(x)[[2]] })
  list(obj_canon[[1]], c(obj_canon[[2]], canon_constr))
})

setMethod("variables", "Problem", function(object) {
  vars_ <- variables(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { variables(constr) })
  unique(flatten_list(c(vars_, constrs_)))
})

setMethod("parameters", "Problem", function(object) {
  params <- parameters(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { parameters(constr) })
  unique(flatten_list(c(params, constrs_)))
})

setMethod("cvxr_solve", "Problem", function(object, solver = NULL, ignore_dcp = FALSE, warm_start = FALSE, verbose = FALSE, parallel = FALSE, ...) {
  if(!is_dcp(object)) {
    if(ignore_dcp)
      print("Problem does not follow DCP rules. Solving a convex relaxation.")
    else
      stop("Problem does not follow DCP rules.")
  }

  canon <- canonicalize(object)
  objective <- canon[[1]]
  constraints <- canon[[2]]

  # TODO: Solve in parallel

  print("Calling CVXcanon")
  if(is(object@objective, "Minimize")) {
    sense <- "Minimize"
    canon_objective <- objective
  } else {
    sense <- "Maximize"
    canon_objective <- neg_expr(objective)  # preserve sense
  }

  solve(sense, canon_objective, constraints, verbose, ...)

  ## Start of section commented out by Naras

  # Choose a solver/check the chosen solver.
  ## if(is.null(solver))
  ##   solver <- Solver.choose_solver(constraints)
  ## else if(is(solver, "Solver") && name(solver) %in% SOLVERS)
  ##   validate_solver(solver, constraints)
  ## else if(is.character(solver) && solver %in% SOLVERS) {
  ##   solver <- new(solver)
  ##   validate_solver(solver, constraints)
  ## } else
  ##   stop("Unknown solver.")

  ## sym_data <- get_sym_data(solver, objective, constraints, object@.cached_data)

  ## # Presolve couldn't solve the problem.
  ## if(is.na(sym_data@.presolve_status)) {
  ##   results_dict <- cvxr_solve_int(solver, objective, constraints, object@.cached_data, warm_start, verbose, ...)
  ## # Presolve determined the problem was unbounded or infeasible.
  ## } else {
  ##   results_dict <- list()
  ##   results_dict[[STATUS]] <- sym_data@.presolve_status
  ## }

  ## object@.update_problem_state(results_dict, sym_data, solver)
  ## object@value

  ## End of section commented out by Naras

})
