BisectionData <- list("BisectionData" = list("feas_problem", "param", "tighten_lower", "tighten_upper"))

Dqcp2Dcp.get_lazy_and_real_constraints <- function(constraints) {
  lazy_constraints <- list()
  real_constraints <- list()
  for(c in constraints) {
    if(is.function(c))
      lazy_constraints <- c(lazy_constraints, c)
    else
      real_constraints <- c(real_constraints, c)
  }
  return(list(lazy_constraints, real_constraints))
}

#'
#' Reduce DQCP Problem to Parametrized DCP Problem
#'
#' This reduction takes as input a DQCP problem and returns a parameterized
#' DCP problem that can be solved by bisection. Some of the constraints might
#' be lazy, i.e., callables that return a constraint when called. The problem
#' will only be DCP once the lazy constraints are replaced with actual
#' constraints.
#' 
#' Problems emitted by this reduction can be solved with the bisect function.
#'
#' @rdname Dqcp2Cone-class
.Dqcp2Cone <- setClass("Dqcp2Cone", prototype(.bisection_data = NULL), contains = "Canonicalization")
Dqcp2Cone <- function(problem = NULL) { .Dqcp2Cone(problem = problem) }

setMethod("initialize", "Dqcp2Cone", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = Dcp2Cone.CANON_METHODS)
})

#' @param object A \linkS4class{Dqcp2Cone} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Dqcp2Cone A problem is accepted if it is a minimization and is DQCP.
setMethod("accepts", signature(object = "Dqcp2Cone", problem = "Problem"), function(object, problem) {
  inherits(problem@objective, "Minimize") && is_dqcp(problem)
})

#' @param solution A \linkS4class{Solution} object to invert.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn D2cp2Cone Returns a solution to the original problem given the inverse data.
setMethod("invert", signature(object = "Dqcp2Cone", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  pvars <- list()
  for(vid in inverse_data@id_map) {
    if(vid in solution@primal_vars)
      pvars[[vid]] <- solution@primal_vars[[vid]]
    else
      # Variable was optimized out because it was unconstrained.
      pvars[[vid]] <- 0.0
  }
  return(Solution(solution@status, solution@opt_val, pvars, list(), solution@attr))
}

#' @describeIn Dqcp2Cone Recursively canonicalize the objective and every constraint.
setMethod("perform", signature(object = "Dqcp2Cone", problem = "Problem"), function(object, problem) {
  constraints <- list()
  for(constr in problem@constraints)
    constraints <- c(constraints, Dqcp2Dcp.canonicalize_constraint(object, constr))
  tmp <- Dqcp2Dcp.get_lazy_and_real_constraints(constraints)
  lazy <- tmp[[1]]
  real <- tmp[[2]]
  feas_problem <- Problem(Minimize(0), real)
  feas_problem@.lazy_constraints <- lazy
  
  objective <- problem@objective@expr
  if(is_nonneg(objective))
    t <- Parameter(nonneg = TRUE)
  else if(is_nonpos(objective))
    t <- Parameter(nonpos = TRUE)
  else
    t <- Parameter()
  constraints <- c(constraints, Dqcp2Dcp.canonicalize_constraint(objective <= t))
  
  tmp <- Dqcp2Dcp.get_lazy_and_real_constraints(constraints)
  lazy <- tmp[[1]]
  real <- tmp[[2]]
  param_problem <- Problem(Minimize(0), real)
  param_problem@.lazy_constraints <- lazy
  param_problem@.bisection_data <- BisectionData(feas_problem, t, tighten_fns(objective))
  return(list(param_problem, InverseData(problem)))
})

# TODO: Finish converting .canonicalize_tree, .canon_args, .canonicalize_constraint, 
# and other files in reductions/dqcp2dcp folder.
