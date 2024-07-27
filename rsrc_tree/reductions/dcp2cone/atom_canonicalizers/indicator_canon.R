## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/indicator_canon.py
#'
#' Dcp2Cone canonicalizer for the indicator atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from an indicator atom and
#' where 0 is the objective function with the given constraints
#' in the function.
Dcp2Cone.indicator_canon <- function(expr, args) {
  return(list(Constant(0), args))
}

