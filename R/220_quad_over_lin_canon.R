## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/quad_over_lin_canon.py

#'
#' Dcp2Cone canonicalizer for the quadratic over linear term atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a quadratic over linear
#' term atom where the objective function consists of a one
#' dimensional variable t with SOC constraints.
Dcp2Cone.quad_over_lin_canon <- function(expr, args) {
  # quad_over_lin := sum_{ij} X^2_{ij} / y
  x <- args[[1]]
  y <- flatten(args[[2]])

  # Pre-condition: dim = c()
  t <- Variable(1)

  # (y+t, y-t, 2*x) must lie in the second-order cone, where y+t is the scalar part
  # of the second-order cone constraint
  # BUG: In Python, flatten produces single dimension (n,), but in R, we always treat
  # these as column vectors with dimension (n,1), necessitating the use of VStack.
  # constraints <- list(SOC(t = y+t, X = HStack(y-t, 2*flatten(x)), axis = 2))
  constraints <- list(SOC(t = y+t, X = VStack(y-t, 2*flatten(x)), axis = 2))
  return(list(t, constraints))
}
