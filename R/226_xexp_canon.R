## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/xexp_canon.py

#'
#' Dcp2Cone canonicalizer for the xexp atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a xexp atom
#' where the objective function consists of the variable t
#' that is of the same dimension as the original expression
#' with specified constraints in the function.
Dcp2Cone.xexp_canon <- function(expr, args) {
  x <- args[[1]]
  u <- new("Variable", dim = dim(expr), nonneg = TRUE)
  t <- new("Variable", dim = dim(expr), nonneg = TRUE)
  power_expr <- Power(x, 2)
  canon <- Dcp2Cone.power_canon(power_expr, power_expr@args)
  power_obj <- canon[[1]]
  constraints <- canon[[2]]

  constraints <- c(constraints, list(ExpCone(u, x, t), u >= power_obj, x >= 0))
  return(list(t, constraints))
}
