## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/kl_div_canon.py
#'
#' Dcp2Cone canonicalizer for the KL Divergence atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a KL divergence atom
#' where t is the objective function with the ExpCone constraints.
Dcp2Cone.kl_div_canon <- function(expr, args) {
  expr_dim <- dim(expr)
  x <- promote(args[[1]], expr_dim)
  y <- promote(args[[2]], expr_dim)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)
  constraints <- list(ExpCone(t, x, y))
  obj <- y - x - t
  return(list(obj, constraints))
}

