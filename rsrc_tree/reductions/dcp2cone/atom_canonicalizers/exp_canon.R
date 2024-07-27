## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/exp_canon.py
#'
#' Dcp2Cone canonicalizer for the exponential atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from an exponential atom
#' where the objective function is the variable t with an ExpCone constraint.
Dcp2Cone.exp_canon <- function(expr, args) {
  expr_dim <- dim(expr)
  x <- promote(args[[1]], expr_dim)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)
  if(is.null(expr_dim))
    ones <- Constant(1)
  else
    ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- list(ExpCone(x, ones, t))
  return(list(t, constraints))
}

