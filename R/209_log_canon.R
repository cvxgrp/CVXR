## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/log_canon.py
#'
#' Dcp2Cone canonicalizer for the log atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a log atom where
#' t is the objective function and the constraints consist of
#' ExpCone constraints
Dcp2Cone.log_canon <- function(expr, args) {
  x <- args[[1]]
  expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)
  if(is.null(expr_dim))
    ones <- Constant(1)
  else
    ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  # TODO: ExpCone requires each of its inputs to be a Variable; is this something that we want to change?
  constraints <- list(ExpCone(t, ones, x))
  return(list(t, constraints))
}
