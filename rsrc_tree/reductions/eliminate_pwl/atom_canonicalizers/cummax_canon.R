## CVXPY SOURCE: cvxpy/reductions/eliminate_pwl/cummax_canon.py
#'
#' EliminatePwl canonicalizer for the cumulative max atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed from a cumulative max atom where the objective
#' function consists of the variable Y which is of the same
#' dimension as the original expression and the constraints
#' consist of row/column constraints depending on the axis
EliminatePwl.cummax_canon <- function(expr, args) {
  X <- args[[1]]
  axis <- expr@axis

  # Implicit O(n) definition:
  # Y_{k} = maximum(Y_{k-1}, X_k)
  # Y <- Variable(dim(expr))
  Y <- new("Variable", dim = dim(expr))
  constr <- list(X <= Y)
  if(axis == 2) {
    if(nrow(expr) == 1)
      return(list(X, list()))
    else
      constr <- c(constr, list(Y[1:(nrow(Y) - 1), ] <= Y[2:nrow(Y),]))
  } else {
    if(ncol(expr) == 1)
      return(list(X, list()))
    else
      constr <- c(constr, list(Y[, 1:(ncol(Y) - 1)], Y[, 2:ncol(Y)]))
  }
  return(list(Y, constr))
}
