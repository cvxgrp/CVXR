## CVXPY SOURCE: cvxpy/reductions/eliminate_pwl/norm_inf_canon.py
#'
#' EliminatePwl canonicalizer for the infinite norm atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed by the infinite norm atom where the objective
#' function consists variable t of the same dimension as the
#' expression and the constraints consist of a vector
#' constructed by multiplying t to a vector of 1's
EliminatePwl.norm_inf_canon <- function(expr, args) {
  x <- args[[1]]
  axis <- expr@axis
  # expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = dim(expr))

  if(is.na(axis))   # dim(expr) = c(1,1)
    promoted_t <- promote(t, dim(x))
  else if(axis == 2)   # dim(expr) = c(1,n)
    promoted_t <- Constant(matrix(1, nrow = nrow(x), ncol = 1) %*% reshape_expr(t, c(1, ncol(x))))
  else   # shape = c(m,1)
    promoted_t <- reshape_expr(t, c(nrow(x), 1)) %*% Constant(matrix(1, nrow = 1, ncol = ncol(x)))

  return(list(t, list(x <= promoted_t, x + promoted_t >= 0)))
}

