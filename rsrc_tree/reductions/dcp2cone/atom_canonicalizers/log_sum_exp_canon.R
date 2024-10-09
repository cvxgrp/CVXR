## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/log_sum_exp_canon.py

#'
#' Dcp2Cone canonicalizer for the log sum of the exp atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from the log sum
#' of the exp atom where the objective is the t variable
#' and the constraints consist of the ExpCone constraints and
#' requiring t to be less than a matrix of ones of the same size.
Dcp2Cone.log_sum_exp_canon <- function(expr, args) {
  x <- args[[1]]
  x_dim <- dim(x)
  expr_dim <- dim(expr)
  axis <- expr@axis
  keepdims <- expr@keepdims
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)

  # log(sum(exp(x))) <= t <=> sum(exp(x-t)) <= 1.
  if(is.na(axis))   # shape = c(1,1)
    promoted_t <- promote(t, x_dim)
  else if(axis == 2)   # shape = c(1,n)
    promoted_t <- Constant(matrix(1, nrow = x_dim[1], ncol = 1)) %*% reshape_expr(t, c(1, x_dim[2:length(x_dim)]))
  else   # shape = c(m,1)
    promoted_t <- reshape_expr(t, c(x_dim[1:(length(x_dim)-1)], 1)) %*% Constant(matrix(1, nrow = 1, ncol = x_dim[2]))

  exp_expr <- Exp(x - promoted_t)
  canon <- Dcp2Cone.exp_canon(exp_expr, exp_expr@args)
  obj <- lu.sum_entries(canon[[1]], axis = axis, keepdims = keepdims)
  if(is.null(expr_dim))
    ones <- Constant(1)
  else
    ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- c(canon[[2]], obj <= ones)
  return(list(t, constraints))
}

