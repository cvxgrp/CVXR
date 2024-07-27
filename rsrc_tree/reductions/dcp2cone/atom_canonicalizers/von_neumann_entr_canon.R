## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/von_neumann_entr_canon.py

#'
#' Dcp2Cone canonicalizer for the von Neumann entropy atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a von Neumann entropy atom
#' where the objective function consists of the variable t
#' that is of the same dimension as the original expression
#' with specified constraints in the function.
Dcp2Cone.von_neumann_entr_canon <- function(expr, args) {
  N <- args[[1]]
  if(!is_real(N))
    stop("N must be real")

  n <- nrow(N)
  x <- Variable(n)
  t <- Variable()

  # START code that applies to all spectral functions
  constrs <- list()
  if(n > 1) {
    for(r in seq(2, n)) {
      # lambda_sum_largest(N, r) <= sum(x[1:(r-1)])
      expr_r <- lambda_sum_largest(N, r)
      canon <- Dcp2Cone.lambda_sum_largest_canon(expr_r, expr_r@args)
      epi <- canon[[1]]
      cons <- canon[[2]]
      constrs <- c(constrs, cons)
      con <- NonPosConstraint(epi - sum(x[1:(r-1)]))
      constrs <- c(constrs, list(con))
    }
  }

  # trace(N) \leq sum(x)
  con <- (matrix_trace(N) == sum(x))
  constrs <- c(constrs, list(con))

  # trace(N) == sum(x)
  con <- ZeroConstraint(matrix_trace(N) - sum(x))
  constrs <- c(constrs, list(con))

  # x[1:(n-1)] >= x[2:n]
  #   x[1] >= x[2],  x[2] >= x[3], ...
  if(n > 1) {
    con <- NonPosConstraint(x[2:n] - x[1:(n - 1)])
    constrs <- c(constrs, list(con))
  }

  # END code that applies to all spectral functions

  # sum(entr(x)) >= t
  canon <- Dcp2Cone.entr_canon(x, list(x))
  hypos <- canon[[1]]
  entr_cons <- canon[[2]]
  constrs <- c(constrs, entr_cons)
  con <- NonPosConstraint(t - sum(hypos))
  constrs <- c(constrs, list(con))

  return(list(t, constrs))
}
