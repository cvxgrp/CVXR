## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/normNuc_canon.py

#'
#' Dcp2Cone canonicalizer for the nuclear norm atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a nuclear norm atom,
#' where the objective function consists of .5 times the trace of
#' a matrix X of size m+n by m+n where the constraint consist of
#' the top right corner of the matrix being the original matrix.
Dcp2Cone.normNuc_canon <- function(expr, args) {
  A <- args[[1]]
  A_dim <- dim(A)
  m <- A_dim[1]
  n <- A_dim[2]

  # Create the equivalent problem:
  #   minimize (lo.trace(U) + lo.trace(V))/2
  #   subject to:
  #            [U A; t(A) V] is positive semidefinite
  U <- Variable(m, m, symmetric = TRUE)
  V <- Variable(n, n, symmetric = TRUE)
  X <- bmat(list(list(U, A),
                 list(t(A), V)))
  constraints <- list(X %>>% 0)
  trace_value <- 0.5 * (matrix_trace(U) + matrix_trace(V))
  return(list(trace_value, constraints))
}
