## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/lambda_sum_largest_canon.py
#'
#' Dcp2Cone canonicalizer for the largest lambda sum atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a lambda sum of the k
#' largest elements atom where k*t + trace(Z) is the objective function.
#' t denotes the variable subject to constraints and Z is a PSD matrix variable
#' whose dimensions consist of the length of the vector at hand. The constraints
#' require the the diagonal matrix of the vector to be symmetric and PSD.
Dcp2Cone.lambda_sum_largest_canon <- function(expr, args) {
  # S_k(X) denotes lambda_sum_largest(X, k)
  # t >= k S_k(X - Z) + trace(Z), Z is PSD
  # implies
  # t >= ks + trace(Z)
  # Z is PSD
  # sI >= X - Z (PSD sense)
  # which implies
  # t >= ks + trace(Z) >= S_k(sI + Z) >= S_k(X)
  # We use the fact that
  # S_k(X) = sup_{sets of k orthonormal vectors u_i}\sum_{i}u_i^T X u_i
  # and if Z >= X in PSD sense then
  # \sum_{i}u_i^T Z u_i >= \sum_{i}u_i^T X u_i
  #
  # We have equality when s = lambda_k and Z diagonal
  #  with Z_{ii} = (lambda_i - lambda_k)_+

  X <- expr@args[[1]]
  k <- expr@k
  # Z <- Variable(c(nrow(X), nrow(X)), PSD = TRUE)
  Z <- Variable(nrow(X), ncol(X), PSD = TRUE)
  canon <- Dcp2Cone.lambda_max_canon(expr, list(X - Z))
  obj <- canon[[1]]
  constr <- canon[[2]]
  obj <- k*obj + matrix_trace(Z)
  return(list(obj, constr))
}

