## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/pf_eigenvalue_canon.py
#'
#' Dgp2Dcp canonicalizer for the spectral radius atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the spectral radius atom of a DGP expression,
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.pf_eigenvalue_canon <- function(expr, args) {
  X <- args[[1]]
  # rho(X) <= lambda iff there exists v s.t. Xv <= lambda v.
  # v and lambd represent log variables, hence no positivity constraints.
  lambd <- Variable()
  v <- Variable(nrow(X))
  lhs <- X %*% v
  rhs <- lambd*v
  lhs <- Dgp2Dcp.mulexpression_canon(lhs, lhs@args)[[1]]
  rhs <- Dgp2Dcp.mul_canon(rhs, rhs@args)[[1]]
  return(list(lambd, list(lhs <= rhs)))
}

