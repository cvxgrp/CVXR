## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/power_canon.py

#'
#' Dgp2Dcp canonicalizer for the power atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the power atom of a DGP expression,
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.power_canon <- function(expr, args) {
  # y = log(x); x^p --> exp(y^p) --> p*log(exp(y)) = p*y.
  return(list(expr@p*args[[1]], list()))
}

