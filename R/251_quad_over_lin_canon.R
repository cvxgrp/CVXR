## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/quad_over_lin_canon.py
#'
#' Dgp2Dcp canonicalizer for the quadratic over linear term atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the quadratic over linear atom of a
#' DGP expression, where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.quad_over_lin_canon <- function(expr, args) {
  summed <- Dgp2Dcp.explicit_sum(2*args[[1]])
  numerator <- Dgp2Dcp.add_canon(summed, summed@args)
  return(list(numerator - args[[2]], list()))
}

