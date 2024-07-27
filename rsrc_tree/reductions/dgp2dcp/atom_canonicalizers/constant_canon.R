## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/constant_canon.py

#'
#' Dgp2Dcp canonicalizer for the constant atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the constant atom of a DGP expression,
#' where the returned expression is the DCP equivalent resulting
#' from the log of the expression.
Dgp2Dcp.constant_canon <- function(expr, args) {
  # args <- list()
  return(list(Constant(log(value(expr))), list()))
}

