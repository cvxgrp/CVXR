## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/div_canon.py
#'
#' Dgp2Dcp canonicalizer for the division atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the division atom of a DGP expression,
#' where the returned expression is the log transformed DCP equivalent.
Dgp2Dcp.div_canon <- function(expr, args) {
  # expr <- NULL
  # x / y == x * y^(-1)
  return(list(args[[1]] - args[[2]], list()))
}
