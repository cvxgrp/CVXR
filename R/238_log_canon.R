## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/log_canon.py
#'
#' Dgp2Dcp canonicalizer for the log atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the log atom of a DGP expression,
#' where the returned expression is the log of the original expression..
Dgp2Dcp.log_canon <- function(expr, args) {
  return(list(Log(args[[1]]), list()))
}

