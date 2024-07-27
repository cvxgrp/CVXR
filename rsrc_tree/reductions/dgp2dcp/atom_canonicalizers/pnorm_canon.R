## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/pnorm_canon.py
#'
#' Dgp2Dcp canonicalizer for the p-norm atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the pnorm atom of a DGP expression,
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.pnorm_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@original_p
  if(is.null(dim(x)))
    x <- promote(p, c(1))
  if(is.na(expr@axis) || length(dim(x)) == 1) {
    x <- Vec(x)
    # hstack_args <- lapply(seq_len(size(x)), function(j) { x[j]^p })
    # return(list((1.0/p) * log_sum_exp(do.call("HStack", hstack_args)), list()))
    return(list((1.0/p) * log_sum_exp(x^p), list()))
  }

  if(expr@axis == 2)
    x <- t(x)

  rows <- list()
  for(i in seq_len(nrow(x))) {
    row <- x[i,]
    # hstack_args <- lapply(seq_len(size(row)), function(j) { row[j]^p })
    # rows <- c(rows, list((1.0/p)*log_sum_exp(do.call("HStack", hstack_args))))
    rows <- c(rows, list((1.0/p) * log_sum_exp(row^p)))
  }
  return(list(do.call("VStack", rows), list()))
}
