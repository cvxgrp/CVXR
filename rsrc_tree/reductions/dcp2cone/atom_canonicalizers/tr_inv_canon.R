## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/tr_inv_canon.py

#'
#' Dcp2Cone canonicalizer for the trinv atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a trinv atom
#' where the objective function consists of the variable t
#' that is of the same dimension as the original expression
#' with specified constraints in the function.
Dcp2Cone.tr_inv_canon <- function(expr, args) {
  X <- args[[1]]
  n <- nrow(X)
  su <- NULL

  constraints <- list()
  for(i in seq_len(n)) {
    ei <- matrix(0, nrow = n, ncol = 1)
    ei[i] <- 1.0
    ui <- Variable(1, 1)
    R <- bmat(list(list(X, ei),
                   list(t(ei), ui)))
    constraints <- c(constraints, list(R %>>% 0))
    if(is.null(su))
      su <- ui
    else
      su <- su + ui
  }
  return(list(su, constraints))
}
