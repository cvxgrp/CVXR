## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/pnorm_canon.py

#'
#' Dcp2Cone canonicalizer for the p norm atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a pnorm atom, where
#' the objective is a variable t of dimension of the original
#' vector in the problem and the constraints consist of geometric
#' mean constraints.
Dcp2Cone.pnorm_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@p
  axis <- expr@axis
  expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)

  if(p == 2) {
    if(is.na(axis)) {
      # if(!is.null(expr_dim))
      #  stop("Dimensions should be NULL")
      if(!all(expr_dim == c(1,1)))
        stop("Dimensions should be c(1,1)")
      return(list(t, list(SOC(t, vec(x)))))
    } else
      return(list(t, list(SOC(vec(t), x, axis))))
  }

  # We need an absolute value constraint for the symmetric convex branches (p > 1)
  constraints <- list()
  if(p > 1) {
    # TODO: Express this more naturally (recursively) in terms of the other atoms
    abs_expr <- abs(x)
    canon <- EliminatePwl.abs_canon(abs_expr, abs_expr@args)
    x <- canon[[1]]
    abs_constraints <- canon[[2]]
    constraints <- c(constraints, abs_constraints)
  }

  # Now, we take care of the remaining convex and concave branches to create the
  # rational powers. We need a new variable, r, and the constraint sum(r) == t
  # r <- Variable(dim(x))
  r <- new("Variable", dim = dim(x))
  constraints <- c(constraints, list(sum(r) == t))

  # TODO: No need to run gm_constr to form the tree each time.
  # We only need to form the tree once.
  promoted_t <- Constant(matrix(1, nrow = nrow(x), ncol = ncol(x))) %*% t
  p <- gmp::as.bigq(p)
  if(p < 0)
    constraints <- c(constraints, gm_constrs(promoted_t, list(x, r), c(-p/(1-p), 1/(1-p))))
  else if(p > 0 && p < 1)
    constraints <- c(constraints, gm_constrs(r, list(x, promoted_t), c(p, 1-p)))
  else if(p > 1)
    constraints <- c(constraints, gm_constrs(x, list(r, promoted_t), c(1/p, 1-1/p)))
  return(list(t, constraints))
}
