## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/log_det_canon.py

#'
#' Dcp2Cone canonicalizer for the log determinant atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a log determinant atom where
#' the objective function is the sum of the log of the vector D
#' and the constraints consist of requiring the matrix Z to be
#' diagonal and the diagonal Z to equal D, Z to be upper triangular
#' and DZ; t(Z)A to be positive semidefinite, where A is a n by n
#' matrix.
Dcp2Cone.log_det_canon <- function(expr, args) {
  # Reduces the atom to an affine expression and list of constraints.
  #
  # Creates the equivalent problem::
  #
  # maximize    sum(log(D[i, i]))
  # subject to: D diagonal
  # diag(D) = diag(Z)
  # Z is upper triangular.
  # [D Z; t(Z) A] is positive semidefinite
  #
  # The problem computes the LDL factorization:
  #
  # A = (Z^TD^{-1})D(D^{-1}Z)
  #
  # This follows from the inequality:
  #
  # \det(A) >= \det(D) + \det([D, Z; Z^T, A])/\det(D) >= \det(D)
  #
  # because (Z^TD^{-1})D(D^{-1}Z) is a feasible D, Z that achieves
  # det(A) = det(D) and the objective maximizes det(D).
  #
  # Parameters
  # ----------
  # expr : log_det
  # args : list of arguments for the expression
  #
  # Returns
  # -------
  # (Variable for objective, list of constraints)

  A <- args[[1]]  # n by n matrix.
  n <- nrow(A)

  z_small <- Variable(floor(n*(n+1)/2))
  Z <- UpperTri.vec_to_upper_tri(z_small, strict = FALSE)
  d_small <- DiagMat(Z)  # a vector
  D <- DiagVec(d_small)  # a matrix
  X <- bmat(list(list(D, Z),
                 list(t(Z), A)))

  constraints <- list(PSDConstraint(X))
  log_expr <- log(d)
  canon <- Dcp2Cone.log_canon(log_expr, log_expr@args)
  obj <- canon[[1]]
  constr <- canon[[2]]
  constraints <- c(constraints, constr)
  return(list(sum(obj), constraints))
}
