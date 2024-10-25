## CVXPY SOURCE: cvxpy/atoms/affine/partial_transpose.py
#'
#' The PartialTranspose atom.
#' 
#' Assume \eqn{expr = X_1 \otimes \cdots \otimes X_n} is a 2D Kronecker product 
#' composed of \eqn{n = length(dims)} implicit subsystems. 
#' Letting \eqn{k = axis}, this class represents the partial transpose of \eqn{expr}, 
#' with the transpose applied to its \eqn{k}-th implicit subsystem:
#' 
#' \deqn{X_1 \otimes ...\otimes X_k^T \otimes ... \otimes X_n}.
#'
#' @param expr An \linkS4class{Expression} representing a 2D expression of which to take the partial transpose.
#' @param dims A vector of integers encoding the dimensions of each subsystem.
#' @param axis The index of the subsystem to be transposed out from the tensor product that defines \code{expr}.
#' @return The partial transpose of \code{expr}.
PartialTranspose <- function(expr, dims, axis) {
  expr <- as.Constant(expr)
  dims <- as.integer(dims)
  axis <- as.integer(axis)
  
  if(any(dims) <= 0)
    stop("dims must have positive integer entries")
  if(axis <= 0)
    stop("axis must be a positive integer")
  if(ndim(expr) < 2 || dim(expr)[1] != dim(expr)[2])
    stop("Only supports square matrices")
  if(axis <= 0 || axis > length(dims))
    stop("Invalid axis argument, should be between 1 and ", length(dims))
  if(dim(expr)[1] != base::prod(dims))
    stop("Dimension of system doesn't correspond to dimension of subsystems")
  
  term_list <- list()
  for(i in seq_len(dims[axis])) {
    for(j in seq_len(dims[axis]))
      term_list <- c(term_list, PartialTranspose.term(expr, i, j, dims, axis))
  }
  return(AddExpression(arg_groups = term_list))
}

#
# Helper function for PartialTranspose.
#
# Parameters
# -----------
# expr: The 2D expression of which to take the partial transpose.
# i: The term in the partial transpose sum.
# j: The term in the partial transpose sum.
# dims: A vector of integers encoding the dimensions of each subsystem.
# axis: The index of the subsystem to be transposed from the tensor product that defines expr.
#
# (I ⊗ |i><j| ⊗ I) x (I ⊗ |i><j| ⊗ I) for all (i,j)'s
# in the system we want to transpose.
# This function returns the (i,j)-th term in the sum, namely
# (I ⊗ |i><j| ⊗ I) x (I ⊗ |i><j| ⊗ I).
#
PartialTranspose.term <- function(expr, i, j, dims, axis) {
  a <- Matrix(1, 1, 1, sparse = TRUE)
  
  ndims <- length(dims)
  for(i_axis in seq_len(ndims)) {
    dim <- dims[i_axis]
    if(i_axis == axis) {
      v <- sparseMatrix(i, j, x = 1, dims = c(dim, dim))
      a <- kronecker(a, v)
    } else {
      eye_mat <- sparseMatrix(seq_len(dim), seq_len(dim), x = 1)
      a <- kronecker(a, eye_mat)
    }
  }
  return(a %*% expr %*% a)
}
