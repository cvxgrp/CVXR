# Get the dimensions of the constant.
intf_dim <- function(constant) {
  if((is.null(dim(constant)) && length(constant) == 1) ||
     (!is.null(dim(constant)) && all(dim(constant) == c(1,1))))
    return(c(1,1))
  else if(is.vector(constant) && !is.list(constant))
    return(c(1,length(constant)))
  else if(is.matrix(constant) || is(constant, "Matrix"))
    return(dim(constant))
  else
    stop("Unknown class: ", class(constant))
}

# Is the constant a column vector?
intf_is_vector <- function(constant) {
  dim(constant)[2] == 1
}

# Is the constant a scalar?
intf_is_scalar <- function(constant) {
  all(dim(constant) == c(1,1))
}

# Return the collective sign of the matrix entries.
intf_sign <- function(constant) {
  c(min(constant) >= 0, max(constant) <= 0)
}

# Return (is real, is imaginary).
intf_is_complex <- function(constant, tol = 1e-5) {
  real_max <- max(abs(Re(constant)))
  imag_max <- max(abs(Im(constant)))
  c(real_max >= tol, imag_max >= tol)
}

# Check if a matrix is Hermitian and/or symmetric
intf_is_hermitian <- function(constant) {
  # I'm replicating np.allclose here to match CVXPY. May want to use base::isSymmetric later.
  # is_herm <- isSymmetric(constant)   # This returns TRUE iff a complex matrix is Hermitian.
  # TODO: Catch complex symmetric, but not Hermitian?
  atol <- 1e-8
  rtol <- 1e-5
  is_symm <- all(abs(constant - t(constant)) <= atol + rtol*abs(t(constant)))
  is_herm <- all(abs(constant - Conj(t(constant))) <= atol + rtol*Conj(t(constant)))
  c(is_symm, is_herm)
}

intf_scalar_value <- function(constant) {
  if(is.list(constant))
    constant[[1]]
  else
    as.numeric(constant)
}

intf_convert_if_scalar <- function(constant) {
  if(all(intf_dim(constant) == c(1,1)))
    intf_scalar_value(constant)
  else
    constant
}

intf_block_add <- function(mat, block, vert_offset, horiz_offset, rows, cols, vert_step = 1, horiz_step = 1) {
  block <- matrix(block, nrow = rows, ncol = cols)
  vert_idx <- seq(vert_offset + 1, vert_offset + rows, vert_step)
  horiz_idx <- seq(horiz_offset + 1, horiz_offset + cols, horiz_step)
  mat[vert_idx, horiz_idx] <- mat[vert_idx, horiz_idx] + block
  mat
}
