## CVXPY SOURCE: cvxpy/constraints/utilities.py
#####################################
#                                   #
# Utility functions for constraints #
#                                   #
#####################################
# Formats all the row/column cones for the solver.
format_axis <- function(t, X, axis) {
  # Reduce to norms of columns
  if(axis == 1)
    X <- lu.transpose(X)

  # Create matrices Tmat, Xmat such that Tmat*t + Xmat*X
  # gives the format for the elementwise cone constraints.
  cone_size <- 1 + nrow(X)
  terms <- list()

  # Make t_mat.
  mat_dim <- c(cone_size, 1)
  t_mat <- sparseMatrix(i = 1, j = 1, x = 1.0, dims = mat_dim)
  t_mat <- lu.create_const(t_mat, mat_dim, sparse = TRUE)
  t_vec <- t
  if(is.null(dim(t)))   # t is scalar.
    t_vec <- lu.reshape(t, c(1,1))
  else   # t is 1-D.
    t_vec <- lu.reshape(t, c(1, nrow(t)))
  mul_dim <- c(cone_size, ncol(t_vec))
  terms <- c(terms, list(lu.mul_expr(t_mat, t_vec, mul_dim)))

  # Make X_mat.
  if(length(dim(X)) == 1)
    X <- lu.reshape(X, c(nrow(X), 1))
  mat_dim <- c(cone_size, nrow(X))
  val_arr <- rep(1.0, cone_size - 1)
  row_arr <- 2:cone_size
  col_arr <- 1:(cone_size - 1)
  X_mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = mat_dim)
  X_mat <- lu.create_const(X_mat, mat_dim, sparse = TRUE)
  mul_dim <- c(cone_size, ncol(X))
  terms <- c(terms, list(lu.mul_expr(X_mat, X, mul_dim)))
  list(lu.create_geq(lu.sum_expr(terms)))
}

# Formats all the elementwise cones for the solver.
format_elemwise <- function(vars_) {
  # Create matrices A_i such that 0 <= A_0*x_0 + ... + A_n*x_n
  # gives the format for the elementwise cone constraints.
  spacing <- length(vars_)

  # Matrix spaces out columns of the LinOp expressions
  var_dim <- vars_[[1]]$dim
  mat_dim <- c(spacing*var_dim[1], var_dim[1])

  mats <- lapply(0:(spacing-1), function(offset) { get_spacing_matrix(mat_dim, spacing, offset) })
  terms <- mapply(function(var, mat) { list(lu.mul_expr(mat, var)) }, vars_, mats)
  list(lu.create_geq(lu.sum_expr(terms)))
}

# Returns a sparse matrix LinOp that spaces out an expression.
get_spacing_matrix <- function(dim, spacing, offset) {
  col_arr <- seq_len(dim[2])
  row_arr <- spacing * (col_arr - 1) + 1 + offset
  val_arr <- rep(1.0, dim[2])
  mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = dim)
  lu.create_const(mat, dim, sparse = TRUE)
}
