#' Make a CSC sparse diagonal matrix
#' @param size number
#' @return a dgCMatrix diagonal matrix
make_sparse_diagonal_matrix <- function(osize) {
  indices <- seq_len(size)
  ind <- cbind(indices, indices)
  Matrix::sparseMatrix(p = c(0L, indices),
                       i = indices - 1L,
                       x = rep(1.0, size), dims = c(size, size), index1 = FALSE)
}
