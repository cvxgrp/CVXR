#' Make a CSC sparse diagonal matrix
#' @param size number of rows or columns
#' @param diagonal if specified, the diagonal values, in which case size is ignored
#' @return a compressed sparse column diagonal matrix
make_sparse_diagonal_matrix <- function(size, diagonal = NULL) {
  if (!is.null(diagonal)) {
    size <- length(diagonal)
    values <- diagonal
  } else {
    values <- rep(1.0, size)
  }
  indices <- seq_len(size)
  Matrix::sparseMatrix(p = c(0L, indices),
                       i = indices - 1L,
                       x = values,
                       dims = c(size, size), index1 = FALSE)
}
