#'
#' Given a problem returns a PSD constrain
#'
#' @param problem A \linkS4class{Problem} object.
#' @param c A vector of coefficients.
#' @return Returns an array G and vector h such that the given constraint is
#' equivalent to G*z <=_{PSD} h.
psd_coeff_offset <- function(problem, c) {
  extractor <- CoeffExtractor(InverseData(problem))
  tmp <- affine(extractor, expr(c))
  A_vec <- tmp[[1]]
  b_vec <- tmp[[2]]
  G <- -A_vec
  h <- b_vec
  dim <- nrow(expr(c))
  return(list(G, h, dim))
}

#'
#' Utility methods for special handling of semidefinite constraints.
#'
#' @param matrix The matrix to get the lower triangular matrix for
#' @return The lower triangular part of the matrix, stacked in column-major order
scaled_lower_tri <- function(matrix) {
  # Returns an expression representing the lower triangular entries.
  # Scales the strictly lower triangular entries by sqrt(2), as
  # required by SCS.
  rows <- cols <- nrow(matrix)
  entries <- floor(rows * (cols + 1)/2)

  row_arr <- seq_len(entries)

  col_arr <- matrix(1:(rows*cols), nrow = rows, ncol = cols)
  col_arr <- col_arr[lower.tri(col_arr, diag = TRUE)]

  val_arr <- matrix(0, nrow = rows, ncol = cols)
  val_arr[lower.tri(val_arr, diag = TRUE)] <- sqrt(2)
  diag(val_arr) <- 1
  val_arr <- as.vector(val_arr)
  val_arr <- val_arr[val_arr != 0]

  coeff <- Constant(sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(entries, rows*cols)))
  vectorized_matrix <- reshape_expr(matrix, c(rows*cols, 1))
  return(coeff %*% vectorized_matrix)
}

#'
#' Is the constraint a stuffed cone constraint?
#'
#' @param constraint A \linkS4class{Constraint} object.
#' @return Is the constraint a stuffed-cone constraint?
is_stuffed_cone_constraint <- function(constraint) {
  # Conic solvers require constraints to be stuffed in the following way.
  if(length(variables(constraint)) != 1)
    return(FALSE)
  for(arg in constraint@args) {
    if(inherits(arg, "Reshape"))
      arg <- arg@args[[1]]
    if(inherits(arg, "AddExpression")) {
      if(!inherits(arg@args[[1]], c("MulExpression", "Multiply")))
        return(FALSE)
      if(!inherits(arg@args[[1]]@args[[1]], "Constant"))
        return(FALSE)
      if(!inherits(arg@args[[2]], "Constant"))
        return(FALSE)
    } else if(inherits(arg, c("MulExpression", "Multiply"))) {
      if(!inherits(arg@args[[1]], "Constant"))
        return(FALSE)
    } else
      return(FALSE)
  }
  return(TRUE)
}

#'
#' Is the objective a stuffed cone objective?
#'
#' @param objective An \linkS4class{Objective} object.
#' @return Is the objective a stuffed-cone objective?
is_stuffed_cone_objective <- function(objective) {
  # Conic solvers require objectives to be stuffed in the following way.
  expr <- expr(objective)
  return(is_affine(expr) && length(variables(expr)) == 1 && inherits(expr, "AddExpression") && length(expr@args) == 2
                         && inherits(expr@args[[1]], c("MulExpression", "Multiply")) && inherits(expr@args[[2]], "Constant"))
}

