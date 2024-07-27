## CVXPY SOURCE: cvxpy/utilities/key_utils.py
#################
#               #
# Key utilities #
#               #
#################
Key <- function(row, col) {
  if(missing(row)) row <- NULL   # Missing row/col index implies that we select all rows/cols
    if(missing(col)) col <- NULL
    list(row = row, col = col, class = "key")
}

ku_validate_key <- function(key, dim) {   # TODO: This may need to be reassessed for consistency in handling keys.
  # Check if the key is a valid index.
  if(length(key) > 3)
    stop("Invalid index/slice")

  nrow <- dim[1]
  ncol <- dim[2]
  row <- ku_format_slice(key$row, nrow)
  col <- ku_format_slice(key$col, ncol)

  # Change single indices for vectors into double indices
  if(!is.null(row) && !is.null(col))
    key <- Key(row = row, col = col)
  else if(is.null(row) && !is.null(col))
    key <- Key(row = seq_len(nrow), col = col)
  else if(!is.null(row) && is.null(col))
    key <- Key(row = row, col = seq_len(ncol))
  else
    stop("A key cannot be empty")
  return(key)
}

ku_format_slice <- function(key_val, dim, axis) {
  # Converts part of a key into a slice with a start and step.
  if(is.null(key_val))
    return(NULL)
  orig_key_val <- as.integer(key_val)

  # Return if all zero indices.
  if(all(orig_key_val == 0))
    return(orig_key_val)

  # Convert negative indices to positive indices.
  if(all(orig_key_val >= 0))
    key_val <- orig_key_val
  else if(all(orig_key_val <= 0))
    key_val <- setdiff(seq_len(dim), -orig_key_val)
  else
    stop("Only 0's may be mixed with negative subscripts")

  if(all(key_val >= 0 & key_val <= dim))
    return(key_val)
  else
    stop("Index is out of bounds for axis with size ", dim)
}

ku_slice_mat <- function(mat, key) {
  if(is.vector(mat))
    mat <- matrix(mat, ncol = 1)

  if(is.matrix(key$row) && is.null(key$col))
    select_mat  <- matrix(mat[key$row], ncol = 1)
  else if(is.null(key$row) && is.null(key$col))
    select_mat <- mat
  else if(is.null(key$row) && !is.null(key$col))
    select_mat <- mat[, key$col, drop = FALSE]
  else if(!is.null(key$row) && is.null(key$col))
    select_mat <- mat[key$row, , drop = FALSE]
  else
    select_mat <- mat[key$row, key$col, drop = FALSE]
  select_mat
}

ku_dim <- function(key, dim) {
  dims <- c()

  for(i in 1:2) {
    idx <- key[[i]]
    if(is.null(idx))
      size <- dim[i]
    else {
      selection <- (1:dim[i])[idx]
      size <- length(selection)
    }
    dims <- c(dims, size)
  }
  dims
}

