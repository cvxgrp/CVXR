#' Simple object system that will be sufficient to backing C++ classes for us.
#' @param ... named args to object 
make_LinOp <- function(...) {
  self <- environment()
  class(self) <- c("LinOp", class(self))
  self$self <- self  

  ## Create C++ object
  self$xptr <- .Call(paste0("_CVXR_LinOp__new"), PACKAGE = "CVXR")

  ## Define Methods
  self$get_type <- function() {
    index <- .Call("_CVXR_LinOp__get_type", self$xptr, PACKAGE = "CVXR")
    ## make 1-based index
    LINOP_TYPES[index + 1]
  }
  
  self$is_constant <- function() {
    index <- .Call("_CVXR_LinOp__get_type", self$xptr, PACKAGE = "CVXR")
    ## make 1-based index
    LINOP_TYPES[index + 1] %in% c("SCALAR_CONST", "DENSE_CONST", "SPARSE_CONST")
  }

  self$get_shape <- function() {
    .Call("_CVXR_LinOp__get_shape", self$xptr, PACKAGE = "CVXR")
  }

  self$get_args <- function() {
    .Call("_CVXR_LinOp__get_args", self$xptr, PACKAGE = "CVXR")
    ## TODO The args returned might have to be converted to R structure    
  }
  
  self$get_slice <- function() {
    ## value is a list of integer vectors
    .Call("_CVXR_LinOp__get_slice", self$xptr, PACKAGE = "CVXR")
  }

  self$push_back_slice <- function(integer_vector) {
    .Call("_CVXR_LinOp__push_back_slice", self$xptr, integer_vector, PACKAGE = "CVXR")
  }

  self$has_numerical_data <- function() {
    .Call("_CVXR_LinOp__data_has_been_set", self$xptr, PACKAGE = "CVXR")
  }

  self$get_linOp_data <- function() {
    .Call("_CVXR_LinOp__get_linOp_data", self$xptr, PACKAGE = "CVXR")
    ## TODO: The returned value will have to be converted to R structure
  }

  self$set_linOp_data <- function(linOp_tree) {
    ## TODO: The LinOp_tree will have to be converted to a C++ structure
    .Call("_CVXR_LinOp__set_linOp_data", self$xptr, linOp_tree, PACKAGE = "CVXR")
  }

  self$get_data_ndim <- function() {
    .Call("_CVXR_LinOp__get_data_ndim", self$xptr, PACKAGE = "CVXR")
  }

  self$set_data_ndim <- function(value) { ## int value
    .Call("_CVXR_LinOp__set_data_ndim", self$xptr, value, PACKAGE = "CVXR")
  }
  
  self$is_sparse <- function() {
    ## value should be a boolean
    .Call("_CVXR_LinOp__is_sparse",  self$xptr, PACKAGE = "CVXR")
  }

  self$get_sparse_data <- function() {
    ## In our C++ code Matrix is typedef'd as Eigen::SparseMatrix<double> so a dgCMatrix!
    ## Return value will be a dgCMatrix
    .Call("_CVXR_LinOp__get_sparse_data", self$xptr, PACKAGE = "CVXR")
  }

  self$get_dense_data <- function() {
    ## return value will be a plain matrix
    .Call("_CVXR_LinOp__get_sparse_data", self$xptr, PACKAGE = "CVXR")
  }
  
  self$set_sparse_data <- function(data, row_idx, col_idx, nrow, ncol) {
    ## Input is triplet form components
    ## In our C++ code Matrix is typedef'd as Eigen::SparseMatrix<double> so a dgCMatrix!
    .Call("_CVXR_LinOp__set_sparse_data", self$xptr, data, row_idx, col_idx, nrow, ncol, PACKAGE = "CVXR")
  }
  
  self$set_dense_data <- function(matrix) {
    .Call("_CVXR_LinOp__set_dense_data", self$xptr, matrix, PACKAGE = "CVXR")
  }
  
  self
}

#' @method print LinOp
print.LinOp <- function(x, ...) {
  comps <- ls(x)
  methods <- sapply(comps, function(name) is.function(x[[name]]))
  others <- comps[!methods]
  out <- "LinOp class"
  if (length(others) > 0) {
    out <- c(
      out,
      " " = paste0("Slots: ", paste(others, collapse = ", "))
    )
  }
  if (length(methods) > 0) {
    out <- c(
      out,
      " " = paste0("Methods: ", paste(comps[methods], collapse = ", "))
    )
  }
  cli::cli_bullets(out)
}

