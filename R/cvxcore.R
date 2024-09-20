#' This is a class that backs a C++ LinOp class
#' @param type the type of LinOp, one of the types above
#' @param dim the shape of the LinOp, a tuple, so for us a vector of integers
#' @param args the arguments of the LinOp
#' @param data the data for the LinOp, which is later set to C++ LinOp objects' linOp_data_ field
cvxcore_LinOp <- function(type, dim, args, data = NULL) {
  self <- environment()
  class(self) <- c("cvxcore_LinOp", class(self))
  self$self <- self  

  ## Create C++ object
  self$xptr <- .Call(paste0("_CVXR_LinOp__new"), type, dim, args, PACKAGE = "CVXR")

  ## Define Methods
  self$get_type <- function() {
    index <- .Call("_CVXR_LinOp__get_type", self$xptr, PACKAGE = "CVXR")
    ## make 1-based index
    LINOP_TYPES[index + 1L]
  }
  
  self$is_constant <- function() {
    index <- .Call("_CVXR_LinOp__get_type", self$xptr, PACKAGE = "CVXR")
    ## make 1-based index
    LINOP_TYPES[index + 1L] %in% c("SCALAR_CONST", "DENSE_CONST", "SPARSE_CONST")
  }

  self$get_dim <- function() {
    .Call("_CVXR_LinOp__get_dim", self$xptr, PACKAGE = "CVXR")
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

#' @method print cvxcore.LinOp
print.cvxcore_LinOp <- function(x, ...) {
  comps <- ls(x)
  methods <- sapply(comps, function(name) is.function(x[[name]]))
  others <- comps[!methods]
  out <- class(x)[[1L]]
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

## ## cvxcore class that basically exposes the build_matrix C++ function
## ## build_matrix methods
## #' @param ... named args to object 
## make_cvxcore <- function(...) {
##   self <- environment()
##   class(self) <- c("cvxcore", class(self))
##   self$self <- self  

##   ## Create interfaces to C++ functions
##   self$build_matrix <- function(constraints, var_length,
##                                 id_to_col, param_to_size,
##                                 num_threads) {
##     ## constraints is an R object backing a vector of ConstLinOps (see make_ConstLinOpVector() in LinOpVector.R)
##     ## id_to_col is an integer vector with names that are integers converted to chacracters
##     .Call('_CVXR_build_matrix_0',
##           constraints$getXPtr(),
##           id_to_col,
##           PACKAGE = 'CVXR')
##     } else {
##       objPtr <- .Call('_CVXR_build_matrix_1',
##                       constraints$getXPtr(),
##                       id_to_col,
##                       constr_offsets,
##                       PACKAGE = 'CVXR')
##     }
##     ##cat("Instantiating ProblemData-R6", "\n")
##     CVXcanon.ProblemData$new(objPtr)
##   }
##   ))
