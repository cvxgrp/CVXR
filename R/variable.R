#
# Upper Triangle to Full Matrix
#
# Returns a coefficient matrix to create a symmetric matrix.
#
# @param n The width/height of the matrix
# @return The coefficient matrix.
upper_tri_to_full <- function(n) {
  if(n == 0)
    return(sparseMatrix(i = c(), j = c(), dims = c(0, 0)))

  entries <- floor(n*(n+1)/2)

  val_arr <- c()
  row_arr <- c()
  col_arr <- c()
  count <- 1
  for(i in 1:n) {
    for(j in i:n) {
      # Index in the original matrix
      col_arr <- c(col_arr, count)

      # Index in the filled matrix
      row_arr <- c(row_arr, (j-1)*n + i)
      val_arr <- c(val_arr, 1.0)
      if(i != j) {
        # Index in the original matrix
        col_arr <- c(col_arr, count)

        # Index in the filled matrix
        row_arr <- c(row_arr, (i-1)*n + j)
        val_arr <- c(val_arr, 1.0)
      }
      count <- count + 1
    }
  }
  sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(n^2, entries))
}

#'
#' The Variable class.
#'
#' This class represents an optimization variable.
#'
#' @slot dim The dimensions of the variable.
#' @slot name (Optional) A character string representing the name of the variable.
#' @name Variable-class
#' @aliases Variable
#' @rdname Variable-class
.Variable <- setClass("Variable", representation(dim = "NumORNULL", name = "character", variable_with_attributes = "ANY", delta = "numeric", gradient = "numeric", .is_vector = "logical"),
                                  prototype(dim = NULL, name = NA_character_, variable_with_attributes = NULL, delta = NA_real_, gradient = NA_real_, .is_vector = NA),
                                  validity = function(object) {
                                    ## if(!is.null(object@variable_with_attributes))
                                    ##   stop("[Variable: validation] variable_with_attributes is an internal slot and should not be set by the user")
                                    ## if(!is.na(object@.is_vector))
                                    ##   stop("[Variable: validation] .is_vector is an internal slot and should not be set by the user")
                                    return(TRUE)
                                  }, contains = "Leaf")

#' @param rows The number of rows in the variable.
#' @param cols The number of columns in the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @param var_id (Optional) A unique ID for the variable.
#' @param ... (Optional) Additional attribute arguments. See \linkS4class{Leaf} for details.
#' @rdname Variable-class
#' @examples
#' x <- Variable(3, name = "x0") ## 3-int variable
#' y <- Variable(3, 3, name = "y0") # Matrix variable
#' as.character(y)
#' id(y)
#' is_nonneg(x)
#' is_nonpos(x)
#' size(y)
#' name(y)
#' value(y) <- matrix(1:9, nrow = 3)
#' value(y)
#' grad(y)
#' variables(y)
#' canonicalize(y)
#' @export
Variable <- function(rows = 1, cols = 1, name = NA_character_, var_id = NA_integer_, ...) { .Variable(dim = c(rows, cols), name = name, id = var_id, ...) }

setMethod("initialize", "Variable", function(.Object, ..., dim = NULL, name = NA_character_, value = NA_real_, variable_with_attributes = NULL, delta = NA_real_, gradient = NA_real_, .is_vector = NA) {
  if(length(dim) == 0 || is.null(dim)) {  # Force constants to default to c(1,1).
    dim <- c(1,1)
    .Object@.is_vector <- TRUE
  } else if(length(dim) == 1) {  # Treat as a column vector.
    dim <- c(dim,1)
    .Object@.is_vector <- TRUE
  } else if(length(dim) == 2)
    .Object@.is_vector <- FALSE
  else if(length(dim) > 2)   # TODO: Tensors are currently unimplemented.
    stop("Unimplemented")

  # This is handled in Canonical: .Object@id <- ifelse(is.na(id), get_id(), id)
  if(is.na(name))
    .Object@name <- sprintf("%s%d", VAR_PREFIX, .Object@id)
  else if(is.character(name))
    .Object@name <- name
  else
    stop("Variable name must be a string")

  .Object@variable_with_attributes <- variable_with_attributes
  .Object@value <- value
  .Object@delta <- delta
  .Object@gradient <- gradient
  callNextMethod(.Object, ..., dim = dim, value = .Object@value)
})

#' @docType methods
#' @param x A \linkS4class{Variable} object.
#' @describeIn Variable The name of the variable
#' @export
setMethod("name", "Variable", function(x) { x@name })

setMethod("is_constant", "Variable", function(object) { FALSE })

#' @describeIn Variable The sub/super-gradient of the variable represented as a sparse matrix
setMethod("grad", "Variable", function(object) {
  # TODO: Do not assume dim is 2-D.
  len <- size(object)
  if(len == 0)
    G <- sparseMatrix(i = c(), j = c())
  else
    G <- sparseMatrix(i = 1:len, j = 1:len, x = rep(1, len))

  result <- list()
  result[[as.character(id(object))]] <- G
  return(result)
})

#' @docType methods
#' @describeIn Variable Returns itself as a variable
setMethod("variables", "Variable", function(object) { list(object) })

#' @docType methods
#' @describeIn Variable The graph implementation of the object
setMethod("canonicalize", "Variable", function(object) {
  obj <- lu.create_var(dim(object), id(object))
  list(obj, list())
})

#' @docType methods
#' @describeIn Variable Returns TRUE iff variable generated when lowering a variable with attributes
setMethod("attributes_were_lowered", "Variable", function(object) {
  !is.null(object@variable_with_attributes)
})

#' @docType methods
#' @describeIn Variable Returns the variable of provenance
setMethod("set_variable_of_provenance", signature(object = "Variable", variable = "Variable"), function(object, variable) {
  if(is.null(variable@attributes) || length(variable@attributes) == 0)
    stop("variable attributes must be defined")
  object@variable_with_attributes <- variable
  object
})

#' @docType methods
#' @describeIn Variable Returns a variable with attributes from which this variable was generated
setMethod("variable_of_provenance", "Variable", function(object) {
  object@variable_with_attributes
})

setMethod("show", "Variable", function(object) {
  attr_str <- get_attr_str(object)
  if(length(attr_str) > 0)
    print(paste("Variable((", paste(dim(object), collapse = ", "), "), ", attr_str, ")", sep = ""))
  else
    print(paste("Variable(", paste(dim(object), collapse = ", "), ")", sep = ""))
})

#' @describeIn Variable Convert to a character representation
setMethod("as.character", "Variable", function(x) {
  attr_str <- get_attr_str(x)
  if(length(attr_str) > 0)
    cat("Variable((", paste(dim(x), collapse = ", "), "), ", attr_str, ")", sep = "")
  else
    cat("Variable(", paste(dim(x), collapse = ", "), ")", sep = "")
})
