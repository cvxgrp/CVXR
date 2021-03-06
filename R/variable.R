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
.Variable <- setClass("Variable", representation(dim = "NumORNULL", name = "character", value = "ConstVal", .is_vector = "logical"),
                                  prototype(dim = NULL, name = NA_character_, value = NA_real_, .is_vector = NA), 
                                  validity = function(object) {
                                    if(!is.na(object@value))
                                      stop("[Variable: validation] value is an internal slot and should not be set by the user")
                                    return(TRUE)
                                  }, contains = "Leaf")

#' @param rows The number of rows in the variable.
#' @param cols The number of columns in the variable.
#' @param name (Optional) A character string representing the name of the variable.
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
# Variable <- function(dim = NULL, name = NA_character_, id = NA_integer_, ...) { .Variable(dim = dim, name = name, id = id, ...) }
# Variable <- function(rows = 1, cols = 1, name = NA_character_, ...) { .Variable(dim = c(rows, cols), name = name, ...) }
Variable <- function(rows = NULL, cols = NULL, name = NA_character_, ...) { .Variable(dim = c(rows, cols), name = name, ...) }

setMethod("initialize", "Variable", function(.Object, ..., dim = NULL, name = NA_character_, value = NA_real_, .is_vector = NA) {
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
  
  # .Object@id <- ifelse(is.na(id), get_id(), id)
  if(is.na(name))
    .Object@name <- sprintf("%s%d", VAR_PREFIX, .Object@id)
  else
    .Object@name <- name
  .Object@value <- NA_real_
  callNextMethod(.Object, ..., dim = dim)
})

setMethod("show", "Variable", function(object) {
  attr_str <- get_attr_str(object)
  if(length(attr_str) > 0)
    print(paste("Variable((", paste(dim(object), collapse = ", "), "), ", attr_str, ")", sep = ""))
  else
    print(paste("Variable(", paste(dim(object), collapse = ", "), ")", sep = ""))
})

#' @param x,object A \linkS4class{Variable} object.
#' @rdname Variable-class
setMethod("as.character", "Variable", function(x) {
  attr_str <- get_attr_str(x)
  if(length(attr_str) > 0)
    paste("Variable((", paste(dim(x), collapse = ", "), "), ", attr_str, ")", sep = "")
  else
    paste("Variable(", paste(dim(x), collapse = ", "), ")", sep = "")
})

#' @describeIn Variable The name of the variable.
#' @export
setMethod("name", "Variable", function(x) { x@name })

#' @describeIn Variable Get the value of the variable.
setMethod("value", "Variable", function(object) {
  # if(object@.is_vector)
  #  return(as.vector(object@value))
  return(object@value)
})

#' @describeIn Variable The sub/super-gradient of the variable represented as a sparse matrix.
setMethod("grad", "Variable", function(object) {
  # TODO: Do not assume dim is 2-D.
  len <- size(object)
  result <- list(sparseMatrix(i = 1:len, j = 1:len, x = rep(1, len)))
  names(result) <- as.character(id(object))
  result
})

#' @describeIn Variable Returns itself as a variable.
setMethod("variables", "Variable", function(object) { list(object) })

#' @describeIn Variable The canonical form of the variable.
setMethod("canonicalize", "Variable", function(object) {
  obj <- create_var(dim(object), id(object))
  list(obj, list())
})

# This differs from value(object) <- value because there is no validation.
setMethod("save_value", "Variable", function(object, value) {
  object@value <- value
  object
})

# Deprecated constructors
Bool <- function(rows = 1, cols = 1, name = NA_character_) { Variable(rows = rows, cols = cols, name = name, boolean = TRUE) }
Int <- function(rows = 1, cols = 1, name = NA_character_) { Variable(rows = rows, cols = cols, name = name, integer = TRUE) }
NonNegative <- function(rows = 1, cols = 1, name = NA_character_) { Variable(rows = rows, cols = cols, name = name, nonneg = TRUE) }
Symmetric <- function(n, name = NA_character_) { Variable(rows = n, cols = n, name = name, symmetric = TRUE) }
Semidef <- function(n, name = NA_character_) { Variable(rows = n, cols = n, name = name, PSD = TRUE) }

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
