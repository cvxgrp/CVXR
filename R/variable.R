#'
#' The Variable class.
#'
#' This class represents an optimization variable.
#'
#' @slot dim The dimensions of the variable.
#' @slot name (Optional) A character string representing the name of the variable.
#' @slot var_id (Internal) A unique identification number used internally.
#' @name Variable-class
#' @aliases Variable
#' @rdname Variable-class
.Variable <- setClass("Variable", representation(dim = "NumORNULL", name = "character", var_id = "integer", value = "ConstVal"),
                                  prototype(dim = NULL, name = NA_character_, var_id = NA_integer_, value = NA_real_), 
                                  validity = function(object) {
                                    if(!is.na(object@value))
                                      stop("[Variable: validation] value is an internal slot and should not be set by the user")
                                    return(TRUE)
                                  }, contains = "Leaf")

#' @param dim The dimensions of the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @rdname Variable-class
#' @examples
#' x <- Variable(3, name = "x0") ## 3-int variable
#' y <- Variable(3, 3, name = "y0") # Matrix variable
#' as.character(y)
#' id(y)
#' is_positive(x)
#' is_negative(x)
#' size(y)
#' name(y)
#' value(y) <- matrix(1:9, nrow = 3)
#' value(y)
#' grad(y)
#' variables(y)
#' canonicalize(y)
#' @export
Variable <- function(dim = NULL, name = NA_character_, var_id = NA_integer_, ...) { .Variable(dim = dim, name = name, var_id = var_id, ...) }

setMethod("initialize", "Variable", function(.Object, ..., dim = NULL, name = NA_character_, var_id = get_id(), value = NA_real_) {
  .Object@var_id <- var_id
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
    paste("Variable((", paste(dim(object), collapse = ", "), "), ", attr_str, ")", sep = "")
  else
    paste("Variable(", paste(dim(object), collapse = ", "), ")", sep = "")
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

#' @describeIn Variable The unique ID of the variable.
setMethod("id", "Variable", function(object) { object@id })

#' @describeIn Variable The name of the variable.
#' @export
setMethod("name", "Variable", function(x) { x@name })

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

# Deprecated constructors
Bool <- function(...) { Variable(..., boolean = TRUE) }
Int <- function(...) { Variable(..., integer = TRUE) }
NonNegative <- function(...) { Variable(..., nonneg = TRUE) }
Symmetric <- function(...) { Variable(..., symmetric = TRUE) }
Semidef <- function(...) { Variable(..., semidef = TRUE) }
