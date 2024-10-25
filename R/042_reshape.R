## CVXPY SOURCE: cvxpy/atoms/affine/reshape.py

#'
#' The Reshape class.
#'
#' This class represents the reshaping of an expression. The operator vectorizes the expression,
#' then unvectorizes it into the new dimensions. Entries are stored in column-major order.
#'
#' @slot expr An \linkS4class{Expression} or numeric matrix.
#' @slot new_dim The new dimensions.
#' @name Reshape-class
#' @aliases Reshape
#' @rdname Reshape-class
.Reshape <- setClass("Reshape", representation(expr = "ConstValORExpr", new_dim = "numeric", byrow = "logical"), 
                                prototype(byrow = FALSE), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} or numeric matrix.
#' @param new_dim The new dimensions.
#' @param byrow A logical value indicating whether the matrix is filled by rows. If \code{FALSE} (default), the matrix is filled by columns.
#' @rdname Reshape-class
Reshape <- function(expr, new_dim, byrow = FALSE) { .Reshape(expr = expr, new_dim = new_dim, byrow = byrow) }

setMethod("initialize", "Reshape", function(.Object, ..., expr, new_dim, byrow) {
  if(length(new_dim) > 2)
    stop("Expressions of dimension greater than 2 are not supported")
  .Object@new_dim <- new_dim
  .Object@expr <- expr
  .Object@byrow <- byrow
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Reshape} object.
#' @param values A list of arguments to the atom.
#' @describeIn Reshape Reshape the value into the specified dimensions.
setMethod("to_numeric", "Reshape", function(object, values) {
  if(object@byrow) {
    mat <- values[[1]]
    dim(mat) <- rev(object@new_dim)
    return(t(mat))
  } else {
    dim(values[[1]]) <- object@new_dim
    return(values[[1]])
  }

})

#' @describeIn Reshape Check the new shape has the same number of entries as the old.
setMethod("validate_args", "Reshape", function(object) {
  old_len <- size(object@args[[1]])
  new_len <- prod(object@new_dim)
  if(old_len != new_len)
    stop("Invalid reshape dimensions (", paste(object@new_dim, sep = ","), ")")
})

#' @describeIn Reshape The \code{c(rows, cols)} dimensions of the new expression.
setMethod("dim_from_args", "Reshape", function(object) { object@new_dim })

#' @describeIn Reshape Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Reshape", function(object) { TRUE })

#' @describeIn Reshape Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Reshape", function(object) { TRUE })

#' @describeIn Reshape Returns a list containing the new shape.
setMethod("get_data", "Reshape", function(object) { list(object@new_dim, object@byrow) })

Reshape.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  arg <- arg_objs[[1]]
  if(data[2]) {   # byrow = TRUE
    arg <- lu.transpose(arg)
    if(length(dim) <= 1)
      list(lu.reshape(arg, dim), list())
    else {
      result <- lu.reshape(arg, c(dim[2], dim[1]))
      list(lu.transpose(result), list())
    }
  } else   # byrow = FALSE
    list(lu.reshape(arg, dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Reshape The graph implementation of the atom.
setMethod("graph_implementation", "Reshape", function(object, arg_objs, dim, data = NA_real_) {
  Reshape.graph_implementation(arg_objs, dim, data)
})

Reshape.deep_flatten <- function(x)
{
  # Base cases
  if(is(x, "Expression")) {
    if(length(dim(x)) == 1 || all(dim(x) == c(1,1)))
      return(x)
    else
      return(flatten(x))
  } else if(is.numeric(x)) {
    x <- as.Constant(x)
    return(flatten(x))
  }
  
  # Recursion
  if(is(x, "list")) {
    y <- list()
    for(xo in x) {
      x1 <- deep_flatten(x0)
      y <- c(y, x1)
    }
    y <- HStack(y)
    return(y)
  }
  
  stop("The input to deep_flatten must be an Expression, numeric array or scalar, or a nested list thereof. Received input of type ", class(x))
}
