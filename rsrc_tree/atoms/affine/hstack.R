## CVXPY SOURCE: cvxpy/atoms/affine/HStack.py
#'
#' The HStack class.
#'
#' Horizontal concatenation of values.
#'
#' @slot ... \linkS4class{Expression} objects or matrices. All arguments must have the same dimensions except for axis 2 (columns).
#' @name HStack-class
#' @aliases HStack
#' @rdname HStack-class
.HStack <- setClass("HStack", contains = "AffAtom")

#' @param ... \linkS4class{Expression} objects or matrices. All arguments must have the same dimensions except for axis 2 (columns).
#' @rdname HStack-class
HStack <- function(...) {
  arg_list <- lapply(list(...), as.Constant)
  for(idx in seq_along(arg_list)) {
    arg <- arg_list[[idx]]
    if(ndim(arg) == 0)
      arg_list[[idx]] <- flatten(arg)
  }
  .HStack(atom_args = arg_list)
}

#' @param object A \linkS4class{HStack} object.
#' @param values A list of arguments to the atom.
#' @describeIn HStack Horizontally concatenate the values using \code{cbind}.
setMethod("to_numeric", "HStack", function(object, values) {
  # do.call("cbind", values)   # Doesn't work on some objects like xts.
  mat <- Reduce("cbind", values)
  colnames(mat) <- NULL   # Get rid of init column name.
  return(mat)
})

#' @describeIn HStack The dimensions of the atom.
setMethod("dim_from_args", "HStack", function(object) {
  if(ndim(object@args[[1]]) == 1)
    # return(c(sum(sapply(object@args, size)), NA))
    return(c(sum(sapply(object@args, size)), 1))
  else {
    cols <- sum(sapply(object@args, function(arg) { dim(arg)[2] }))
    arg_dim <- dim(object@args[[1]])
    dims <- c(arg_dim[1], cols)
    if(length(arg_dim) >= 3)
      dims <- c(dims, arg_dim[3:length(arg_dim)])
    return(dims)
  }
})

#' @describeIn HStack Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "HStack", function(object) { TRUE })

#' @describeIn HStack Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "HStack", function(object) { TRUE })

#' @describeIn HStack Check all arguments have the same height.
setMethod("validate_args", "HStack", function(object) {
  model <- dim(object@args[[1]])
  error <- "All the input dimensions except for axis 2 (columns) must match exactly."
  len <- length(object@args)
  
  if(len >= 2) {
    for(arg in object@args[2:len]) {
      if(length(dim(arg)) != length(model))
        stop(error)
      else if(length(model) > 1) {
        for(i in 1:length(model)) {
          if(i != 2 && dim(arg)[i] != model[i])
            stop(error)
        }
      }
    }
  }
})

HStack.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lu.hstack(arg_objs, dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn HStack The graph implementation of the atom.
setMethod("graph_implementation", "HStack", function(object, arg_objs, dim, data = NA_real_) {
  HStack.graph_implementation(arg_objs, dim, data)
})

setMethod("cbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { HStack(x, y) })
setMethod("cbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { HStack(x, y) })

