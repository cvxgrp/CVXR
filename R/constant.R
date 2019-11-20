#'
#' The Constant class.
#'
#' This class represents a constant.
#'
#' @slot value A numeric element, vector, matrix, or data.frame. Vectors are automatically cast into a matrix column.
#' @slot sparse (Internal) A logical value indicating whether the value is a sparse matrix.
#' @slot is_pos (Internal) A logical value indicating whether all elements are non-negative.
#' @slot is_neg (Internal) A logical value indicating whether all elements are non-positive.
#' @name Constant-class
#' @aliases Constant
#' @rdname Constant-class
.Constant <- setClass("Constant", representation(value = "ConstVal", sparse = "logical", imag = "logical", nonneg = "logical", nonpos = "logical", symm = "logical", herm = "logical", eigvals = "numeric", .is_vector = "logical", .cached_is_pos = "logical"),
                                  prototype(value = NA_real_, sparse = NA, imag = NA, nonneg = NA, nonpos = NA, symm = NA, herm = NA, eigvals = NA_real_, .is_vector = NA, .cached_is_pos = NA), contains = "Leaf")

#' @param value A numeric element, vector, matrix, or data.frame. Vectors are automatically cast into a matrix column.
#' @rdname Constant-class
#' @examples
#' x <- Constant(5)
#' y <- Constant(diag(3))
#' get_data(y)
#' value(y)
#' is_nonneg(y)
#' size(y)
#' as.Constant(y)
#' @export
Constant <- function(value) { .Constant(value = value) }

setMethod("initialize", "Constant", function(.Object, ..., value = NA_real_, sparse = NA, imag = NA, nonneg = NA, nonpos = NA, symm = NA, herm = NA, eigvals = NA_real_, .is_vector = NA, .cached_is_pos = NA) {
  # Keep sparse matrices sparse.
  if(is(value, "ConstSparseVal")) {
    .Object@value <- Matrix(value, sparse = TRUE)
    .Object@sparse <- TRUE
  } else {
    .Object@value <- as.matrix(value)
    .Object@sparse <- FALSE
  }
  .Object@imag <- imag
  .Object@nonneg <- nonneg
  .Object@nonpos <- nonpos
  .Object@symm <- symm
  .Object@herm <- herm
  .Object@eigvals <- eigvals
  .Object@.is_vector <- is.vector(value)
  .Object@.cached_is_pos <- .cached_is_pos
  callNextMethod(.Object, ..., dim = intf_dim(.Object@value))
})

#' @param x,object A \linkS4class{Constant} object.
#' @rdname Constant-class
setMethod("show", "Constant", function(object) {
  cat("Constant(", curvature(object), ", ", sign(object), ", (", paste(dim(object), collapse = ","), "))", sep = "")
})

#' @describeIn Constant The name of the constant.
setMethod("name", "Constant", function(x) { as.character(head(x@value)) })

#' @describeIn Constant Returns itself as a constant.
setMethod("constants", "Constant", function(object) { list(object) })

#' @describeIn Constant The value of the constant.
setMethod("value", "Constant", function(object) {
  # if(object@.is_vector)
  #  return(as.vector(object@value))
  return(object@value)
})

setMethod("is_pos", "Constant", function(object) {
  if(is.na(object@.cached_is_pos))
    object@.cached_is_pos <- all(object@value > 0)
  object@.cached_is_pos
})

#' @describeIn Constant An empty list since the gradient of a constant is zero.
setMethod("grad", "Constant", function(object) { list() })

#' @describeIn Constant The \code{c(row, col)} dimensions of the constant.
setMethod("dim", "Constant", function(x) { x@dim })

#' @describeIn Constant The canonical form of the constant.
setMethod("canonicalize", "Constant", function(object) {
  obj <- create_const(value(object), dim(object), object@sparse)
  list(obj, list())
})

#' @describeIn Constant A logical value indicating whether all elements of the constant are non-negative.
setMethod("is_nonneg", "Constant", function(object) { 
  if(is.na(object@nonneg))
    object <- .compute_attr(object)
  object@nonneg
})

#' @describeIn Constant A logical value indicating whether all elements of the constant are non-positive.
setMethod("is_nonpos", "Constant", function(object) {
  if(is.na(object@nonpos))
    object <- .compute_attr(object)
  object@nonpos
})

#' @describeIn Constant A logical value indicating whether the constant is imaginary.
setMethod("is_imag", "Constant", function(object) {
  if(is.na(object@imag))
    object <- .compute_attr(object)
  object@imag
})

#' @describeIn Constant A logical value indicating whether the constant is complex-valued.
setMethod("is_complex", "Constant", function(object) { is.complex(value(object)) })

#' @describeIn Constant A logical value indicating whether the constant is symmetric.
setMethod("is_symmetric", "Constant", function(object) {
  if(is_scalar(object))
    return(TRUE)
  else if(ndim(object) == 2 && nrow(object) == ncol(object)) {
    if(is.na(object@symm))
      object <- .compute_symm_attr(object)
    return(object@symm)
  } else
    return(FALSE)
})

#' @describeIn Constant A logical value indicating whether the constant is a Hermitian matrix.
setMethod("is_hermitian", "Constant", function(object) {
  if(is_scalar(object) && is_real(object))
    return(TRUE)
  else if(ndim(object) == 2 && nrow(object) == ncol(object)) {
    if(is.na(object@herm))
      object <- .compute_symm_attr(object)
    return(object@herm)
  } else
    return(FALSE)
})

# Compute the attributes of the constant related to complex/real, sign.
.compute_attr <- function(object) {
  # Set DCP attributes.
  res <- intf_is_complex(value(object))
  is_real <- res[[1]]
  is_imag <- res[[2]]
  
  if(is_complex(object)) {
    is_nonneg <- FALSE
    is_nonpos <- FALSE
  } else {
    sign <- intf_sign(value(object))
    is_nonneg <- sign[[1]]
    is_nonpos <- sign[[2]]
  }
  
  object@imag <- is_imag && !is_real
  object@nonpos <- is_nonpos
  object@nonneg <- is_nonneg
  object
}

# Determine whether the constant is symmetric/Hermitian.
.compute_symm_attr <- function(object) {
  # Set DCP attributes.
  res <- intf_is_hermitian(value(object))
  object@symm <- res[[1]]
  object@herm <- res[[2]]
  object
}

# Compute the eigenvalues of the Hermitian or symmetric matrix represented by this constant.
.compute_eigvals <- function(object) {
  object@eigvals <- eigen(value(object), only.values = TRUE)$values
  object
}

#' @describeIn Constant A logical value indicating whether the constant is a positive semidefinite matrix.
setMethod("is_psd", "Constant", function(object) {
  # Symbolic only cases.
  if(is_scalar(object) && is_nonneg(object))
    return(TRUE)
  else if(is_scalar(object))
    return(FALSE)
  else if(ndim(object) == 1)
    return(FALSE)
  else if(ndim(object) == 2 && nrow(object) != ncol(object))
    return(FALSE)
  else if(!is_hermitian(object))
    return(FALSE)
  
  # Compute eigenvalues if absent.
  if(is.na(object@eigvals))
    object <- .compute_eigvals(object)
  return(all(Re(object@eigvals) >= -EIGVAL_TOL))
})

#' @describeIn Constant A logical value indicating whether the constant is a negative semidefinite matrix.
setMethod("is_nsd", "Constant", function(object) {
  # Symbolic only cases.
  if(is_scalar(object) && is_nonpos(object))
    return(TRUE)
  else if(is_scalar(object))
    return(FALSE)
  else if(ndim(object) == 1)
    return(FALSE)
  else if(ndim(object) == 2 && nrow(object) != ncol(object))
    return(FALSE)
  else if(!is_hermitian(object))
    return(FALSE)
  
  # Compute eigenvalues if absent.
  if(is.na(object@eigvals))
    object <- .compute_eigvals(object)
  return(all(Re(object@eigvals) <= EIGVAL_TOL))
})

#'
#' Cast to a Constant
#'
#' Coerce an R object or expression into the \linkS4class{Constant} class.
#'
#' @param expr An \linkS4class{Expression}, numeric element, vector, matrix, or data.frame.
#' @return A \linkS4class{Constant} representing the input as a constant.
#' @docType methods
#' @rdname Constant-class
#' @export
as.Constant <- function(expr) {
  if(is(expr, "Expression"))
    expr
  else
    Constant(value = expr)
}

#'
#' The Parameter class.
#'
#' This class represents a parameter, either scalar or a matrix.
#'
#' @slot rows The number of rows in the parameter.
#' @slot cols The number of columns in the parameter.
#' @slot name (Optional) A character string representing the name of the parameter.
#' @slot value (Optional) A numeric element, vector, matrix, or data.frame. Defaults to \code{NA} and may be changed with \code{value<-} later.
#' @name Parameter-class
#' @aliases Parameter
#' @rdname Parameter-class
.Parameter <- setClass("Parameter", representation(dim = "numeric", name = "character", value = "ConstVal", .is_vector = "logical"),
                                    prototype(dim = NULL, name = NA_character_, value = NA_real_, .is_vector = NA), contains = "Leaf")

#' @param rows The number of rows in the parameter.
#' @param cols The number of columns in the parameter.
#' @param name (Optional) A character string representing the name of the parameter.
#' @param value (Optional) A numeric element, vector, matrix, or data.frame. Defaults to \code{NA} and may be changed with \code{value<-} later.
#' @rdname Parameter-class
#' @examples
#' x <- Parameter(3, name = "x0", sign="NONPOSITIVE") ## 3-vec negative
#' is_nonneg(x)
#' is_nonpos(x)
#' size(x)
#' @export
# Parameter <- function(dim = NULL, name = NA_character_, value = NA_real_, ...) { .Parameter(dim = dim, name = name, value = value, ...) }
# Parameter <- function(rows = 1, cols = 1, name = NA_character_, value = NA_real_, ...) { .Parameter(dim = c(rows, cols), name = name, value = value, ...) }
Parameter <- function(rows = NULL, cols = NULL, name = NA_character_, value = NA_real_, ...) { .Parameter(dim = c(rows, cols), name = name, value = value, ...) }

setMethod("initialize", "Parameter", function(.Object, ..., dim = NULL, name = NA_character_, value = NA_real_, .is_vector = NA) {
  # .Object@id <- get_id()
  if(is.na(name))
    .Object@name <- sprintf("%s%s", PARAM_PREFIX, .Object@id)
  else
    .Object@name <- name
  
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

  # Initialize with value if provided
  # .Object@value <- value
  # callNextMethod(.Object, ..., id = .Object@id, dim = dim, value = value)
  .Object@value <- NA_real_
  callNextMethod(.Object, ..., dim = dim, value = value)
})

#' @param object,x A \linkS4class{Parameter} object.
#' @describeIn Parameter Returns \code{list(dim, name, value, attributes)}.
setMethod("get_data", "Parameter", function(object) {
  list(dim = dim(object), name = object@name, value = value(object), attributes = attributes(object))
})

#' @describeIn Parameter The name of the parameter.
#' @export
setMethod("name", "Parameter", function(x) { x@name })

#' @describeIn Parameter The value of the parameter.
setMethod("value", "Parameter", function(object) {
  # if(object@.is_vector)
  #  return(as.vector(object@value))
  return(object@value)
})

#' @describeIn Parameter Set the value of the parameter.
setReplaceMethod("value", "Parameter", function(object, value) {
  object@value <- validate_val(object, value)
  object
})

#' @describeIn Parameter An empty list since the gradient of a parameter is zero.
setMethod("grad", "Parameter", function(object) { list() })

#' @describeIn Parameter Returns itself as a parameter.
setMethod("parameters", "Parameter", function(object) { list(object) })

#' @describeIn Parameter The canonical form of the parameter.
setMethod("canonicalize", "Parameter", function(object) {
  obj <- create_param(object, dim(object))
  list(obj, list())
})

setMethod("show", "Parameter", function(object) {
  attr_str <- get_attr_str(object)
  if(length(attr_str) > 0)
    cat("Parameter(", paste(dim(object), collapse = ", "), ", ", attr_str, ")", sep = "")
  else
    cat("Parameter(", paste(dim(object), collapse = ", "), ")", sep = "")
})

#'
#' The CallbackParam class.
#'
#' This class represents a parameter whose value is obtained by evaluating a function.
#'
#' @slot callback A numeric element, vector, matrix, or data.frame.
#' @name CallbackParam-class
#' @aliases CallbackParam
#' @rdname CallbackParam-class
.CallbackParam <- setClass("CallbackParam", representation(callback = "function", dim = "numeric"), 
                                            prototype(dim = NULL), contains = "Parameter")

#' @param callback A callback function that generates the parameter value.
#' @param dim The dimensions of the parameter.
#' @param name (Optional) A character string representing the name of the parameter.
#' @param sign A character string indicating the sign of the parameter. Must be "ZERO", "NONNEGATIVE", "NONPOSITIVE", or "UNKNOWN".
#' @rdname CallbackParam-class
#' @examples
#' x <- Variable(2)
#' y <- CallbackParam(value(x), dim(x), sign = "NONNEGATIVE")
#' get_data(y)
#' @export
CallbackParam <- function(callback, dim = NULL, ...) {
  .CallbackParam(callback = callback, dim = dim, ...)
}

setMethod("initialize", "CallbackParam", function(.Object, ..., callback, dim = NULL) {
  .Object@callback <- callback
  callNextMethod(.Object, ..., dim = dim)
})

#' @param object A \linkS4class{CallbackParam} object.
#' @rdname CallbackParam-class
setMethod("value", "CallbackParam", function(object) { validate_val(object, object@callback()) })   # TODO: Cast to vector if object@.is_vector == TRUE.
