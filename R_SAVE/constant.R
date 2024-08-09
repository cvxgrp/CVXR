#'
#' The Constant class.
#'
#' This class represents a constant value.
#'
#' Raw numerical constants (R primitive types, Matrix dense and sparse matrices,
#' et) are implicitly cast to constants via Expression operator overloading.
#' For example, if \code{x} is an expression and \code{c} is a raw constant,
#' then \code{x + c} creates an expression by casting \code{c} to a Constant.
#'
#' @slot value A numeric element, vector, matrix, or data.frame. Vectors are automatically cast into a matrix column.
#' @slot sparse (Internal) A logical value indicating whether the value is a sparse matrix.
#' @name Constant-class
#' @aliases Constant
#' @rdname Constant-class
.Constant <- setClass("Constant", representation(sparse = "logical", symm = "logical", herm = "logical", psd_test = "logical", nsd_test = "logical", skew_symm = "logical", .cached_is_pos = "logical"),
                                  prototype(sparse = NA, symm = NA, herm = NA, psd_test = NA, nsd_test = NA, skew_symm = NA, .cached_is_pos = NA), contains = "Leaf")

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

setMethod("initialize", "Constant", function(.Object, ..., value = NA_real_, sparse = NA, symm = NA, herm = NA, psd_test = NA, nsd_test = NA, skew_symm = NA, .cached_is_pos = NA) {
  # Keep sparse matrices sparse.
  if(is(value, "ConstSparseVal")) {
    .Object@value <- Matrix(value, sparse = TRUE)
    .Object@sparse <- TRUE
  } else {
    if(any(is.na(value)))
      .Object@value <- value
    else
      .Object@value <- as.matrix(value)
    .Object@sparse <- FALSE
  }
  .Object@symm <- NA
  .Object@herm <- NA
  .Object@psd_test <- NA
  .Object@nsd_test <- NA
  .Object@skew_symm <- NA
  .Object@.cached_is_pos <- NA
  callNextMethod(.Object, ..., dim = intf_dim(.Object@value), value = .Object@value)
})

#' @param x,object A \linkS4class{Constant} object.
#' @rdname Constant-class
setMethod("show", "Constant", function(object) {
  cat("Constant(", curvature(object), ", ", sign(object), ", (", paste(dim(object), collapse = ","), "))", sep = "")
})

#' @describeIn Constant The name of the constant.
setMethod("name", "Constant", function(x) { as.character(head(value(x))) })

#' @describeIn Constant Returns itself as a constant.
setMethod("constants", "Constant", function(object) { list(object) })

#' @describeIn Constant Is the expression a constant?
setMethod("is_constant", "Constant", function(object) { TRUE })

#' @describeIn Constant The value of the constant.
setMethod("value", "Constant", function(object) {
  # if(object@.is_vector)
  #  return(as.vector(object@value))
  return(object@value)
})

#' @describeIn Constant A logical value indicating whether all elements of the constant are positive.
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
    object <- .compute_attr(object)   # TODO: Modify slot in place.
  list(object, object@nonneg)
})

#' @describeIn Constant A logical value indicating whether all elements of the constant are non-positive.
setMethod("is_nonpos", "Constant", function(object) {
  if(is.na(object@nonpos))
    object <- .compute_attr(object)   # TODO: Modify slot in place.
  list(object, object@nonpos)
})

#' @describeIn Constant A logical value indicating whether the constant is imaginary.
setMethod("is_imag", "Constant", function(object) {
  if(is.na(object@imag))
    object <- .compute_attr(object)   # TODO: Modify slot in place.
  list(object, object@imag)
})

#' @describeIn Constant A logical value indicating whether the constant is complex-valued.
setMethod("is_complex", "Constant", function(object) { is.complex(value(object)) })

#' @describeIn Constant A logical value indicating whether the constant is symmetric.
setMethod("is_symmetric", "Constant", function(object) {
  if(is_scalar(object))
    symm <- TRUE
  else if(ndim(object) == 2 && nrow(object) == ncol(object)) {
    if(is.na(object@symm))
      object <- .compute_symm_attr(object)   # TODO: Modify slot in place.
    symm <- object@symm
  } else
    symm <- FALSE
  list(object, symm)
})

#' @describeIn Constant A logical value indicating whether the constant is a Hermitian matrix.
setMethod("is_hermitian", "Constant", function(object) {
  if(is_scalar(object) && is_real(object))
    herm <- TRUE
  else if(ndim(object) == 2 && nrow(object) == ncol(object)) {
    if(is.na(object@herm))
      object <- .compute_symm_attr(object)   # TODO: Modify slot in place.
    herm <- object@herm
  } else
    herm <- FALSE
  list(object, herm)
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

#' @describeIn Constant Is the expression skew symmetric?
setMethod("is_skew_symmetric", "Constant", function(object) {
  if(is.na(object@skew_symm))
    object@skew_symm <- intf_is_skew_symmetric(value(object))   # TODO: Modify slot in place.
  list(object, object@skew_symm)
})

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

  # Compute sign of bottom eigenvalue if absent.
  if(is.na(object@psd_test))
    object@psd_test <- is_psd_within_tol(value(object), EIGVAL_TOL)   # TODO: Modify slot in place.
  list(object, object@psd_test)
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

  # Compute sign of top eigenvalue if absent.
  if(is.na(object@nsd_test))
    object@nsd_test <- is_psd_within_tol(-value(object), EIGVAL_TOL)
  list(object, object@nsd_test)
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
  if(is(expr, "list")) {
    for(elem in expr) {
      if(is(elem, "Expression"))
        stop("The input must be a single CVXR Expression, not a list. Combine Expressions using atoms such as bmat, hstack, and vstack.")
    }
  }
  if(is(expr, "Expression"))
    expr
  else
    Constant(value = expr)
}

#'
#' The Parameter class.
#'
#' This class represents a parameter in an optimization problem.
#'
#' Parameters are constant expressions whose value may be specified
#' after problem creation. The only way to modify a problem after its
#' creation is through parameters. For example, you might choose to declare
#' the hyper-parameters of a machine learning model to be Parameter objects;
#' more generally, Parameters are useful for computing trade-off curves.
#'
#' @slot dim the dimensions of the parameter.
#' @slot name the name of this parameter
#' @slot PARAM_COUNT count of parameters, (currently unused?)
#' @slot venv an environment used for reference semantics in setting values for parameter
#' @slot delta the delta parameter
#' @slot gradient the gradient parameter
#' @slot is_constant flag indicating the parameter is constant
#' @slot is_vector flag indicating the parameter is a vector
#' @rdname Parameter-class
#' @examples
#' x <- Parameter(3, name = "x0", nonpos = TRUE) ## 3-vec negative
#' is_nonneg(x)
#' is_nonpos(x)
#' size(x)
#' @export
Parameter <- setClass("Parameter",
                      slots =
                        list(dim = "integer",
                             name = "character",
                             PARAM_COUNT = "integer",
                             venv = "environment", # so that we may have reference sematics for value!
                             delta = "numeric",
                             gradient = "numeric",
                             is_constant = "logical",
                             is_vector = "logical"
                             ),
                      contains = "Leaf")

setMethod("initialize", "Parameter",
          function(.Object, ..., dim = integer(0), name = NA_character_, value = NA_real_,
                   PARAM_COUNT = 0L, delta = NA_real_, gradient = NA_real_) {

            # This is handled in Canonical: .Object@id <- ifelse(is.na(id), get_id(), id)
            if (is.na(name)) {
              name <- sprintf("%s%s", PARAM_PREFIX, .Object@id)
            }
            .Object@name <- name
            d <- length(dim)
            if (d == 0L) {  # Force constants to default to c(1,1).
              dim <- c(1L, 1L)
              .Object@is_vector <- TRUE
            } else if (d == 1L) {  # Treat as a column vector.
              dim <- c(dim, 1L)
              .Object@is_vector <- TRUE
            } else if (d == 2L) {
              .Object@is_vector <- FALSE
            } else if (d > 2) {  # TODO: Tensors are currently unimplemented.
              stop("Unimplemented")
            } else {
              stop("Unimplemented")
            }
            .Object@dim <- dim
            .Object@delta <- delta
            .Object@gradient <- gradient
            .Object <- callNextMethod(.Object, ..., dim = dim, value = value)
            .Object@is_constant <- TRUE
            # Cache in virtual environment.
            .Object@venv <- new.env(parent = emptyenv())
            # Initialize with value if provided
            .Object@venv$value <- .Object@value
            return(.Object)
          })

#' @param object,x A \linkS4class{Parameter} object.
#' @describeIn Parameter Returns information needed to reconstruct the expression besides args: \code{list(dim, name, value, id, attributes)}.
setMethod("get_data", "Parameter", function(object) {
  list(dim(object), object@name, value(object), id(object), attributes(object))
})

#' @describeIn Parameter The name of the parameter.
#' @export
setMethod("name", "Parameter", function(x) { x@name })

#' @describeIn Parameter Is the expression constant?
setMethod("is_constant", "Parameter", function(object) {
  if(dpp_scope_active())   # TODO: Implement DPP scoping.
    return(FALSE)
  return(TRUE)
})

## # We also need a value_impl
setMethod("value_impl", "Parameter", function(object) {
  object@venv$value
})

#' @describeIn Parameter The value of the parameter.
setMethod("value", "Parameter", function(object) {
  # if(object@.is_vector)
  #  return(as.vector(object@value))
    ## return(object@value)
 object@venv$value
})

#' @describeIn Parameter Set the value of the parameter.
setReplaceMethod("value", "Parameter", function(object, value) {
  ## object@value <- validate_val(object, value)
  object@venv$value <- validate_val(object, value)
  object
})

#' @describeIn Parameter An empty list since the gradient of a parameter is zero.
setMethod("grad", "Parameter", function(object) { list() })

#' @describeIn Parameter Returns itself as a parameter.
setMethod("parameters", "Parameter", function(object) { list(object) })

#' @describeIn Parameter The canonical form of the parameter.
setMethod("canonicalize", "Parameter", function(object) {
  obj <- create_param(dim(object), id(object))
  list(obj, list())
})

setMethod("show", "Parameter", function(object) {
  attr_str <- get_attr_str(object)
  if(length(attr_str) > 0)
    print(paste("Parameter(", paste(dim(object), collapse = ", "), ", ", attr_str, ")", sep = ""))
  else
    print(paste("Parameter(", paste(dim(object), collapse = ", "), ")", sep = ""))
})

setMethod("as.character", "Parameter", function(x) {
  attr_str <- get_attr_str(x)
  if(length(attr_str) > 0)
    cat("Parameter(", paste(dim(x), collapse = ", "), ", ", attr_str, ")", sep = "")
  else
    cat("Parameter(", paste(dim(x), collapse = ", "), ")", sep = "")
})

is_param_affine <- function(expr) {
  # Returns TRUE if expression is parameters-affine (and variable-free).
  dpp_scope()   # TODO: Implement DPP scoping.
  length(variables(expr)) == 0 && is_affine(expr)
}

is_param_free <- function(expr) {
  # Returns TRUE if expression is not parametrized.
  length(parameters(expr)) == 0
}

#'
#' The CallbackParam class.
#'
#' This class represents a parameter whose value is obtained by evaluating a function.
#'
#' @slot callback A callback function that generates the parameter value.
#' @slot dim The dimensions of the parameter.
#' @name CallbackParam-class
#' @aliases CallbackParam
#' @rdname CallbackParam-class
#' @examples
#' x <- Variable(2)
#' fun <- function() { value(x) }
#' y <- CallbackParam(fun, dim = dim(x), nonneg = TRUE)
#' get_data(y)
#' @export
CallbackParam <- setClass("CallbackParam",
                          slots = list(callback = "function"),
                          contains = "Parameter")

#' @param object A \linkS4class{CallbackParam} object.
#' @rdname CallbackParam-class
setMethod("value", "CallbackParam", function(object) {
  object@venv@value <- validate_val(object, object@callback())
  # if(object@.is_vector)
  #   val <- as.vector(val)
  # return(val)
})

