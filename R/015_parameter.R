## CVXPY SOURCE: cvxpy/expressions/constants/parameter.py
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
  obj <- lu.create_param(dim(object), id(object))
  list(obj, list())
})

setMethod("as.character", "Parameter", function(x) {
  dims <- glue::glue_collapse("{dim(x)}", sep = ", ")
  attr_str <- get_attr_str(object)
  if(length(attr_str) > 0)
    glue::glue("Parameter[{dims}]({attr_str})")    
  else
    glue::glue("Parameter[{dims}]")        
})

setMethod("show", "Parameter", function(object) {
  cli::text(as.character(object))
})

is_param_affine <- function(expr) {
  # Returns TRUE if expression is parameters-affine (and variable-free).
  saved_scope <- .CVXR_options$dpp_scope_active
  .CVXR_options$dpp_scope_active <- TRUE
  on.exit({
    .CVXR_options$dpp_scope_active <- saved_scope
  })
  return (length(variables(expr)) == 0 && is_affine(expr) )
}

is_param_free <- function(expr) {
  # Returns TRUE if expression is not parametrized.
  length(parameters(expr)) == 0
}
