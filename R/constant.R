#'
#' The Constant class.
#'
#' This class represents a constant.
#'
#' @slot value A numeric element, data.frame, matrix, or vector.
#' @aliases Constant
#' @export
.Constant <- setClass("Constant", representation(value = "ConstVal", is_1D_array = "logical", sparse = "logical", size = "numeric", is_pos = "logical", is_neg = "logical"), 
                                 prototype(value = NA_real_, is_1D_array = FALSE, sparse = NA, size = NA_real_, is_pos = NA, is_neg = NA), 
                      validity = function(object) {
                        if((!is(object@value, "ConstSparseVal") && !is.data.frame(object@value) && !is.numeric(object@value)) ||
                           ((is(object@value, "ConstSparseVal") || is.data.frame(object@value)) && !all(sapply(object@value, is.numeric))))
                          stop("[Constant: validation] value must be a data.frame, matrix (CsparseMatrix, TsparseMatrix, or R default), vector, or atomic element containing only numeric entries")
                        return(TRUE)
                      }, contains = "Leaf")
Constant <- function(value) { .Constant(value = value) }

setMethod("initialize", "Constant", function(.Object, ..., value = NA_real_, is_1D_array = FALSE, .sparse = NA, .size = NA_real_, .is_pos = NA, .is_neg = NA) {
  .Object@is_1D_array <- is_1D_array
  .Object@value <- value
  if(is(value, "ConstSparseVal")) {
    .Object@value <- Matrix(value, sparse = TRUE)
    .Object@sparse <- TRUE
  } else {
    if(is.vector(value) && length(vector) == 1)
      .Object@is_1D_array <- TRUE
    .Object@value <- as.matrix(value)
    .Object@sparse <- FALSE
  }
  .Object@size <- intf_size(.Object@value)
  sign <- intf_sign(.Object@value)
  .Object@is_pos <- sign[1]
  .Object@is_neg <- sign[2]
  callNextMethod(.Object, ...)
})

setMethod("show", "Constant", function(object) {
  cat("Constant(", curvature(object), ", ", sign(object), ", (", paste(size(object), collapse = ","), "))", sep = "")
})

setMethod("as.character", "Constant", function(x) {
  paste("Constant(", curvature(x), ", ", sign(x), ", (", paste(size(x), collapse = ","), "))", sep = "")
})

setMethod("constants", "Constant", function(object) { list(object) })
setMethod("get_data", "Constant", function(object) { list(value(object)) })
setMethod("value", "Constant", function(object) { object@value })
setMethod("grad", "Constant", function(object) { list() })
setMethod("size", "Constant", function(object) { object@size })
setMethod("is_positive", "Constant", function(object) { object@is_pos })
setMethod("is_negative", "Constant", function(object) { object@is_neg })

setMethod("canonicalize", "Constant", function(object) {
  obj <- create_const(value(object), size(object), object@sparse)
  list(obj, list())
})

as.Constant <- function(expr) {
  if(is(expr, "Expression"))
    expr
  else
    Constant(value = expr)
}

get_sign <- function(constant) {
  if(!is(constant, "ConstVal"))
    stop("constant must be a data.frame, matrix, vector, or atomic value")
  if((!is(constant, "ConstSparseVal") && !is.data.frame(constant) && !is.numeric(constant)) ||
     ((is(constant, "ConstSparseVal") || is.data.frame(constant)) && !all(sapply(constant, is.numeric))))
    stop("constant must contain only numeric values")
  max_sign <- val_to_sign(max(constant))
  min_sign <- val_to_sign(min(constant))
  max_sign + min_sign
}

#'
#' The Parameter class.
#'
#' This class represents a parameter.
#'
#' @slot rows The number of rows in the parameter.
#' @slot cols The number of columns in the parameter.
#' @slot name (Optional) A character string indicating the name of the parameter
#' @slot sign A character string indicating the sign of the parameter. Must be "POSITIVE" or "NEGATIVE"
#' @slot value Holds the solved numeric value of the parameter.
#' @aliases Parameter
#' @export
.Parameter <- setClass("Parameter", representation(id = "integer", rows = "numeric", cols = "numeric", name = "character", sign_str = "character", value = "ConstVal"),
                                    prototype(rows = 1, cols = 1, name = NA_character_, sign_str = UNKNOWN, value = NA_real_), 
                      validity = function(object) {
                        if(!(object@sign_str %in% SIGN_STRINGS))
                          stop("[Sign: validation] sign_str must be in ", paste(SIGN_STRINGS, collapse = ", "))
                        else
                          return(TRUE)
                        }, contains = "Leaf")
Parameter <- function(rows = 1, cols = 1, name = NA_character_, sign = UNKNOWN, value = NA_real_) {
  .Parameter(rows = rows, cols = cols, name = name, sign_str = toupper(sign), value = value)
}

setMethod("initialize", "Parameter", function(.Object, ..., id = get_id(), rows = 1, cols = 1, name = NA_character_, sign_str = UNKNOWN, value = NA_real_) {
  .Object@id <- id
  .Object@rows <- rows
  .Object@cols <- cols
  .Object@sign_str <- sign_str
  if(is.na(name))
    .Object@name <- sprintf("%s%s", PARAM_PREFIX, .Object@id)
  else
    .Object@name <- name
  
  # Initialize with value if provided
  .Object@value <- NA_real_
  if(!(length(value) == 1 && is.na(value)))
    value(.Object) <- value
  callNextMethod(.Object, ...)
})

setMethod("show", "Parameter", function(object) {
  cat("Parameter(", object@rows, ", ", object@cols, ", sign = ", sign(object), ")", sep = "")
})

setMethod("as.character", "Parameter", function(x) {
  paste("Expression(", x@rows, ", ", x@cols, ", sign = ", sign(x), ")", sep = "")
})

setMethod("get_data", "Parameter", function(object) {
  list(rows = object@rows, cols = object@cols, name = object@name, sign_str = object@sign_str, value = object@value)
})

setMethod("name", "Parameter", function(object) { object@name })
setMethod("size", "Parameter", function(object) { c(object@rows, object@cols) })
setMethod("is_positive", "Parameter", function(object) { object@sign_str == ZERO || toupper(object@sign_str) == POSITIVE })
setMethod("is_negative", "Parameter", function(object) { object@sign_str == ZERO || toupper(object@sign_str) == NEGATIVE })
setMethod("value", "Parameter", function(object) { object@value })
setReplaceMethod("value", "Parameter", function(object, value) {
  object@value <- validate_val(object, value)
  object
})
setMethod("grad", "Parameter", function(object) { list() })
setMethod("parameters", "Parameter", function(object) { list(object) })
setMethod("canonicalize", "Parameter", function(object) {
  obj <- create_param(object, size(object))
  list(obj, list())
})

.CallbackParam <- setClass("CallbackParam", representation(callback = "ConstVal"), contains = "Parameter")
CallbackParam <- function(callback, rows = 1, cols = 1, name = NA_character_, sign = UNKNOWN) {
  .CallbackParam(callback = callback, rows = rows, cols = cols, name = name, sign_str = sign)
}

setMethod("get_data", "CallbackParam", function(object) {
  list(callback = object@callback, rows = object@rows, cols = object@cols, name = object@name, sign_str = object@sign_str)
})
