#'
#' The Constant class.
#'
#' This class represents a constant.
#'
#' @slot value A numeric element, data.frame, matrix, or vector.
#' @aliases Constant
#' @export
.Constant <- setClass("Constant", representation(value = "ConstVal", .sparse = "logical"), 
                                 prototype(value = NA_real_, .sparse = NA), 
                      validity = function(object) {
                        if((!is(object@value, "ConstSparseVal") && !is.data.frame(object@value) && !is.numeric(object@value)) ||
                           ((is(object@value, "ConstSparseVal") || is.data.frame(object@value)) && !all(sapply(object@value, is.numeric))))
                          stop("[Constant: validation] value must be a data.frame, matrix (CsparseMatrix, TsparseMatrix, or R default), vector, or atomic element containing only numeric entries")
                        return(TRUE)
                      }, contains = "Leaf")
Constant <- function(value) { .Constant(value = value) }

setMethod("init_dcp_attr", "Constant", function(object) {
  val_dim = dim(object@value)
  if(is.null(val_dim))
    shape = Shape(rows = 1, cols = length(object@value))
  else
    shape = Shape(rows = val_dim[1], cols = val_dim[2])
  sign = get_sign(object@value)
  DCPAttr(sign = sign, curvature = Curvature(curvature = CURV_CONSTANT_KEY), shape = shape)
})

setMethod("initialize", "Constant", function(.Object, ..., dcp_attr, value = NA_real_, .sparse = NA) {
  .Object@value <- value
  .Object@.sparse <- is(value, "ConstSparseVal")
  .Object@dcp_attr <- init_dcp_attr(.Object)
  callNextMethod(.Object, ..., dcp_attr = .Object@dcp_attr)
})

setMethod("get_data", "Constant", function(object) { list(object@value) })
setMethod("value", "Constant", function(object) { object@value })
setMethod("canonicalize", "Constant", function(object) {
  obj <- create_const(object@value, size(object), object@.sparse)
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
.Parameter <- setClass("Parameter", representation(rows = "numeric", cols = "numeric", name = "character", sign = "character", value = "ConstVal"),
                                    prototype(rows = 1, cols = 1, name = NA_character_, sign = SIGN_UNKNOWN_KEY, value = NA_real_), 
                      validity = function(object) {
                        if(!(object@sign %in% SIGN_STRINGS))
                          stop("[Sign: validation] sign must be in ", paste(SIGN_STRINGS, collapse = ", "))
                        else
                          return(TRUE)
                        }, contains = "Leaf")
Parameter <- function(rows = 1, cols = 1, name = NA_character_, sign = SIGN_UNKNOWN_KEY, value = NA_real_) {
  .Parameter(rows = rows, cols = cols, name = name, sign = toupper(sign), value = value)
}

setMethod("init_dcp_attr", "Parameter", function(object) {
  shape <- Shape(rows = object@rows, cols = object@cols)
  DCPAttr(sign = Sign(sign = object@sign), curvature = Curvature.CONSTANT, shape = shape)
})

setMethod("initialize", "Parameter", function(.Object, ..., dcp_attr, rows = 1, cols = 1, name = NA_character_, sign = SIGN_UNKNOWN_KEY, value = NA_real_) {
  .Object@rows <- rows
  .Object@cols <- cols
  .Object@sign <- sign
  .Object@name <- name
  .Object@value <- value
  .Object@dcp_attr <- init_dcp_attr(.Object)
  callNextMethod(.Object, ..., dcp_attr = .Object@dcp_attr)
})

setMethod("get_data", "Parameter", function(object) {
  list(rows = object@rows, cols = object@cols, name = object@name, sign = object@sign, value = object@value)
})
setMethod("parameters", "Parameter", function(object) { list(object) })
setMethod("canonicalize", "Parameter", function(object) {
  obj <- create_param(object, c(object@rows, object@cols))
  list(obj, list())
})

.CallbackParam <- setClass("CallbackParam", representation(.callback = "ConstVal"), contains = "Parameter")
CallbackParam <- function(.callback, rows = 1, cols = 1, name = NA_character_, sign = SIGN_UNKNOWN_KEY) {
  .CallbackParam(.callback = .callback, rows = rows, cols = cols, name = name, sign = sign)
}

setMethod("initialize", "CallbackParam", function(.Object, ..., .callback) {
  .Object@.callback <- .callback
  callNextMethod(.Object, ...)
})
