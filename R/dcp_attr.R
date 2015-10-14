#'
#' The Sign class.
#'
#' This class represents the sign of an expression.
#'
#' @slot sign A \code{character} string specifying the type of sign.
#' @aliases Sign
#' @export
.Sign <- setClass("Sign", 
                 representation(sign = "character"), 
                 prototype(sign = SIGN_UNKNOWN_KEY),
                 validity = function(object) {
                   if(!(object@sign %in% SIGN_STRINGS))
                     stop("[Sign: validation] sign must be in ", paste(SIGN_STRINGS, collapse = ", "))
                   else
                     return(TRUE)
                 })

Sign <- function(sign = SIGN_UNKNOWN_KEY) { .Sign(sign = toupper(sign)) }
Sign.POSITIVE <- Sign(sign = SIGN_POSITIVE_KEY)
Sign.NEGATIVE <- Sign(sign = SIGN_NEGATIVE_KEY)
Sign.ZERO <- Sign(sign = SIGN_ZERO_KEY)
Sign.UNKNOWN <- Sign(sign = SIGN_UNKNOWN_KEY)

setMethod("show", "Sign", function(object) { cat("Sign(", object@sign, ")", sep = "") })
setMethod("as.character", "Sign", function(x) { paste("Sign(", x@sign, ")", sep = "") })
setMethod("is_zero", "Sign", function(object) { object == Sign.ZERO })
setMethod("is_positive", "Sign", function(object) { is_zero(object) || object == Sign.POSITIVE })
setMethod("is_negative", "Sign", function(object) { is_zero(object) || object == Sign.NEGATIVE })
setMethod("is_unknown", "Sign", function(object) { object == Sign.UNKNOWN })

setMethod("+", c("Sign", "missing"), function(e1, e2) { e1 * Sign.POSITIVE })
setMethod("-", c("Sign", "missing"), function(e1, e2) { e1 * Sign.NEGATIVE })
setMethod("+", signature(e1 = "Sign", e2 = "Sign"), function(e1, e2) {
  if(is_zero(e1))
    e2
  else if(is_zero(e2))
    e1
  else if(is_positive(e1) && is_positive(e2))
    e1
  else if(is_negative(e1) && is_negative(e2))
    e1
  else
    Sign.UNKNOWN
})
setMethod("-", signature(e1 = "Sign", e2 = "Sign"), function(e1, e2) { e1 + -e2 })
setMethod("*", signature(e1 = "Sign", e2 = "Sign"), function(e1, e2) {
  if(is_zero(e1) || is_zero(e2))
    Sign.ZERO
  else if(is_unknown(e1) || is_unknown(e2))
    Sign.UNKNOWN
  else if(e1 != e2)
    Sign.NEGATIVE
  else
    Sign.POSITIVE
})
setMethod("==", signature(e1 = "Sign", e2 = "Sign"), function(e1, e2) { e1@sign == e2@sign })
setMethod("!=", signature(e1 = "Sign", e2 = "Sign"), function(e1, e2) { e1@sign != e2@sign })

val_to_sign <- function(val) {
  if(val > 0)
    Sign.POSITIVE
  else if(val == 0)
    Sign.ZERO
  else
    Sign.NEGATIVE
}

#'
#' The Curvature class.
#'
#' This class represents the curvature of an expression.
#'
#' @slot curvature A \code{character} string specifying the type of curvature.
#' @aliases Curvature
#' @export
Curvature <- setClass("Curvature", 
                       representation(curvature = "character"), 
                       prototype(curvature = CURV_UNKNOWN_KEY),
                       validity = function(object) {
                         if(!(object@curvature %in% CURVATURE_STRINGS))
                           stop("[Curvature: validation] curvature must be in ", paste(CURVATURE_STRINGS, collapse = ", "))
                         else
                           return(TRUE)
                        })
setMethod("show", "Curvature", function(object) { cat("Curvature(", object@curvature, ")", sep = "") })
setMethod("as.character", "Curvature", function(x) { paste("Curvature(", x@curvature, ")", sep = "") })
setMethod("is_constant", "Curvature", function(object) { object@curvature == CURV_CONSTANT_KEY })
setMethod("is_affine", "Curvature", function(object) { object@curvature == CURV_AFFINE_KEY })
setMethod("is_convex", "Curvature", function(object) { object@curvature == CURV_CONVEX_KEY })
setMethod("is_concave", "Curvature", function(object) { object@curvature == CURV_CONCAVE_KEY })
setMethod("is_unknown", "Curvature", function(object) { object@curvature == CURV_UNKNOWN_KEY })
setMethod("is_dcp", "Curvature", function(object) { object@curvature != CURV_UNKNOWN_KEY })

setMethod("+", c("Curvature", "missing"), function(e1, e2) { e1@curvature })
setMethod("-", c("Curvature", "missing"), function(e1, e2) {
  if(is_convex(e1))
    Curvature(curvature = CURV_CONCAVE_KEY)
  else if(is_concave(e1))
    Curvature(curvature = CURV_CONVEX_KEY)
  else
    e1
})
setMethod("+", signature(e1 = "Curvature", e2 = "Curvature"),
               function(e1, e2) {
                 if(is_constant(e1))
                   e2
                 else if(is_affine(e1) && is_affine(e2))
                   Curvature(curvature = CURV_AFFINE_KEY)
                 else if(is_convex(e1) && is_convex(e2))
                   Curvature(curvature = CURV_CONVEX_KEY)
                 else if(is_concave(e1) && is_concave(e2))
                   Curvature(curvature = CURV_CONCAVE_KEY)
                 else
                   Curvature(curvature = CURV_UNKNOWN_KEY)
               })
setMethod("-", signature(e1 = "Curvature", e2 = "Curvature"), function(e1, e2) { e1 + -e2 })
setMethod("==", signature(e1 = "Curvature", e2 = "Curvature"), function(e1, e2) { e1@curvature == e2@curvature })
setMethod("!=", signature(e1 = "Curvature", e2 = "Curvature"), function(e1, e2) { e1@curvature != e2@curvature })

setMethod("sign_mul", signature(sign = "Sign", curv = "Curvature"), function(sign, curv) {
  if(is_zero(sign))
    Curvature(curvature = CURV_CONSTANT_KEY)
  else if(is_positive(sign) || is_affine(curv))
    curv
  else if(is_negative(sign)) {
    new_curv = CURVATURE_NEGATION_MAP(curv@curvature)
    Curvature(curvature = new_curv)
  } else
    Curvature(curvature = CURV_UNKNOWN_KEY)
})

#'
#' The Shape class.
#'
#' This class represents the sign of an expression.
#'
#' @slot rows The number of rows in the matrix
#' @slot cols The number of columns in the matrix
#' @aliases Shape
#' @export
Shape <- setClass("Shape", representation(rows = "numeric", cols = "numeric"), prototype(rows = NA_integer_, cols = NA_integer_))
setMethod("show", "Shape", function(object) { cat("Shape(", object@rows, ", ", object@cols, ")", sep = "") })
setMethod("as.character", "Shape", function(x) { paste("Shape(", x@rows, ", ", x@cols, ")", sep = "") })

setMethod("size", "Shape", function(object) { c(object@rows, object@cols) })
setMethod("+", signature(e1 = "Shape", e2 = "Shape"), function(e1, e2) {
  if(all(size(e1) == c(1,1)))
    e2
  else if(all(size(e2) == c(1,1)))
    e1
  else if(all(size(e1) == size(e2)))
    e1
  else
    stop("Incompatible dimensions: ", as.character(e1), " vs. ", as.character(e2))
})
setMethod("-", signature(e1 = "Shape", e2 = "Shape"), function(e1, e2) { e1 + e2 })
setMethod("*", signature(e1 = "Shape", e2 = "Shape"), function(e1, e2) {
  if(all(size(e1) == c(1,1)))
    e2
  else if(all(size(e2) == c(1,1)))
    e1
  else if(e1@cols == e2@rows)
    Shape(rows = e1@rows, cols = e2@cols)
  else
    stop("Incompatible dimensions: ", as.character(e1), " vs. ", as.character(e2))
})
setMethod("/", signature(e1 = "Shape", e2 = "Shape"), function(e1, e2) { e1 })

#'
#' The DCPAttr class.
#'
#' This class represents the attributes of a DCP (disciplined convex program).
#'
#' @slot sign The signs of the entries in the matrix expression.
#' @slot curvature The curvatures of the entries in the matrix expression.
#' @slot shape The dimensions of the matrix expression.
#' @aliases DCPAttr
#' @export
DCPAttr <- setClass("DCPAttr", representation(sign = "Sign", curvature = "Curvature", shape = "Shape"),
                    prototype(sign = new("Sign"), curvature = new("Curvature"), shape = new("Shape")))
setMethod("show", "DCPAttr", function(object) {
  cat("DCPAttr(", as.character(object@sign), ", ", as.character(object@curvature), ", ", as.character(object@shape), ")", sep = "")
})

setMethod("+", signature(e1 = "DCPAttr", e2 = "missing"), function(e1, e2) { e1 })
setMethod("-", signature(e1 = "DCPAttr", e2 = "missing"), function(e1, e2) {
  DCPAttr(sign = -e1@sign, curvature = -e1@curvature, shape = e1@shape)
})
setMethod("+", signature(e1 = "DCPAttr", e2 = "DCPAttr"), function(e1, e2) {
  shape <- e1@shape + e2@shape
  sign <- e1@sign + e2@sign
  curvature <- e1@curvature + e2@curvature
  DCPAttr(sign = sign, curvature = curvature, shape = shape)
})
setMethod("-", signature(e1 = "DCPAttr", e2 = "DCPAttr"), function(e1, e2) {
  shape <- e1@shape - e2@shape
  sign <- e1@sign - e2@sign
  curvature <- e1@curvature - e2@curvature
  DCPAttr(sign = sign, curvature = curvature, shape = shape)
})
setMethod("*", signature(e1 = "DCPAttr", e2 = "DCPAttr"), function(e1, e2) {
  shape <- e1@shape * e2@shape
  sign <- e1@sign * e2@sign
  if(is_constant(e1@curvature))
    curvature <- sign_mul(e1@sign, e2@curvature)
  else
    curvature <- sign_mul(e2@sign, e1@curvature)
  DCPAttr(sign = sign, curvature = curvature, shape = shape)
})
setMethod("/", signature(e1 = "DCPAttr", e2 = "DCPAttr"), function(e1, e2) { e2 * e1 })

setMethod("DCPAttr.mul_elemwise", signature(lh_exp = "DCPAttr", rh_exp = "DCPAttr"), function(lh_exp, rh_exp) {
  shape <- lh_exp@shape + rh_exp@shape
  sign <- lh_exp@sign * rh_exp@sign
  curvature <- sign_mul(lh_exp@sign, rh_exp@curvature)
  DCPAttr(sign = sign, curvature = curvature, shape = shape)
})

#'
#' Applies DCP composition rules to determine curvature in each argument
#' 
#' Composition rules:
#' Key: Function curvature + monotonicity + argument curvature == curvature in argument
#'   anything + anything + constant == constant
#'   anything + anything + affine == original curvature
#'   convex/affine + increasing + convex == convex
#'   convex/affine + decreasing + concave == convex
#'   concave/affine + increasing + concave == concave
#'   concave/affine + decreasing + convex == concave
#' Notes: Increasing (decreasing) means non-decreasing (non-increasing).
#' Any combinations not covered by the rules result in a nonconvex expression.
#'
#' @param monotonicity: The monotonicity of the function in the given argument.
#' @param func_curvature: The curvature of the function.
#' @param arg_sign: The sign of the given argument.
#' @param arg_curvature: The curvature of the given argument.
#' @export
setMethod("dcp_curvature", signature(monotonicity = "character", func_curvature = "Curvature", arg_sign = "Sign", arg_curvature = "Curvature"),
            function(monotonicity, func_curvature, arg_sign, arg_curvature) {
              if(!(monotonicity %in% MONOTONICITY_STRINGS))
                stop("monotonicity must be in ", paste(MONOTONICITY_STRINGS, collapse = ", "))
              
              if(is_constant(arg_curvature))
                result_curv <- Curvature(curvature = CURV_CONSTANT_KEY)
              else if(is_affine(arg_curvature))
                result_curv <- func_curvature
              else if(monotonicity == INCREASING)
                result_curv <- func_curvature + arg_curvature
              else if(monotonicity == DECREASING)
                result_curv <- func_curvature - arg_curvature
              else if(monotonicity == SIGNED && is_convex(func_curvature)) {
                if((is_convex(arg_curvature) && is_positive(arg_sign)) || (is_concave(arg_curvature) && is_negative(arg_sign)))
                  result_curv <- func_curvature
                else
                  result_curv <- Curvature(curvature = CURV_UNKNOWN_KEY)
              } else   # Non-monotonic
                result_curv <- func_curvature + arg_curvature - arg_curvature  # TODO: Is this correct? Why add/subtract?
              result_curv
            })
