## CVXPY SOURCE: cvxpy/problems/param_prob.py

#'
#' The ParamProb class.
#'
#' This virtual class represents parametrized problems.
#'
#' Parameterized problems are produced during the first canonicalization
#' and allow canonicalization to be short-circuited for future solves.
#'
#' @name ParamProb-class
#' @aliases ParamProb
#' @rdname ParamProb-class
ParamProb <- setClass("ParamProb", contains = "VIRTUAL")

#' @param object A \linkS4class{ParamProb} object.
#' @describeIn ParamProb Is the problem mixed-integer?
setMethod("is_mixed_integer", "ParamProb", function(object) { stop("Unimplemented") })

#' @param id_to_param_value (Optional) List mapping parameter IDs to values.
#' @param zero_offset (Optional) If TRUE, zero out the constant offset in the parameter vector.
#' @param keep_zeros (Optional) If TRUE, store explicit zeros in A where parameters are affected.
#' @describeIn ParamProb Returns A, b after applying parameters (and reshaping).
setMethod("apply_parameters", "ParamProb", function(object, id_to_param_value = NULL, zero_offset = FALSE, keep_zeros = FALSE) {
  stop("Unimplemented")
})

