#'
#' CVXR Package Constants
#'
#' Curvature types as strings
#'
CURV_CONSTANT_KEY = "CONSTANT"
CURV_AFFINE_KEY = "AFFINE"
CURV_CONVEX_KEY = "CONVEX"
CURV_CONCAVE_KEY = "CONCAVE"
CURV_UNKNOWN_KEY = "UNKNOWN"
CURVATURE_STRINGS = c(CURV_CONSTANT_KEY, CURV_AFFINE_KEY, CURV_CONVEX_KEY, CURV_CONCAVE_KEY, CURV_UNKNOWN_KEY)
CURVATURE_NEGATION_MAP <- function(curvature) {
  if(curvature == CURV_CONVEX_KEY)
    CONCAVE_KEY
  else if(curvature == CURV_CONCAVE_KEY)
    CONVEX_KEY
  else if(curvature %in% CURVATURE_STRINGS)
    curvature
  else
    stop("Curvature type ", curvature, " not recognized")
}
  
#'
#' Sign types as strings
#'
SIGN_POSITIVE_KEY = "POSITIVE"
SIGN_NEGATIVE_KEY = "NEGATIVE"
SIGN_UNKNOWN_KEY = "UNKNOWN"
SIGN_ZERO_KEY = "ZERO"
SIGN_STRINGS = c(SIGN_POSITIVE_KEY, SIGN_NEGATIVE_KEY, SIGN_UNKNOWN_KEY, SIGN_ZERO_KEY)

#'
#' Monotonicity types as strings
#' 
INCREASING = "INCREASING"
DECREASING = "DECREASING"
SIGNED = "SIGNED"
NONMONOTONIC = "NONMONOTONIC"
MONOTONICITY_STRINGS = c(INCREASING, DECREASING, SIGNED, NONMONOTONIC)

flatten_list <- function(x) {
  y <- list()
  rapply(x, function(x) y <<- c(y,x))
  y
}
