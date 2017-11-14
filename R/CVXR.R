#' CVXR: Disciplined Convex Optimization in R
#'
#' CVXR is an R package that provides an object-oriented modeling
#' language for convex optimization, similar to CVX, CVXPY, YALMIP,
#' and Convex.jl. This domain specific language (DSL) allows the user
#' to formulate convex optimization problems in a natural mathematical
#' syntax rather than the restrictive standard form required by most
#' solvers. The user specifies an objective and set of constraints by
#' combining constants, variables, and parameters using a library of
#' functions with known mathematical properties. CVXR then applies
#' signed disciplined convex programming (DCP) to verify the problem's
#' convexity. Once verified, the problem is converted into standard
#' conic form using graph implementations and passed to a cone solver
#' such as ECOS or SCS.
#'
#' @name CVXR-package
#' @useDynLib CVXR
#' @importFrom Rcpp evalCpp
#' @import methods
#' @import Matrix
#' @importFrom stats rnorm runif
#' @aliases CVXR-package CVXR
#' @docType package
#' @author Anqi Fu, Balasubramanian Narasimhan, John Miller, Steven Diamond, Stephen Boyd
#'
#' Maintainer: Anqi Fu<anqif@stanford.edu>
#' @keywords package
NULL



