## CVXPY SOURCE: cvxpy/reductions/reduction.py

## TODO: Is this classunion a repeat?
setClassUnion("ProblemORNULL", c("Problem", "NULL"))

#'
#' The Reduction class.
#'
#' This virtual class represents a reduction, an actor that transforms a problem
#' into an equivalent problem. By equivalent, we mean that there exists a mapping
#' between solutions of either problem: if we reduce a problem \eqn{A} to another
#' problem \eqn{B} and then proceed to find a solution to \eqn{B}, we can convert
#' it to a solution of \eqn{A} with at most a moderate amount of effort.
#'
#' A reduction that is instantiated with a non-NULL problem offers two key methods:
#' reduce and retrieve. The reduce method converts the problem the reduction
#' was instantiated with to an equivalent problem. The retrieve method takes
#' as an argument a Solution for the equivalent problem and returns a Solution
#' for the problem owned by the reduction.
#'
#' Every reduction supports three low-level methods: accepts, perform, and invert.
#' The accepts method of a particular reduction specifies the types of problems
#' that it is applicable to, the perform method takes a problem and reduces it to
#' an equivalent form, and the invert method maps solutions from reduced-to problems
#' to their problems of provenance.
#'
#' @rdname Reduction-class
setClass("Reduction", representation(problem = "ProblemORNULL", .emitted_problem = "ANY", .retrieval_data = "ANY"),
                      prototype(problem = NULL, .emitted_problem = NULL, .retrieval_data = NULL), contains = "VIRTUAL")

#' @param object A \linkS4class{Reduction} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Reduction States whether the reduction accepts a problem.
setMethod("accepts", signature(object = "Reduction", problem = "Problem"), function(object, problem) { stop("Unimplemented") })

#' @describeIn Reduction Reduces the owned problem to an equivalent problem.
setMethod("reduce", "Reduction", function(object) {
  if(!is.null(object@.emitted_problem))
    return(list(object, object@.emitted_problem))

  if(is.null(object@problem))
    stop("The reduction was constructed without a Problem")

  tmp <- perform(object, object@problem)
  object <- tmp[[1]]
  object@.emitted_problem <- tmp[[2]]
  object@.retrieval_data <- tmp[[3]]
  return(list(object, object@.emitted_problem))
})

#' @param object A \linkS4class{Reduction} object.
#' @param solution A \linkS4class{Solution} object.
#' @describeIn Reduction Retrieves a solution to the owned problem.
setMethod("retrieve", signature(object = "Reduction", solution = "Solution"), function(object, solution) {
  if(is.null(object@.retrieval_data))
    stop("reduce must be called before retrieve")
  return(invert(object, solution, object@.retrieval_data))
})

#' @param object A \linkS4class{Reduction} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Reduction Performs the reduction on a problem and returns an equivalent problem.
setMethod("perform", signature(object = "Reduction", problem = "Problem"), function(object, problem) { stop("Unimplemented") })

#' @param object A \linkS4class{Reduction} object.
#' @param solution A \linkS4class{Solution} to a problem that generated the inverse data.
#' @param inverse_data The data encoding the original problem.
#' @describeIn Reduction Returns a solution to the original problem given the inverse data.
setMethod("invert", signature(object = "Reduction", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) { stop("Unimplemented") })

