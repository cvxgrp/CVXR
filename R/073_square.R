## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/square.py

#'
#' Square
#'
#' The elementwise square.
#'
#' @param x An \linkS4class{Expression}.
#' @return An \linkS4class{Expression} representing the square of the input.
#' A <- Variable(2,2)
#' val <- cbind(c(2,4), c(16,1))
#' prob <- Problem(Minimize(square(A)[1,2]), list(A == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @aliases square
#' @rdname square
#' @export
setMethod("square", "Expression", function(x) { Power(x = x, p = 2) })
