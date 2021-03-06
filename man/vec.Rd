% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exports.R
\docType{methods}
\name{vec}
\alias{vec}
\title{Vectorization of a Matrix}
\usage{
vec(X)
}
\arguments{
\item{X}{An \linkS4class{Expression} or matrix.}
}
\value{
An \linkS4class{Expression} representing the vectorized matrix.
}
\description{
Flattens a matrix into a vector in column-major order.
}
\examples{
A <- Variable(2,2)
c <- 1:4
expr <- vec(A)
obj <- Minimize(t(expr) \%*\% c)
constraints <- list(A == cbind(c(-1,-2), c(3,4)))
prob <- Problem(obj, constraints)
result <- solve(prob)
result$value
result$getValue(expr)
}
