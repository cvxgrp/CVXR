% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constant.R
\docType{class}
\name{Constant-class}
\alias{Constant-class}
\alias{.Constant}
\alias{Constant}
\alias{Constant}
\alias{as.character,Constant-method}
\alias{constants,Constant-method}
\alias{get_data,Constant-method}
\alias{value,Constant-method}
\alias{grad,Constant-method}
\alias{size,Constant-method}
\alias{is_positive,Constant-method}
\alias{is_negative,Constant-method}
\alias{canonicalize,Constant-method}
\alias{as.Constant}
\title{The Constant class.}
\usage{
Constant(value)

\S4method{as.character}{Constant}(x)

\S4method{constants}{Constant}(object)

\S4method{get_data}{Constant}(object)

\S4method{value}{Constant}(object)

\S4method{grad}{Constant}(object)

\S4method{size}{Constant}(object)

\S4method{is_positive}{Constant}(object)

\S4method{is_negative}{Constant}(object)

\S4method{canonicalize}{Constant}(object)

as.Constant(expr)
}
\arguments{
\item{value}{A numeric element, vector, matrix, or data.frame. Vectors are automatically cast into a matrix column.}

\item{x, object}{A \linkS4class{Constant} object.}

\item{expr}{An \linkS4class{Expression}, numeric element, vector, matrix, or data.frame.}
}
\value{
A \linkS4class{Constant} representing the input as a constant.
}
\description{
This class represents a constant.

Coerce an R object or expression into the \linkS4class{Constant} class.
}
\section{Methods (by generic)}{
\itemize{
\item \code{constants}: Returns itself as a constant.

\item \code{get_data}: A list with the value of the constant.

\item \code{value}: The value of the constant.

\item \code{grad}: An empty list since the gradient of a constant is zero.

\item \code{size}: The \code{c(row, col)} dimensions of the constant.

\item \code{is_positive}: A logical value indicating whether all elemenets of the constant are non-negative.

\item \code{is_negative}: A logical value indicating whether all elemenets of the constant are non-positive.

\item \code{canonicalize}: The canonical form of the constant.
}}

\section{Slots}{

\describe{
\item{\code{value}}{A numeric element, vector, matrix, or data.frame. Vectors are automatically cast into a matrix column.}

\item{\code{is_1D_array}}{(Internal) A logical value indicating whether the value is a vector or 1-D matrix.}

\item{\code{sparse}}{(Internal) A logical value indicating whether the value is a sparse matrix.}

\item{\code{size}}{(Internal) A vector of containing the number of rows and columns.}

\item{\code{is_pos}}{(Internal) A logical value indicating whether all elements are non-negative.}

\item{\code{is_neg}}{(Internal) A logical value indicating whether all elements are non-positive.}
}}

\examples{
x <- Constant(5)
y <- Constant(diag(3))
get_data(y)
value(y)
is_positive(y)
size(y)
as.Constant(y)
}