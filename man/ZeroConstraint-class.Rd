% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constraints.R
\docType{class}
\name{ZeroConstraint-class}
\alias{ZeroConstraint-class}
\alias{.ZeroConstraint}
\alias{name,ZeroConstraint-method}
\alias{dim,ZeroConstraint-method}
\alias{is_dcp,ZeroConstraint-method}
\alias{is_dgp,ZeroConstraint-method}
\alias{residual,ZeroConstraint-method}
\alias{canonicalize,ZeroConstraint-method}
\title{The ZeroConstraint class}
\usage{
\S4method{name}{ZeroConstraint}(x)

\S4method{dim}{ZeroConstraint}(x)

\S4method{is_dcp}{ZeroConstraint}(object)

\S4method{is_dgp}{ZeroConstraint}(object)

\S4method{residual}{ZeroConstraint}(object)

\S4method{canonicalize}{ZeroConstraint}(object)
}
\arguments{
\item{x, object}{A \linkS4class{ZeroConstraint} object.}
}
\description{
The ZeroConstraint class
}
\section{Methods (by generic)}{
\itemize{
\item \code{name(ZeroConstraint)}: The string representation of the constraint.

\item \code{dim(ZeroConstraint)}: The dimensions of the constrained expression.

\item \code{is_dcp(ZeroConstraint)}: Is the constraint DCP?

\item \code{is_dgp(ZeroConstraint)}: Is the constraint DGP?

\item \code{residual(ZeroConstraint)}: The residual of a constraint

\item \code{canonicalize(ZeroConstraint)}: The graph implementation of the object.

}}
