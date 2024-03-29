% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/atoms.R
\docType{class}
\name{NormNuc-class}
\alias{NormNuc-class}
\alias{.NormNuc}
\alias{NormNuc}
\alias{to_numeric,NormNuc-method}
\alias{allow_complex,NormNuc-method}
\alias{dim_from_args,NormNuc-method}
\alias{sign_from_args,NormNuc-method}
\alias{is_atom_convex,NormNuc-method}
\alias{is_atom_concave,NormNuc-method}
\alias{is_incr,NormNuc-method}
\alias{is_decr,NormNuc-method}
\alias{.grad,NormNuc-method}
\title{The NormNuc class.}
\usage{
NormNuc(A)

\S4method{to_numeric}{NormNuc}(object, values)

\S4method{allow_complex}{NormNuc}(object)

\S4method{dim_from_args}{NormNuc}(object)

\S4method{sign_from_args}{NormNuc}(object)

\S4method{is_atom_convex}{NormNuc}(object)

\S4method{is_atom_concave}{NormNuc}(object)

\S4method{is_incr}{NormNuc}(object, idx)

\S4method{is_decr}{NormNuc}(object, idx)

\S4method{.grad}{NormNuc}(object, values)
}
\arguments{
\item{A}{An \linkS4class{Expression} or numeric matrix.}

\item{object}{A \linkS4class{NormNuc} object.}

\item{values}{A list of numeric values for the arguments}

\item{idx}{An index into the atom.}
}
\description{
The nuclear norm, i.e. sum of the singular values of a matrix.
}
\section{Methods (by generic)}{
\itemize{
\item \code{to_numeric(NormNuc)}: The nuclear norm (i.e., the sum of the singular values) of \code{A}.

\item \code{allow_complex(NormNuc)}: Does the atom handle complex numbers?

\item \code{dim_from_args(NormNuc)}: The atom is a scalar.

\item \code{sign_from_args(NormNuc)}: The atom is positive.

\item \code{is_atom_convex(NormNuc)}: The atom is convex.

\item \code{is_atom_concave(NormNuc)}: The atom is not concave.

\item \code{is_incr(NormNuc)}: The atom is not monotonic in any argument.

\item \code{is_decr(NormNuc)}: The atom is not monotonic in any argument.

\item \code{.grad(NormNuc)}: Gives the (sub/super)gradient of the atom w.r.t. each variable

}}
\section{Slots}{

\describe{
\item{\code{A}}{An \linkS4class{Expression} or numeric matrix.}
}}

