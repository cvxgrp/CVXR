% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elementwise.R
\docType{class}
\name{KLDiv-class}
\alias{KLDiv-class}
\alias{.KLDiv}
\alias{KLDiv}
\alias{to_numeric,KLDiv-method}
\alias{sign_from_args,KLDiv-method}
\alias{is_atom_convex,KLDiv-method}
\alias{is_atom_concave,KLDiv-method}
\alias{is_incr,KLDiv-method}
\alias{is_decr,KLDiv-method}
\alias{.grad,KLDiv-method}
\alias{.domain,KLDiv-method}
\title{The KLDiv class.}
\usage{
KLDiv(x, y)

\S4method{to_numeric}{KLDiv}(object, values)

\S4method{sign_from_args}{KLDiv}(object)

\S4method{is_atom_convex}{KLDiv}(object)

\S4method{is_atom_concave}{KLDiv}(object)

\S4method{is_incr}{KLDiv}(object, idx)

\S4method{is_decr}{KLDiv}(object, idx)

\S4method{.grad}{KLDiv}(object, values)

\S4method{.domain}{KLDiv}(object)
}
\arguments{
\item{x}{An \linkS4class{Expression} or numeric constant.}

\item{y}{An \linkS4class{Expression} or numeric constant.}

\item{object}{A \linkS4class{KLDiv} object.}

\item{values}{A list of numeric values for the arguments}

\item{idx}{An index into the atom.}
}
\description{
The elementwise KL-divergence \eqn{x\log(x/y) - x + y}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{to_numeric(KLDiv)}: The KL-divergence evaluted elementwise on the input value.

\item \code{sign_from_args(KLDiv)}: The atom is positive.

\item \code{is_atom_convex(KLDiv)}: The atom is convex.

\item \code{is_atom_concave(KLDiv)}: The atom is not concave.

\item \code{is_incr(KLDiv)}: The atom is not monotonic in any argument.

\item \code{is_decr(KLDiv)}: The atom is not monotonic in any argument.

\item \code{.grad(KLDiv)}: Gives the (sub/super)gradient of the atom w.r.t. each variable

\item \code{.domain(KLDiv)}: Returns constraints describng the domain of the node

}}
\section{Slots}{

\describe{
\item{\code{x}}{An \linkS4class{Expression} or numeric constant.}

\item{\code{y}}{An \linkS4class{Expression} or numeric constant.}
}}

