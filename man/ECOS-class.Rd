% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conic_solvers.R
\docType{class}
\name{ECOS-class}
\alias{ECOS-class}
\alias{ECOS}
\alias{mip_capable,ECOS-method}
\alias{status_map,ECOS-method}
\alias{import_solver,ECOS-method}
\alias{name,ECOS-method}
\alias{perform,ECOS,Problem-method}
\alias{invert,ECOS,list,list-method}
\title{An interface for the ECOS solver}
\usage{
ECOS()

\S4method{mip_capable}{ECOS}(solver)

\S4method{status_map}{ECOS}(solver, status)

\S4method{import_solver}{ECOS}(solver)

\S4method{name}{ECOS}(x)

\S4method{perform}{ECOS,Problem}(object, problem)

\S4method{invert}{ECOS,list,list}(object, solution, inverse_data)
}
\arguments{
\item{solver, object, x}{A \linkS4class{ECOS} object.}

\item{status}{A status code returned by the solver.}

\item{problem}{A \linkS4class{Problem} object.}

\item{solution}{The raw solution returned by the solver.}

\item{inverse_data}{A list containing data necessary for the inversion.}
}
\description{
An interface for the ECOS solver
}
\section{Methods (by generic)}{
\itemize{
\item \code{mip_capable(ECOS)}: Can the solver handle mixed-integer programs?

\item \code{status_map(ECOS)}: Converts status returned by the ECOS solver to its respective CVXPY status.

\item \code{import_solver(ECOS)}: Imports the solver

\item \code{name(ECOS)}: Returns the name of the solver

\item \code{perform(object = ECOS, problem = Problem)}: Returns a new problem and data for inverting the new solution.

\item \code{invert(object = ECOS, solution = list, inverse_data = list)}: Returns the solution to the original problem given the inverse_data.

}}
