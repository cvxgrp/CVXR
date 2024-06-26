% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clarabel.R
\docType{class}
\name{CLARABEL-class}
\alias{CLARABEL-class}
\alias{CLARABEL}
\alias{mip_capable,CLARABEL-method}
\alias{status_map,CLARABEL-method}
\alias{name,CLARABEL-method}
\alias{import_solver,CLARABEL-method}
\alias{reduction_format_constr,CLARABEL-method}
\alias{perform,CLARABEL,Problem-method}
\alias{invert,CLARABEL,list,list-method}
\alias{solve_via_data,CLARABEL-method}
\title{An interface for the CLARABEL solver}
\usage{
CLARABEL()

\S4method{mip_capable}{CLARABEL}(solver)

\S4method{status_map}{CLARABEL}(solver, status)

\S4method{name}{CLARABEL}(x)

\S4method{import_solver}{CLARABEL}(solver)

\S4method{reduction_format_constr}{CLARABEL}(object, problem, constr, exp_cone_order)

\S4method{perform}{CLARABEL,Problem}(object, problem)

\S4method{invert}{CLARABEL,list,list}(object, solution, inverse_data)

\S4method{solve_via_data}{CLARABEL}(
  object,
  data,
  warm_start,
  verbose,
  feastol,
  reltol,
  abstol,
  num_iter,
  solver_opts,
  solver_cache
)
}
\arguments{
\item{solver, object, x}{A \linkS4class{CLARABEL} object.}

\item{status}{A status code returned by the solver.}

\item{problem}{A \linkS4class{Problem} object.}

\item{constr}{A \linkS4class{Constraint} to format.}

\item{exp_cone_order}{A list indicating how the exponential cone arguments are ordered.}

\item{solution}{The raw solution returned by the solver.}

\item{inverse_data}{A list containing data necessary for the inversion.}

\item{data}{Data generated via an apply call.}

\item{warm_start}{A boolean of whether to warm start the solver.}

\item{verbose}{A boolean of whether to enable solver verbosity.}

\item{feastol}{The feasible tolerance on the primal and dual residual.}

\item{reltol}{The relative tolerance on the duality gap.}

\item{abstol}{The absolute tolerance on the duality gap.}

\item{num_iter}{The maximum number of iterations.}

\item{solver_opts}{A list of Solver specific options}

\item{solver_cache}{Cache for the solver.}
}
\description{
An interface for the CLARABEL solver
}
\section{Methods (by generic)}{
\itemize{
\item \code{mip_capable(CLARABEL)}: Can the solver handle mixed-integer programs?

\item \code{status_map(CLARABEL)}: Converts status returned by CLARABEL solver to its respective CVXPY status.

\item \code{name(CLARABEL)}: Returns the name of the solver

\item \code{import_solver(CLARABEL)}: Imports the solver

\item \code{reduction_format_constr(CLARABEL)}: Return a linear operator to multiply by PSD constraint coefficients.

\item \code{perform(object = CLARABEL, problem = Problem)}: Returns a new problem and data for inverting the new solution

\item \code{invert(object = CLARABEL, solution = list, inverse_data = list)}: Returns the solution to the original problem given the inverse_data.

\item \code{solve_via_data(CLARABEL)}: Solve a problem represented by data returned from apply.

}}
