% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\docType{methods}
\name{perform}
\alias{perform}
\title{Perform Reduction}
\usage{
perform(object, problem)
}
\arguments{
\item{object}{A \linkS4class{Reduction} object.}

\item{problem}{A \linkS4class{Problem} on which the reduction will be performed.}
}
\value{
A list containing
\describe{
   \item{"problem"}{A \linkS4class{Problem} or list representing the equivalent problem.}
   \item{"inverse_data"}{A \linkS4class{InverseData} or list containing the data needed to invert this particular reduction.}
}
}
\description{
Performs the reduction on a problem and returns an equivalent problem.
}
