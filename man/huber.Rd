% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exports.R
\docType{methods}
\name{huber}
\alias{huber}
\title{Huber Function}
\usage{
huber(x, M = 1)
}
\arguments{
\item{x}{An \linkS4class{Expression}, vector, or matrix.}

\item{M}{(Optional) A positive scalar value representing the threshold. Defaults to 1.}
}
\value{
An \linkS4class{Expression} representing the Huber function evaluated at the input.
}
\description{
The elementwise Huber function, \eqn{Huber(x, M) = 1}
\describe{
  \item{\eqn{2M|x|-M^2}}{for \eqn{|x| \geq |M|}}
   \item{\eqn{|x|^2}}{for \eqn{|x| \leq |M|.}}
 }
}
\examples{
set.seed(11)
n <- 10
m <- 450
p <- 0.1    # Fraction of responses with sign flipped

# Generate problem data
beta_true <- 5*matrix(stats::rnorm(n), nrow = n)
X <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
y_true <- X \%*\% beta_true
eps <- matrix(stats::rnorm(m), nrow = m)

# Randomly flip sign of some responses
factor <- 2*rbinom(m, size = 1, prob = 1-p) - 1
y <- factor * y_true + eps

# Huber regression
beta <- Variable(n)
obj <- sum(huber(y - X \%*\% beta, 1))
prob <- Problem(Minimize(obj))
result <- solve(prob)
result$getValue(beta)
}
