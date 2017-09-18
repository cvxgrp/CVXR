library(cvxr)
rho <- 1.0
n <- 3
m <- 2
##set.seed(123)
##A <- matrix(rnorm(m * n), nrow=m)
##b <- rnorm(m)

A <- matrix(c(-0.5604756, 1.55870831, 0.1292877,
              -0.2301775, 0.07050839, 1.7150650), byrow=TRUE, ncol=n)
b <- c( 0.4609162, -1.2650612)

x <- Variable(n)
f <- Pnorm(x, 1)

# Solve with CVXPY.
constraint <- list(A %*% x == b)
p <- Problem(Minimize(f), constraint)
##debug(cvxr:::solve.Problem)
##debug(cvxr:::.update_problem_state)
##debug(cvxr:::valuesById)
soln <- solve(p, verbose = TRUE)

sprintf("Optimal value from CVXPY %f", soln$value)

# Solve with method of multipliers.
resid <- A %*% x - b
y <- Parameter(rows = m)
value(y) <- matrix(0, nrow = m)

aug_lagr <- f + t(y) %*% resid + (rho / 2) * SumSquares(resid)
p <- Problem(Minimize(aug_lagr))
soln <- solve(p, verbose = TRUE)

##debug(set_matrix_data)
##data <- get_problem_data(p, "ECOS")


#sol <- solve(p)

## for (i in seq.int(MAX_ITERS)) {
##     sol <- solve(Problem(Minimize(aug_lagr)))
##     value(y) <- value(y)  + rho * resid.value
## }

## print "Optimal value from method of multipliers", f.value
