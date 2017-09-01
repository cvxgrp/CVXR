library(cvxr)
MAX_ITERS <- 10
rho <- 1.0
n <- 20
m <- 10
set.seed(123)
A <- matrix(rnorm(m * n), nrow=m)
b <- rnorm(m)

x <- Variable(n)
f <- Pnorm(x, 1)

# Solve with CVXPY.
p <- Problem(Minimize(f), list(A %*% x == b))
soln <- solve(p, verbose = TRUE)

sprintf("Optimal value from CVXPY %f", soln$value)

# Solve with method of multipliers.
resid <- A %*% x - b
y <- Parameter(rows = m)
value(y) <- matrix(0, nrow = m)

aug_lagr <- f + t(y) %*% resid + (rho / 2) * SumSquares(resid)
debug(set_matrix_data)

sol <- solve(Problem(Minimize(aug_lagr)))

for (i in seq.int(MAX_ITERS)) {
    sol <- solve(Problem(Minimize(aug_lagr)))
    value(y) <- value(y)  + rho * resid.value
}

print "Optimal value from method of multipliers", f.value
