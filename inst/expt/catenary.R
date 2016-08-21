library(cvxr)

n <- 51
L <- 2
h <- L/(n-1)

x <- Variable(2*n)
B <- diag(2*n)
B[1:n, 1:n] <- 0
objective <- Minimize(SumEntries(B*x))

A <- matrix(0, nrow=4, ncol=2*n)
A[1, 1] <- A[2, n] <- A[3, n + 1] <- A[4, 2 * n] <- 1
b <- matrix(c(0, 1, 1, 1), nrow=4)

constraints = list( x >= 0, A * x == b )

for (i in seq.int(n-1)) {
    A <- matrix(numeric(2 * 2 * n), nrow = 2)
    A[1, i] <- -1; A[1, i+1] <- 1
    A[2, n+i] <- -1; A[2, n+i+1] <- 1
    constraints <- c(constraints, Norm2(A*x) <= h)
}
prob <- Problem(objective, constraints)
system.time(sole <- cvxr_solve(prob))

x <- sole$primal_values[[1]]
xs <- x[1:n, 1, drop=TRUE]; ys <- x[(n+1):(2*n), 1, drop=TRUE]
plot(c(0, 1), c(0, 1), type='n')
lines(xs, ys, col="blue", lwd=2)

points(c(0, 1), c(1, 1))
curve(0.22964*cosh((x-0.5)/0.22964)-0.02603, 0, 1,
      col="red", add=TRUE)
grid()

