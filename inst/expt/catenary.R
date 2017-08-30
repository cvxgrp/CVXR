## This is a working example

library(cvxr)

# Generate data
n <- 51
L <- 2
h <- L/(n-1)

# Define objective
x <- Variable(2*n)
B <- diag(2*n)
B[1:n,1:n] <- 0
objective <- Minimize(sum(B %*% x))

# Define constraints
A <- matrix(0, nrow = 4, ncol = 2*n)
A[1, 1] <- A[2, n] <- A[3, n+1] <- A[4, 2*n] <- 1
b <- matrix(c(0, 1, 1, 1), nrow = 4)
constraints <- list(x >= 0, A %*% x == b)

for(i in seq.int(n-1)) {
    A <- matrix(numeric(2 * 2 * n), nrow = 2)
    A[1, i] <- -1; A[1, i+1] <- 1
    A[2, n+i] <- -1; A[2, n+i+1] <- 1
    constraints <- c(constraints, Norm2(A %*% x) <= h)
}
# TODO: Fix SOC_AXIS LinOp type invalid error
# constraints <- list(diff(x[1:n])^2 + diff(x[(n+1):(2*n)])^2 <= h^2)

# Solve problem
prob <- Problem(objective, constraints)
system.time(sole <- solve(prob))
x <- matrix(sole$primal, ncol = 1)


# Plot and compare with ideal catenary
xs <- x[1:n, 1, drop = TRUE]
ys <- x[(n+1):(2*n), 1, drop = TRUE]
plot(c(0,1), c(0,1), type = "n", xlab = "x", ylab = "y")
lines(xs, ys, col = "blue", lwd = 2)

points(c(0,1), c(1,1))
curve(0.22964*cosh((x - 0.5)/0.22964) - 0.02603, 0, 1, col = "red", add = TRUE)
grid()
