library(cvxr)

n <- 51
L <- 2
h <- L/(n-1)

x <- Variable(2*n)
B <- diag(2*n)
B[1:n, 1:n] <- 0
objective <- Minimize(SumEntries(B*x))

A <- matrix(0, nrow=4, ncol=2*n)
A[1, 1] <- 0
A[2, n] <- 1
A[3, n + 1] <- 1
A[4, 2 * n] <- 1
b <- c(0, 1, 1, 1)

constraints <- list(A*x=b)


lapply(1:(n-1),
       function(i) {
           A <- matrix(rep(0, 2 * 2 * n), nrow = 2)
           A[1, i] = -1; A[1, i+1] = 1
           A[2, n+i] = -1; A[2, n+i+1] = 1
           norm(A*x) <= h
       })

prob <- Problem(objective, constraints)
status <- cvxr_solve(prob)


