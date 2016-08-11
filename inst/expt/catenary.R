n = 51
L = 2; h = L/(n-1)

x <- Variable(2*n)
B <- diag(2*n)
B[1:n, 1:n] <- 0
objective <- Minimize(SumEntries(B*x))

constraints <- list(x > 0,
                    x[1] = 0,
                    x[n] = 1,
                    x[n+1] = 1,
                    x[2*n] = 1)

lapply(1:(n-1),
       function(i) {
           A <- matrix(rep(0, 2 * 2 * n), nrow = 2)
           A[1, i] = -1; A[1, i+1] = 1
           A[2, n+i] = -1; A[2, n+i+1] = 1
           norm(A*x) <= h
       })

status = solve(m)


