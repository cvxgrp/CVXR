library(cvxr)
n <- 10
set.seed(123)
mu <- abs(rnorm(n))
Sigma <- matrix(rnorm(n*n), nrow=n)
Sigma <- t(Sigma) %*% Sigma

# Long only portfolio optimization.
w <- Variable(n)
ret <- t(mu) %*% w
risk <- quad_form(w, Sigma)

constraints <- list(SumEntries(w) == 1,
                    w >= 0)

prob <- Problem(Maximize(ret - risk),
                constraints)
soln <- solve(prob, verbose = TRUE)

soln$getDualValue(constraints[[1]])
soln$getValue(w)

soln$getDualValue(constraints[[2]])

