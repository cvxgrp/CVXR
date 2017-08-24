library(cvxr)

set.seed(1234)
x <- matrix(rnorm(100), nrow=25)

S <- cov(x)

delta <- 0.1

alpha <- Variable(rows = 4)
e1 <- c(1, 0 , 0, 0)
e2 <- c(0, 1 , 0, 0)
e3 <- c(0, 0 , 1, 0)
e4 <- c(0, 0 , 0, 1)

objective <- Minimize(0.5 * quad_form(x = alpha, P = S)  + t(alpha) %*% e2 + delta * Pnorm(alpha, 1))
p <- Problem(objective)
solution <- solve(p)
