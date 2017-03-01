library(cvxr)
library(glmnet)
library(microbenchmark)

# Problem data
n <- 1000
m <- 500
A <- matrix(rnorm(n*p), nrow = n, ncol = m)
b <- matrix(rnorm(n), nrow = n, ncol = 1)

# Construct the least-squares problem
x <- Variable(rows = m)
objective <- Minimize((1/2)*SumSquares(A %*% x - b)/n)
prob <- Problem(objective)
time_cvxr <- microbenchmark(sol_cvxr <- solve(prob, solver = ECOS()), times = 10)
xsol_cvxr <- sol_cvxr$primal_values[[1]]

# Compare results with glmnet
time_glmnet <- microbenchmark(sol_glmnet <- glmnet(A, b, family = "gaussian", lambda = 0, standardize = FALSE, intercept = FALSE), times = 10)
xsol_glmnet <- sol_glmnet$beta

cat("GLMnet Runtime (s):", mean(time_glmnet$time)/10^9)
cat("CVXR Runtime (s):", mean(time_cvxr$time)/10^9)
cat("MSE Between Results:", mean((xsol_cvxr - xsol_glmnet)^2))

