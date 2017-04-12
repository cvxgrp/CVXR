library(cvxr)
library(Matrix)
library(glmnet)
library(microbenchmark)

# Problem data
n <- 1000
m <- 500
mu <- 10
sigma <- 10

A <- matrix(rnorm(n*m, mu, sigma), nrow = n, ncol = m)
x_true <- matrix(seq(0, 100, length.out = m), nrow = m)
eps <- matrix(rnorm(n, 0, 1), nrow = n)
b <- A %*% x_true + eps

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

# Sparse problem data
Asp <- rsparsematrix(n, m, density = 0.15)
bsp <- Asp %*% x_true + eps

# Construct the least-squares problem
xsp <- Variable(rows = m)
objective <- Minimize((1/2)*SumSquares(Asp %*% xsp - bsp)/n)
prob <- Problem(objective)
time_cvxr_sp <- microbenchmark(sol_cvxr_sp <- solve(prob, solver = ECOS()), times = 10)
xsol_cvxr_sp <- sol_cvxr_sp$primal_values[[1]]

# Compare results with glmnet
time_glmnet_sp <- microbenchmark(sol_glmnet_sp <- glmnet(Asp, bsp, family = "gaussian", lambda = 0, standardize = FALSE, intercept = FALSE), times = 10)
xsol_glmnet_sp <- sol_glmnet_sp$beta

cat("GLMnet Sparse Runtime (s):", mean(time_glmnet_sp$time)/10^9)
cat("CVXR Sparse Runtime (s):", mean(time_cvxr_sp$time)/10^9)
cat("MSE Between Results:", mean((xsol_cvxr_sp - xsol_glmnet_sp)^2))
