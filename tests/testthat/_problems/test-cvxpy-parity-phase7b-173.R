# Extracted from test-cvxpy-parity-phase7b.R:173

# test -------------------------------------------------------------------------
set.seed(0)
n <- 3L
alpha <- 0.275
x <- Variable(n)
y <- Variable(n)
z <- Variable(n)
con <- PowCone3D(x, y, z, alpha)
x0 <- 0.1 + runif(n)
y0 <- 0.1 + runif(n)
z0 <- x0^alpha * y0^(1 - alpha)
z0[2] <- -z0[2]
save_leaf_value(x, matrix(x0, n, 1))
