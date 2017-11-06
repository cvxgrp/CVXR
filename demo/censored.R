# Problem data
n <- 30
M <- 50
K <- 200

set.seed(n*M*K)
X <- matrix(rnorm(K*n), nrow = K, ncol = n)
beta_true <- matrix(rnorm(n), nrow = n, ncol = 1)
y <- X %*% beta_true + 0.3*sqrt(n)*rnorm(K)

# Order variables based on y
idx <- order(y, decreasing = FALSE)
y_ordered <- y[idx]
X_ordered <- X[idx,]

# Find cutoff and censor
D <- (y_ordered[M] + y_ordered[M+1])/2
y_censored <- pmin(y_ordered, D)

plot_results <- function(beta_res, bcol = "blue", bpch = 17) {
  plot(1:M, y_censored[1:M], col = "black", xlab = "Observations", ylab = "y", ylim = c(min(y), max(y)), xlim = c(1,K))
  points((M+1):K, y_ordered[(M+1):K], col = "red")
  points(1:K, X_ordered %*% beta_res, col = bcol, pch = bpch)
  abline(a = D, b = 0, col = "black", lty = "dashed")
  legend("topleft", c("Uncensored", "Censored", "Estimate"), col = c("black", "red", bcol), pch = c(1,1,bpch))
}

# Regular OLS
beta <- Variable(n)
obj <- sum((y_censored - X_ordered %*% beta)^2)
prob <- Problem(Minimize(obj))
result <- solve(prob)
beta_ols <- result$getValue(beta)
plot_results(beta_ols)

# OLS using uncensored data
obj <- sum((y_censored[1:M] - X_ordered[1:M,] %*% beta)^2)
prob <- Problem(Minimize(obj))
result <- solve(prob)
beta_unc <- result$getValue(beta)
plot_results(beta_unc)

# Censored regression
obj <- sum((y_censored[1:M] - X_ordered[1:M,] %*% beta)^2)
constr <- list(X_ordered[(M+1):K,] %*% beta >= D)
prob <- Problem(Minimize(obj), constr)
result <- solve(prob)
beta_cens <- result$getValue(beta)
plot_results(beta_cens)
