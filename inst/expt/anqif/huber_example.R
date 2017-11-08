library(cvxr)

n <- 300
SAMPLES <- round(1.5*n)
beta_true <- 5*matrix(rnorm(n), nrow = n, ncol = 1)
X <- matrix(rnorm(n*SAMPLES), nrow = n, ncol = SAMPLES)
Y <- matrix(rep(0, SAMPLES), nrow = SAMPLES, ncol = 1)
v <- matrix(rnorm(SAMPLES), nrow = SAMPLES, ncol = 1)

# Generate the sign changes
p <- 0.10
factor <- 2*rbinom(SAMPLES, 1, 1-p) - 1
Y <- factor * t(X) %*% beta_true + v

# Form and solve the Huber regression problem
beta <- Variable(n)
cost <- SumEntries(Huber(t(X) %*% beta - Y, 1))
prob <- Problem(Minimize(cost))
result <- solve(prob)
result$optimal_value
result$primal_values[[as.character(beta@id)]]
