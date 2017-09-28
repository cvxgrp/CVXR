# SVM classifier with L1 regularization
n <- 20
m <- 1000
TEST <- m
DENSITY <- 0.2

# Generate data
beta_true <- matrix(rnorm(n), nrow = n, ncol = 1)
idxs <- sample(1:n, size = round((1-DENSITY)*n), replace = FALSE)
beta_true[idxs] <- 0

offset <- 0
sigma <- 45
X <- matrix(rnorm(m*n, 0, 5), nrow = m, ncol = n)
y <- sign(X %*% beta_true + offset + matrix(rnorm(m, 0, sigma), nrow = m, ncol = 1))
X_test <- matrix(rnorm(TEST*n, 0, 5), nrow = TEST, ncol = n)
y_test <- sign(X_test %*% beta_true + offset + matrix(rnorm(TEST, 0, sigma), nrow = TEST, ncol = 1))

# Form problem
beta <- Variable(n)
v <- Variable()
loss <- sum(Pos(1 - y * (X %*% beta - v)))
reg <- Norm1(beta)

# Compute a trade-off curve and record train/test error
TRIALS <- 5
train_error <- rep(0, TRIALS)
test_error <- rep(0, TRIALS)
lambda_vals <- 10^seq(-2, 0, length.out = TRIALS)
beta_vals <- matrix(0, nrow = TRIALS, ncol = n)

for(i in 1:TRIALS) {
  lambd <- lambda_vals[i]
  prob <- Problem(Minimize(loss/m + lambd*reg))
  result <- solve(prob)
  beta_res <- result$getValue(beta)
  v_res <- result$getValue(v)
  
  train_error[i] <- sum(sign(X %*% beta_true + offset) != sign(X %*% beta_res - rep(v_res, m)))/m
  test_error[i] <- sum(sign(X_test %*% beta_true + offset) != sign(X_test %*% beta_res - rep(v_res, m)))/TEST
  beta_vals[i,] <- beta_res
}

# Plot the train and test error over the trade-off curve
matplot(lambda_vals, cbind(train_error, test_error), log = "x", xlab = "Lambda", ylab = "Error", col = 1:2, type = "l", lty = 1)
legend("topleft", c("Training Error", "Test Error"), col = 1:2, lty = c(1,1))

# Plot the regularization path for beta
matplot(lambda_vals, beta_vals, log = "x", main = "Regularization Path of Beta", xlab = "Lambda", ylab = "Beta", type = "l")
