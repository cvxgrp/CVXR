# Binary support vector classifier
set.seed(10)
n <- 2
m <- 50

# Generate data
X <- matrix(rnorm(m*n), nrow = m, ncol = n)
y <- c(rep(-1, m/2), rep(1, m/2))
X[y == 1,] = X[y == 1,] + 1

# Form problem
cost <- 10
beta0 <- Variable()
beta <- Variable(n)
slack <- Variable(m)
loss <- (1/2)*sum_squares(vstack(beta, beta0)) + cost*sum(slack)
constr <- list(y * (X %*% beta + beta0) >= 1 - slack, slack >= 0)
prob <- Problem(Minimize(loss), constr)

# Solve and plot result
result <- solve(prob)
result$value
b0 <- result$getValue(beta0)
bres <- result$getValue(beta)

plot(X, col = (3-y), main = "Support Vector Classifier", pch = 19)
legend("topright", paste("y =", unique(y)), col = 3-unique(y), pch = 19)
slope <- -bres[1]/bres[2]
intercept <- b0/bres[2]
margin <- 1/bres[2]
abline(a = intercept, b = slope)
abline(a = intercept + margin, b = slope, lty = 2)
abline(a = intercept - margin, b = slope, lty = 2)
