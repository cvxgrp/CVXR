# Binary support vector classifier
n <- 2
m <- 50

# Generate data
X <- matrix(rnorm(m*n), nrow = m, ncol = n)
y <- c(rep(-1, m/2), rep(1, m/2))
X[y == 1,] = X[y == 1,] + 1

# Form problem
cost <- 10
beta <- Variable(n)
slack <- Variable(m)
loss <- (1/2)*sum(beta) + cost*sum(slack)
constr <- list(y * (X %*% beta) >= 1 - slack, slack >= 0)
prob <- Problem(Minimize(loss), constr)

# Solve and plot result
result <- solve(prob)
result$optimal_value
# TODO: More user-friendly functions to retrieve results
bres <- result$primal_values[[as.character(beta@id)]]
sres <- result$primal_values[[as.character(slack@id)]]

plot(X, col = (3-y), main = "Support Vector Classifier")
legend("topleft", paste("y =", unique(y)), col = 3-unique(y), lty = 1)
abline(a = 1/bres[2], b = -bres[1]/bres[2])