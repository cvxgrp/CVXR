# Skewed sampling distribution of data
skew_sample <- function(data, bias) {
  if(missing(bias))
    bias <- rep(1.0, ncol(data))
  num <- exp(data %*% bias)
  return(num / sum(num))
}

# Plot cumulative distribution function
plot_cdf <- function(data, probs, color = 'k') {
  if(missing(probs))
    probs <- rep(1.0/length(data), length(data))
  distro <- cbind(data, probs)
  dsort <- distro[order(distro[,1]),]
  ecdf <- cumsum(dsort[,2])
  lines(dsort[,1], ecdf, col = color)
}

# Problem data
n <- 2
m <- 1000
msub <- 100

# Generate original distribution
sex <- rbinom(m, 1, 0.5)
age <- sample(10:60, m, replace = TRUE)
mu <- 5 * sex + 0.1 * age
X <- cbind(sex, age)
y <- rnorm(mu, 1.0)
b <- as.matrix(apply(X, 2, mean))

# Generate skewed subsample
skew <- skew_sample(X, c(-0.95, -0.05))
sub <- sample(1:m, msub, replace = TRUE, prob = skew)

# Construct the direct standardization problem
w <- Variable(msub)
objective <- sum(Entr(w))
constraints <- list(w >= 0, sum(w) == 1, t(X[sub,]) %*% w == b)
prob <- Problem(Maximize(objective), constraints)

# Solve for the distribution weights
result <- solve(prob)
result$optimal_value
# TODO: More user-friendly functions to retrieve results
result$primal_values[[as.character(w@id)]]
weights <- result$primal_values[[as.character(w@id)]]

cl <- rainbow(3)
plot(NA, main = "Cumulative Distribution Function", xlab = "y", ylab = NA, xlim = c(-2, 12), ylim = c(0, 1))
plot_cdf(y, color = cl[1])
plot_cdf(y[sub], color = cl[2])
plot_cdf(y[sub], weights, color = cl[3])
legend("topleft", c("True", "Sample", "Estimate"), lty = c(1,1,1), col = cl)
