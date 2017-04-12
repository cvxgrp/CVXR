# Import problem data
data(dspop)   # Population
data(dssamp)  # Skewed sample

y <- dspop[,1]
X <- dspop[,-1]
ysam <- dssamp[,1]
Xsam <- dssamp[,-1]

# Given population mean of features
b <- as.matrix(apply(X, 2, mean))

# Construct the direct standardization problem
w <- Variable(nrow(Xsam))
objective <- sum(Entr(w))
constraints <- list(w >= 0, sum(w) == 1, t(Xsam) %*% w == b)
prob <- Problem(Maximize(objective), constraints)

# Solve for the distribution weights
result <- solve(prob)
result$optimal_value
# TODO: More user-friendly functions to retrieve results
result$primal_values[[as.character(w@id)]]
weights <- result$primal_values[[as.character(w@id)]]

# Plot cumulative distribution function
plot_cdf <- function(data, probs, color = 'k') {
  if(missing(probs))
    probs <- rep(1.0/length(data), length(data))
  distro <- cbind(data, probs)
  dsort <- distro[order(distro[,1]),]
  ecdf <- cumsum(dsort[,2])
  lines(dsort[,1], ecdf, col = color)
}

# Compare weighted with original distribution
cl <- rainbow(3)
plot(NA, main = "Cumulative Distribution Function", xlab = "y", ylab = NA, xlim = c(-2, 12), ylim = c(0, 1))
plot_cdf(y, color = cl[1])
plot_cdf(ysam, color = cl[2])
plot_cdf(ysam, weights, color = cl[3])
legend("topleft", c("True", "Sample", "Estimate"), lty = c(1,1,1), col = cl)
