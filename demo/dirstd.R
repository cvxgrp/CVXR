# Import problem data
data(dspop)   # Population
data(dssamp)  # Skewed sample

ypop <- dspop[,1]
Xpop <- dspop[,-1]
y <- dssamp[,1]
X <- dssamp[,-1]
m <- nrow(X)

# Given population mean of features
b <- as.matrix(apply(Xpop, 2, mean))

# Construct the direct standardization problem
w <- Variable(m)
objective <- sum(entr(w))
constraints <- list(w >= 0, sum(w) == 1, t(X) %*% w == b)
prob <- Problem(Maximize(objective), constraints)

# Solve for the distribution weights
result <- solve(prob)
result$value
result$getValue(w)
weights <- result$getValue(w)

# Plot probability density function
cl <- rainbow(3)
plot(density(ypop), col = cl[1], xlab = "y", ylab = NA, ylim = c(0, 0.2), zero.line = FALSE)
lines(density(y), col = cl[2])
lines(density(y, weights = weights), col = cl[3])
legend("topleft", c("True", "Sample", "Estimate"), lty = c(1,1,1), col = cl)

# Plot cumulative distribution function
plot_cdf <- function(data, probs, color = 'k', lwd = 1) {
  if(missing(probs))
    probs <- rep(1.0/length(data), length(data))
  distro <- cbind(data, probs)
  dsort <- distro[order(distro[,1]),]
  ecdf <- base::cumsum(dsort[,2])
  lines(dsort[,1], ecdf, col = color, lwd = lwd)
}

# Compare weighted with original distribution
plot(NA, main = "Cumulative Distribution Function", xlab = "y", ylab = NA, ylim = c(0,1), xlim = c(-2, 14))
plot_cdf(ypop, color = cl[1])
plot_cdf(y, color = cl[2])
plot_cdf(y, weights, color = cl[3])
legend("topleft", c("True", "Sample", "Estimate"), lty = c(1,1,1), col = cl)
