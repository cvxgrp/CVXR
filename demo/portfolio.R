# Problem data
n <- 10
SAMPLES <- 100
mu <- matrix(abs(rnorm(n)), nrow = n)
Sigma <- matrix(rnorm(n^2), nrow = n, ncol = n)
Sigma <- t(Sigma) %*% Sigma

# Form problem
w <- Variable(n)
ret <- t(mu) %*% w
risk <- QuadForm(w, Sigma)
constraints <- list(w >= 0, sum(w) == 1)

# Risk aversion parameters
gammas <- 10^seq(-2, 3, length.out = SAMPLES)
ret_data <- rep(0, SAMPLES)
risk_data <- rep(0, SAMPLES)
w_data <- matrix(0, nrow = SAMPLES, ncol = n)

# Compute trade-off curve
for(i in 1:length(gammas)) {
  gamma <- gammas[i]
  objective <- ret - gamma * risk
  prob <- Problem(Maximize(objective), constraints)
  result <- solve(prob)
  
  # Evaluate risk/return for current solution
  risk_data[i] <- result$getValue(sqrt(risk))
  ret_data[i] <- result$getValue(ret)
  w_data[i,] <- result$getValue(w)
}

# Plot trade-off curve
plot(risk_data, ret_data, xlab = "Risk (Standard Deviation)", ylab = "Return", xlim = c(0, 4), ylim = c(0, 2), type = "l", lwd = 2, col = "blue")
points(sqrt(diag(Sigma)), mu, col = "red", cex = 1.5, pch = 16)
markers_on <- c(10, 20, 30, 40)
for(marker in markers_on) {
  points(risk_data[marker], ret_data[marker], col = "black", cex = 1.5, pch = 15)
  nstr <- sprintf("%.2f", gammas[marker])
  text(risk_data[marker] + 0.18, ret_data[marker] - 0.03, bquote(paste(gamma, " = ", .(nstr))))
}

# Plot weights for a few gamma
w_plot <- t(w_data[markers_on,])
colnames(w_plot) <- sprintf("%.2f", gammas[markers_on])
if("colorspace" %in% rownames(installed.packages())) {
  require(colorspace)
  barplot(w_plot, xlab = expression(paste("Risk Aversion (", gamma, ")", sep = "")), ylab = "Fraction of Budget", col = sequential_hcl(n))
} else
  barplot(w_plot, xlab = expression(paste("Risk Aversion (", gamma, ")", sep = "")), ylab = "Fraction of Budget")