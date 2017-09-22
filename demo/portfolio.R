# Problem data
set.seed(10)
n <- 10
mu <- matrix(abs(rnorm(n)), nrow = n)
Sigma <- matrix(rnorm(n^2), nrow = n, ncol = n)
Sigma <- t(Sigma) %*% Sigma

# Form problem
w <- Variable(n)
ret <- t(mu) %*% w
risk <- QuadForm(w, Sigma)
constraints <- list(w >= 0, sum(w) == 1)

# Risk aversion parameters
SAMPLES <- 100
gammas <- 10^seq(-2, 3, length.out = SAMPLES)
ret_data <- rep(0, SAMPLES)
risk_data <- rep(0, SAMPLES)

# Compute trade-off curve
for(i in 1:length(gammas)) {
  gamma <- gammas[i]
  objective <- ret - gamma * risk
  prob <- Problem(Maximize(objective), constraints)
  result <- solve(prob)
  
  # Evaluate risk/return for current solution
  risk_data[i] <- result$getValue(sqrt(risk))
  ret_data[i] <- result$getValue(ret)
}

# Plot trade-off curve
plot(risk_data, ret_data, xlab = "Risk (Standard Deviation)", ylab = "Return", xlim = c(0, 4), ylim = c(0, 2), type = "l", lwd = 2, col = "blue")
points(sqrt(diag(Sigma)), mu, col = "red", cex = 1.5, pch = 16)
markers_on <- c(30, 41)
for(marker in markers_on) {
  points(risk_data[marker], ret_data[marker], col = "black", cex = 1.5, pch = 15)
  nstr <- sprintf("%.2f", gammas[marker])
  text(risk_data[marker] + 0.18, ret_data[marker] - 0.03, bquote(paste(gamma, " = ", .(nstr))))
}