set.seed(1)
n <- 20      # Total bets
K <- 100     # Number of possible returns
PERIODS <- 100
TRIALS <- 5

# Generate return probabilities
ps <- runif(K)
ps <- ps/sum(ps)

# Generate matrix of possible returns
rets <- runif(K*(n-1), 0.5, 1.5)
shuff <- sample(1:length(rets), size = length(rets), replace = FALSE)
rets[shuff[1:30]] <- 0      # Set 30 returns to be relatively low
rets[shuff[31:60]] <- 5     # Set 30 returns to be relatively high
rets <- matrix(rets, nrow = K, ncol = n-1)
rets <- cbind(rets, rep(1, K))   # Last column represents not betting

# Solve for Kelly optimal bets
b <- Variable(n)
obj <- Maximize(t(ps) %*% log(rets %*% b))
constraints <- list(sum(b) == 1, b >= 0)
prob <- Problem(obj, constraints)
result <- solve(prob)
bets <- result$getValue(b)
bets

# Naive betting scheme: bet in proportion to expected return
bets_cmp <- matrix(0, nrow = n)
bets_cmp[n] <- 0.15                  # Hold 15% of wealth
rets_avg <- ps %*% rets
# tidx <- order(rets_avg[-n], decreasing = TRUE)[1:9]
tidx <- 1:(n-1)
fracs <- rets_avg[tidx]/sum(rets_avg[tidx])
bets_cmp[tidx] <- fracs*(1-bets_cmp[n])

# Calculate wealth over time
wealth <- matrix(0, nrow = PERIODS, ncol = TRIALS)
wealth_cmp <- matrix(0, nrow = PERIODS, ncol = TRIALS)
for(i in 1:TRIALS) {
  sidx <- sample(1:K, size = PERIODS, replace = TRUE, prob = ps)
  winnings <- rets[sidx,] %*% bets
  wealth[,i] <- cumprod(winnings)
  
  winnings_cmp <- rets[sidx,] %*% bets_cmp
  wealth_cmp[,i] <- cumprod(winnings_cmp)
}

# Plot Kelly optimal growth trajectories
matplot(1:PERIODS, wealth, xlab = "Time", ylab = "Wealth", log = "y", type = "l", col = "red", lty = 1, lwd = 2)
matlines(1:PERIODS, wealth_cmp, col = "blue", lty = 2, lwd = 2)
legend("topleft", c("Kelly Optimal Bets", "Naive Bets"), col = c("red", "blue"), lty = c(1, 2), lwd = 2, bty = "n")
