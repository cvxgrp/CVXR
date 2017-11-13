# Problem data
n <- 5
alpha <- seq(10, n-1+10)/n
beta <- seq(10, n-1+10)/n
P_tot <- 0.5
W_tot <- 1.0

# Form problem
P <- Variable(n)   # Power
W <- Variable(n)   # Bandwidth
R <- kl_div(alpha*W, alpha*(W + beta*P)) - alpha*beta*P   # Bitrate
objective <- Minimize(sum(R))
constraints <- list(P >= 0, W >= 0, sum(P) == P_tot, sum(W) == W_tot)
prob <- Problem(objective, constraints)
result <- solve(prob)

# Optimal utility, power, and bandwidth
-result$value
result$getValue(P)
result$getValue(W)
