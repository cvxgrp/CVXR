# Problem data
w <- matrix(c(0.1, 0.2, -0.05, 0.1))
# Constraint matrix:
# [[0.2, + ,  +,   +-],
#  [+,   0.1, -,   - ],
#  [+,   -,   0.3, + ],
#  [+-,  -,   +,  0.1]]

# Form problem
Sigma <- Semidef(4)
obj <- Maximize(t(w) %*% Sigma %*% w)
constraints <- list(Sigma[1,1] == 0.2, Sigma[2,2] == 0.1, Sigma[3,3] == 0.3, Sigma[4,4] == 0.1,
                    Sigma[1,2] >= 0, Sigma[1,3] >= 0, Sigma[2,3] <= 0, Sigma[2,4] <= 0, Sigma[3,4] >= 0)
prob <- Problem(obj, constraints)
result <- solve(prob, solver = "SCS")
result$getValue(Sigma)