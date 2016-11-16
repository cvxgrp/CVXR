library(cvxr)

# Problem data.

m <- 1
n <- 2

A <- matrix(c(17, 19), nrow=m, byrow=TRUE)

# Construct the problem.
x <- Variable(n)
objective <- Minimize(A %*% x)
constraint <- list(1 <= x)

##base::trace("canonicalize", tracer=browser, exit = browser, signature = c("Variable"))
prob <- Problem(objective, constraint)

##base::trace("solve", tracer=browser, exit = browser, signature = c("Problem"))

##debug(solve)

##debug(build_lin_op_tree)

result <- solve(prob)

cat("Solver Status: ", result$status, "\n")

cat("Primal Solution:\n")
print(result$primal_values)

cat("Dual Solution:\n")
print(result$dual_values)


