library(cvxr)
library(Matrix)

# Problem data
i <- c(1,3:8)
j <- c(2,9,6:10)
x <- 7 * (1:7)
A <- sparseMatrix(i, j, x = x)     # 8 x 10 "dgCMatrix"
aMat <- Constant(A)
canonicalize(aMat)

summary(A)

m <- nrow(A)
n <- ncol(A)

# Construct the problem
y <- Variable(n)
objective <- Minimize(sum(A %*% y))
constraint <- list(1 <= y)
prob <- Problem(objective, constraint)

# debug(solve)
# debug(build_lin_op_tree)
base::trace("solve", tracer=browser, exit = browser, signature = c("Problem"))
result <- solve(prob)

result$optimal_value
result$primal_values[[as.character(y@id)]]
