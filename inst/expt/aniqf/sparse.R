library(cvxr)
library(Matrix)

# Problem data
i <- c(1,3:8)
j <- c(2,9,6:10)
x <- 7 * (1:7)
A <- sparseMatrix(i, j, x = x)     # 8 x 10 "dgCMatrix"
summary(A)

m <- nrow(A)
n <- ncol(A)

# Construct the problem
y <- Variable(n)
objective <- Minimize(SumEntries(A*y))
constraint <- list(1 <= y)
prob <- Problem(objective, constraint)

# debug(cvxr_solve)
# debug(build_lin_op_tree)

result <- cvxr_solve(prob)
result$optimal_value
result$primal_values[[y@id]]
