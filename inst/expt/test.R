## d <- Deque$new()
## d$push_back(1)
## d$push_back(list(x=1, y=2))
## d$pop()
## d


library(cvxr)

# Problem data.


m <- 3
n <- 2

A <- matrix(c(1, 2, 3, 4, 2, 1), nrow=m, byrow=TRUE)
b <- matrix(c(1, 2, 3), nrow=m)

# Construct the problem.
x <- Variable(n)
objective <- Minimize(SumSquares(A %*% x - b))
constraint <- list( 1 <= x)
prob <- Problem(objective, constraint)

debug(get_problem_matrix)
solve(prob)

