library(cvxr)

n <- 2
A <- diag(c(1, 1))
x <- Variable(n)
objective <- Minimize(sum(A %*% x))
constraint <- list(x >=0, x <= 1)
prob <- Problem(objective, constraint)

#debug(get_problem_matrix)
solve(prob, solver = "SCS", verbose = TRUE)

x <- Variable(1)
obj <- Maximize(Log(x))
solve(Problem(obj))


