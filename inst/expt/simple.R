library(cvxr)

# Problem data.

m <- 1
n <- 2

A <- matrix(c(17, 19), nrow=m, byrow=TRUE)

# Construct the problem.
x <- Variable(n)
objective <- Minimize(A*x)
constraint <- list(1 <= x)
prob <- Problem(objective, constraint)

base::trace("cvxr_solve", tracer=browser, exit = browser, signature = c("Problem"))

##debug(cvxr_solve)

##debug(build_lin_op_tree)

cvxr_solve(prob)


