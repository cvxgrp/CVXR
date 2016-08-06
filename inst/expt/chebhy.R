## An example that currently gives wrong answers, for instance!

library(cvxr)
## "Find the largest Euclidean ball in the polyhedron",
## Create the input data


a1 <- matrix(c(2,1))
a2 <- matrix(c(2,-1))
a3 <- matrix(c(-1,2))
a4 <- matrix(c(-1,-2))
b <- rep(1,4)

## Create and solve the model
r <- Variable(name = "r")
x_c <- Variable(2, name = "x_c")
obj <- Maximize(r)
constraints <- list(
    t(a1)*x_c + norm(a1,"F")*r <= b[1],
    t(a2)*x_c + norm(a2,"F")*r <= b[2],
    t(a3)*x_c + norm(a3,"F")*r <= b[3],
    t(a4)*x_c + norm(a4,"F")*r <= b[4]
)

p <- Problem(obj, constraints)
##base::trace("cvxr_solve", tracer=browser, exit = browser, signature = c("Problem"))
result <- cvxr_solve(p)
result$optimal_value
result$primal_values
