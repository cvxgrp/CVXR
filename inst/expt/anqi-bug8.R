# This throws an error: Subscript out of bounds in saveValuesById
library(cvxr)
a <- Variable(1)
obj <- Minimize(0*a)
p <- Problem(obj)
result <- solve(p)

# This runs, but throws a warning: A->p (column pointers) not strictly increasing, column 24 empty.
x <- Variable(2,2)
obj <- Minimize(TotalVariation(x))
prob <- Problem(obj)
base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("SCS"))
result <- solve(prob, solver = "SCS")
