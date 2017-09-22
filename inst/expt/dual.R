library(cvxr)


# Create two scalar optimization variables.
x <- Variable()
y <- Variable()

# Create two constraints.
constraints <- list(x + y == 1,
                    x - y >= 1)

# Form objective.
obj <- Minimize(Square(x - y))

# Form and solve problem.
prob <- Problem(obj, constraints)
result <- solve(prob, verbose = TRUE)

# The optimal dual variable (Lagrange multiplier) for
# a constraint is stored in constraint.dual_value.
print "optimal (x + y == 1) dual variable", constraints[0].dual_value
print "optimal (x - y >= 1) dual variable", constraints[1].dual_value
print "x - y value:", (x - y).value
optimal (x + y == 1) dual variable 6.47610300459e-18
optimal (x - y >= 1) dual variable 2.00025244976
x - y value: 0.999999986374
