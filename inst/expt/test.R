## d <- Deque$new()
## d$push_back(1)
## d$push_back(list(x=1, y=2))
## d$pop()
## d


library(cvxr)

# Problem data.

m <- 3
n <- 2
# numpy.random.seed(1)
# A = numpy.random.randn(m, n)
# b = numpy.random.randn(m)

A <- matrix(c(1, 2, 3, 4, 2, 1), nrow=m, byrow=TRUE)
b <- matrix(c(1, 2, 3), nrow=m)

# Construct the problem.
x <- Variable(n)
objective <- Minimize(SumSquares(A*x - b))
##objective <- canonicalize(objective)
constraint <- list( 1 <= x)
prob <- Problem(objective, constraint)

debug(cvxr_solve)
debug(build_lin_op_tree)

cvxr_solve(prob)

# The optimal objective is returned by prob.solve().
result = prob.solve()
# The optimal value for x is stored in x.value.
print x.value
# The optimal Lagrange multiplier for a constraint
# is stored in constraint.dual_value.
print constraints[0].dual_value

tmp <- list()
##C_objective <- build_lin_op_tree(objective, tmp)

root_linR <- objective
