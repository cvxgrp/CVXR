library(cvxr)

set.seed(12345)

# Initialize data.
MAX_ITERS = 10
rho = 1.0
n = 20
m = 10
A = matrix(rnorm(m * n), nrow = m)
b = matrix(rnorm(m), nrow = m)

# Initialize problem.
x = Variable(n)
f = norm(x, 1)

# Solve with CVXPY.
Problem(Minimize(f), [A*x == b]).solve()
print "Optimal value from CVXPY", f.value

# Solve with method of multipliers.
resid = A*x - b
y = Parameter(m); y.value = np.zeros(m)
aug_lagr = f + y.T*resid + (rho/2)*sum_squares(resid)
for t in range(MAX_ITERS):
    Problem(Minimize(aug_lagr)).solve()
    y.value += rho*resid.value

print "Optimal value from method of multipliers", f.value

                                        # Problem data.

m <- 1
n <- 2

A <- matrix(c(17, 19), nrow=m, byrow=TRUE)

# Construct the problem.
x <- Variable(n)
objective <- Minimize(A*x)
constraint <- list(1 <= x)

##base::trace("canonicalize", tracer=browser, exit = browser, signature = c("Variable"))
prob <- Problem(objective, constraint)

##base::trace("cvxr_solve", tracer=browser, exit = browser, signature = c("Problem"))

##debug(cvxr_solve)

##debug(build_lin_op_tree)


result <- cvxr_solve(prob)

cat("Solver Status: ", result$status, "\n")

cat("Primal Solution: ")
print(result$primal_values)

cat("Dual Solution: ")
print(result$dual_values)


