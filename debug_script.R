library(cvxr)

# TEST
# A <- Variable(2, 2, name = "A")
# obj <- Minimize(0)
# dom <- domain(LogDet(A))
# prob <- Problem(obj, dom)
 
# base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("SCS"))
# base::trace("Solver.get_problem_data", tracer = browser, exit = browser, signature = c("SCS"))
# base::trace("format_constr", tracer = browser, exit = browser, signature = c("SDP"))
# base::trace("get_objective", tracer = browser, exit = browser, signature = c("MatrixData"))
# debug(get_problem_matrix)
# debug("SymData.get_var_offsets")
# result <- solve(prob, solver = "SCS")

# TEST: test_nonlinear_atoms.R
# LinOp data field contains Parameter object rather than its value
# v <- Variable(1)
# p <- Parameter(1, sign = "positive")
# value(p) <- 1
# # p <- 1
#  
# obj <- Minimize(KLDiv(v, p))
# prob <- Problem(obj)
# result <- solve(prob)

# TEST: test_grad.R
a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

value(A) <- diag(rep(1, 2))
value(B) <- diag(rep(1, 2))
expr <- MatrixFrac(A, B)

# base::trace("grad", tracer = browser, exit = browser, signature = c("Atom"))
# base::trace(cvxr:::.axis_grad, tracer = browser, exit = browser, signature = c("AxisAtom"))
# base::trace(cvxr:::.column_grad, tracer = browser, exit = browser, signature = c("LogSumExp"))
grad(expr)
