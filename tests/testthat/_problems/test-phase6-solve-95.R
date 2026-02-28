# Extracted from test-phase6-solve.R:95

# test -------------------------------------------------------------------------
x <- Variable(2)
prob <- Problem(Minimize(sum(x)), list(x >= 1))
psolve(prob, verbose = FALSE)
ss <- solver_stats(prob)
expect_true(S7::S7_inherits(ss, SolverStats))
