# Extracted from test-phase6-solve.R:137

# test -------------------------------------------------------------------------
x <- Variable(2)
prob <- Problem(Minimize(sum(x)), list(x >= 1))
psolve(prob, verbose = FALSE)
sol <- solution(prob)
expect_true(S7::S7_inherits(sol, Solution))
