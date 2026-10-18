## Claims audit, decision 27 follow-up (CVXR 1.9.2.9166): ledger row C0861,
## the Gurobi iteration count.
##
##   Both Gurobi interfaces set NUM_ITERS unconditionally, zero included
##   (gurobi_qpif.py:97-99, gurobi_conif.py:112-116), and pass it to
##   failure_solution too. CVXR set it only when the count was positive, so a
##   problem Gurobi's presolve finished reported NA where CVXPY reports 0.
##   Measured 2026-10-11 on Gurobi 13.0.3 from both sides (notes/audit/
##   claims_audit_item1.md, "Decision 27 follow-up 3").
##
##   Kept, recorded divergence: on the conic interface upstream reads
##   BarIterCount alone (gurobi_conif.py:112), so an LP or MIP that Gurobi
##   solves by simplex reports 0 there (measured: 4 dual-simplex iterations,
##   num_iters 0). CVXR sums barrier and simplex counts on both interfaces,
##   the rule upstream itself applies in gurobi_qpif.py:92-94 and CVXR applies
##   for MOSEK and HiGHS (completeness_known.tsv, gurobi_conif.py NUM_ITERS).
##   The third test pins that sum. It is a pin of existing behavior and cannot
##   be seen failing.
##
## The first two tests were seen failing on 1.9.2.9165 first.

# -- C0861: the count is reported even when it is zero --------------------

## @cvxpy NONE
test_that("Gurobi QP interface: num_iters is 0, not NULL, when presolve solves the problem", {
  skip_if_not_installed("gurobi")
  x <- Variable(1)
  prob <- Problem(Minimize(x), list(x >= 1))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, 1, tolerance = 1e-8)
  ss <- solver_stats(prob)
  es <- ss@extra_stats
  expect_equal(es$itercount + es$baritercount, 0)
  expect_false(is.null(ss@num_iters))
  expect_true(is.numeric(ss@num_iters))
  expect_equal(ss@num_iters, 0)
})

## @cvxpy NONE
test_that("Gurobi conic interface: num_iters is 0, not NULL, when presolve solves the problem", {
  skip_if_not_installed("gurobi")
  b <- Variable(2, boolean = TRUE)
  prob <- Problem(Maximize(sum(b)), list(b[1] + b[2] <= 1))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, 1, tolerance = 1e-8)
  ss <- solver_stats(prob)
  es <- ss@extra_stats
  expect_equal(es$itercount + es$baritercount, 0)
  expect_false(is.null(ss@num_iters))
  expect_true(is.numeric(ss@num_iters))
  expect_equal(ss@num_iters, 0)
})

# -- Recorded divergence: the conic interface sums barrier and simplex ----

## @cvxpy NONE
test_that("Gurobi conic interface: num_iters is the barrier plus simplex sum (pin, recorded divergence)", {
  skip_if_not_installed("gurobi")
  ## A small knapsack with presolve off is solved by simplex at the nodes:
  ## BarIterCount stays 0 and IterCount is positive, the case upstream
  ## reports as 0 iterations.
  w <- c(3, 4, 5, 6)
  v <- c(2, 3, 4, 5)
  b <- Variable(4, boolean = TRUE)
  prob <- Problem(Maximize(sum(v * b)), list(sum(w * b) <= 9))
  psolve(prob, solver = "GUROBI", Presolve = 0)
  expect_equal(status(prob), "optimal")
  ss <- solver_stats(prob)
  es <- ss@extra_stats
  expect_gt(es$itercount, 0)
  expect_equal(ss@num_iters, es$itercount + es$baritercount)
})
