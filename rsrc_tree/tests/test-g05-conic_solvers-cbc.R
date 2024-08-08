context("test-g05-conic_solvers-cbc")
TOL <- 1e-6

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test that all CBC solver options work", {
  prob <- Problem(Minimize(norm2(x)), list(x == Variable(2, boolean = TRUE)))
  if("CBC" %in% installed_solvers()) {
    for(i in 1:2) {
      # Some cut-generators seem to be buggy for now -> set to false
      # prob.solve(solver=cvx.CBC, verbose=True, GomoryCuts=True, MIRCuts=True,
      #            MIRCuts2=True, TwoMIRCuts=True, ResidualCapacityCuts=True,
      #            KnapsackCuts=True, FlowCoverCuts=True, CliqueCuts=True,
      #            LiftProjectCuts=True, AllDifferentCuts=False, OddHoleCuts=True,
      #            RedSplitCuts=False, LandPCuts=False, PreProcessCuts=False,
      #            ProbingCuts=True, SimpleRoundingCuts=True)
      result <- solve(prob, solver = "CBC", verbose = TRUE, maximumSeconds = 100)
    }
    expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
  } else {
    expect_error(solve(prob, solver = "CBC"), 
                 "The solver CBC is not installed", fixed = TRUE)
  }
})

test_that("test cbc lp 0", {
  StandardTestLPs.test_lp_0(solver = "CBC", duals = FALSE)
})

test_that("test cbc lp 1", {
  StandardTestLPs.test_lp_1(solver = "CBC", duals = FALSE)
})

test_that("test cbc lp 2", {
  StandardTestLPs.test_lp_2(solver = "CBC", duals = FALSE)
})

test_that("test cbc lp 3", {
  StandardTestLPs.test_lp_3(solver = "CBC")
})

test_that("test cbc lp 4", {
  StandardTestLPs.test_lp_4(solver = "CBC")
})

test_that("test cbc lp 5", {
  StandardTestLPs.test_lp_5(solver = "CBC", duals = FALSE)
})

test_that("test cbc mi lp 0", {
  StandardTestLPs.test_mi_lp_0(solver = "CBC")
})

test_that("test cbc mi lp 1", {
  StandardTestLPs.test_mi_lp_1(solver = "CBC")
})

test_that("test cbc mi lp 2", {
  StandardTestLPs.test_mi_lp_2(solver = "CBC")
})

test_that("test cbc mi lp 3", {
  StandardTestLPs.test_mi_lp_3(solver = "CBC")
})

test_that("test cbc mi lp 5", {
  # TODO: CyLP <= 0.92.4 has no working integer infeasibility detection, so
  # if detected, this test should be skipped. See CVXPY test_conic_solvers.py under CBC tests.
  # https://github.com/coin-or/CyLP/pull/150
  StandardTestLPs.test_mi_lp_5(solver = "CBC")
})
