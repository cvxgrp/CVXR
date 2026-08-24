## @cvxpy tests/test_dqcp.py
##
## DQCP with a NESTED inverse chain: the bisection's sublevel bound must keep
## depending on the bisection parameter.
##
## Why this file exists: CVXR silently returned a wrong answer -- with
## `status == "optimal"` -- for any DQCP problem whose inverse chain nests two
## or more invertible atoms. `maximize`/`minimize` of e.g.
## `sqrt(inv_pos(power(x, 2)))` inverts to a bound whose left-hand side is a
## function of the bisection PARAMETER alone, here `1/sqrt(1+t^2) <= x`. That
## expression is constant once t is fixed, and CVXPY refuses to canonicalize it
## (`canonicalize_tree(constr, canonicalize_params=False)`, dqcp2dcp.py:167, and
## `canonicalize_expr(..., canonicalize_params=False)`, dqcp2dcp.py:124).
##
## CVXR had no `canonicalize_params` flag, so it expanded the chain into cone
## epigraphs over fresh auxiliary variables:
##     var127 >= t^2 ;  var143 >= 1/(var127+1) ;  var143 >= var158^2 ;  var158 <= x
## Each step is a valid relaxation, but relaxing t^2 UPWARD propagates through
## the decreasing 1/(.) and drives the required bound on x down to nothing:
## var158 = 0 with var143 large satisfies all four for ANY x and ANY t. The
## sublevel set stopped depending on t, so every bisection query reported
## feasible, the search walked t to ~0, and the reported point was whatever the
## final subproblem happened to produce -- a different wrong answer per solver.
##
## Single-atom chains (ceil, floor, a lone power) invert to an AFFINE parameter
## expression, need no epigraph, and were always correct -- which is why the
## rest of the DQCP suite stayed green. The control case below pins that.

## -- minimal nested chain ---------------------------------------------------
## sqrt(inv_pos(x)) = x^(-1/2) is decreasing, so over 0.1 <= x <= 4 the minimum
## is at x = 4, value sqrt(1/4) = 0.5. Measured pre-fix: 0.76466326 at x = 1.71.
test_that("DQCP nested inverse chain: sqrt(inv_pos(x)) minimizes correctly", {
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(sqrt(inv_pos(x))), list(x <= 4, x >= 0.1))
  expect_true(is_dqcp(prob))

  val <- psolve(prob, qcp = TRUE, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(val), 0.5, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 4, tolerance = 1e-4)
})

## -- the hypersonic-shape example (cvxr_docs examples/dqcp/hypersonic-shape) --
## minimize sqrt(inv_pos(x^2) - 1) s.t. a/x - (1-b)*sqrt(1-x^2) <= 0.
## The objective is decreasing in x, so the optimum sits at the largest feasible
## x. Solving a/x = (1-b)*sqrt(1-x^2) numerically gives x* = 0.98952382 and
## objective 0.14589815; CVXPY 1.9.2 returns 0.14589882 on CLARABEL/SCS/ECOS.
## Measured pre-fix in CVXR: 1.5187 (CLARABEL), 1.7711 (SCS), 1.5158 (ECOS),
## 1.5276 (MOSEK) -- all reported `optimal`.
test_that("DQCP nested inverse chain: hypersonic wedge matches CVXPY", {
  a <- 0.05; b <- 0.65
  x <- Variable(pos = TRUE)
  obj <- sqrt(inv_pos(power(x, 2)) - 1)
  con <- list(a * inv_pos(x) - (1 - b) * sqrt(1 - power(x, 2)) <= 0)
  prob <- Problem(Minimize(obj), con)
  expect_true(is_dqcp(prob))

  val <- psolve(prob, qcp = TRUE, solver = "CLARABEL")
  ## optimal_inaccurate is the CORRECT status here, not a weakened expectation:
  ## the bisection's last subproblems sit on the boundary and CLARABEL certifies
  ## them only almost-exactly. CVXPY 1.9.2 reports the same on this problem --
  ## measured: obj=0.14589882, x=0.98952373, status=optimal_inaccurate. Pinning
  ## "optimal" would assert something upstream does not achieve either.
  expect_equal(status(prob), "optimal_inaccurate")
  expect_equal(as.numeric(val), 0.14589815, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 0.98952382, tolerance = 1e-4)

  ## The returned point must satisfy the original constraint.
  xv <- as.numeric(value(x))
  expect_lte(a / xv - (1 - b) * sqrt(1 - xv^2), 1e-6)
})

## The answer must not depend on which solver ran the bisection subproblems.
## Pre-fix these disagreed by more than 15% (1.5187 / 1.7711 / 1.5158) and each
## was ~10x from the truth.
##
## Every solver is held to the SAME tolerance, and that tolerance is set by what
## CVXPY itself achieves on this deliberately degenerate problem -- not by what
## makes CVXR pass. Measured against truth 0.14589815:
##      CLARABEL  CVXR 0.14589882  CVXPY 0.14589882   (identical)
##      ECOS      CVXR 0.14589717  CVXPY 0.14589717   (identical)
##      SCS       CVXR 0.14562171  CVXPY 0.14562028   (1.9e-03 from truth BOTH)
## So 1e-3 would fail upstream too; 3e-3 is the tightest bar CVXPY clears with
## every solver, and it still catches the pre-fix answers, which were 1.5-1.8
## (an order of magnitude out, not a tolerance question).
test_that("DQCP nested inverse chain: answer is solver-independent", {
  a <- 0.05; b <- 0.65
  truth <- 0.14589815
  vals <- list()
  for (s in c("CLARABEL", "SCS", "ECOS")) {
    if (!(s %in% installed_solvers())) next
    x <- Variable(pos = TRUE)
    prob <- Problem(Minimize(sqrt(inv_pos(power(x, 2)) - 1)),
                    list(a * inv_pos(x) - (1 - b) * sqrt(1 - power(x, 2)) <= 0))
    vals[[s]] <- tryCatch(
      as.numeric(psolve(prob, qcp = TRUE, solver = s)),
      error = function(e) NA_real_)
  }
  skip_if(length(vals) < 2, "need at least two conic solvers")
  for (s in names(vals)) {
    expect_false(is.na(vals[[s]]), label = paste(s, "solved"))
    expect_equal(vals[[s]], truth, tolerance = 3e-3,
                 label = paste("objective from", s))
  }
})

## -- control: single-atom chains must be unchanged ---------------------------
## These invert to an affine parameter expression and were always correct. They
## guard the fix against over-reach: canonicalize_params defaults to TRUE, so
## nothing outside Dqcp2Dcp may change.
test_that("DQCP single-atom inverse chains still correct", {
  y <- Variable()
  p1 <- Problem(Minimize(ceiling(y)), list(y >= 2.4))
  expect_equal(as.numeric(psolve(p1, qcp = TRUE, solver = "CLARABEL")), 3,
               tolerance = 1e-6)

  z <- Variable()
  p2 <- Problem(Maximize(floor(z)), list(z <= 5.7))
  expect_equal(as.numeric(psolve(p2, qcp = TRUE, solver = "CLARABEL")), 5,
               tolerance = 1e-6)
})

## -- control: the DPP path is untouched --------------------------------------
## canonicalize_params = TRUE is what preserves parameterized constant subtrees
## for parametrized compilation. A DPP problem must still be DPP, still solve,
## and still track a changed parameter value.
test_that("parameterized (DPP) problems unaffected by canonicalize_params", {
  p <- Parameter(nonneg = TRUE)
  x <- Variable(2)
  A <- matrix(c(1, 0, 0, 1), 2, 2)
  bvec <- c(1, 2)
  prob <- Problem(Minimize(sum_squares(A %*% x - bvec) + p * p_norm(x, 1)))
  expect_true(is_dpp(prob))

  value(p) <- 0
  v0 <- psolve(prob, solver = "CLARABEL")
  x0 <- as.numeric(value(x))
  expect_equal(x0, bvec, tolerance = 1e-5)

  value(p) <- 10
  v10 <- psolve(prob, solver = "CLARABEL")
  x10 <- as.numeric(value(x))
  ## A large L1 weight must shrink the solution toward zero.
  expect_true(sum(abs(x10)) < sum(abs(x0)))
  expect_true(v10 > v0)
})
