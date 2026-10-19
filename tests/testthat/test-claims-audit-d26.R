## Claims audit decision 26 (2026-10-11, CVXR 1.9.2.9162): the fourteen
## "divergence-unrecorded" rows. Ported: C0919 (no EvalParams bake for a DPP
## parametric quadratic objective), C0428(a) (is_dpp rejects an unsupported
## context), C0750 (CVXOPT accepts `max_iters`). Recorded under ADR D_PORT.14:
## C0155/C0199, C0511, C0809 are pinned below. The pins describe behavior CVXR
## already had and cannot be seen failing on 9161; the three port blocks were.

d26_chain_classes <- function(prob) {
  vapply(prob@.cache$compile_chain@reductions, function(r) class(r)[[1L]], "")
}

## @cvxpy cvxpy/reductions/solvers/solving_chain.py:202-207
test_that("C0919: a DPP parametric quadratic objective keeps its parameter program (CLARABEL)", {
  m <- 6
  y <- Variable(m); q <- Parameter(m); value(q) <- sin(seq_len(m))
  prob <- Problem(Minimize(sum_squares(y - q)), list(y >= 0))
  expect_true(is_dpp(prob))
  v1 <- psolve(prob, solver = "CLARABEL")
  expect_false(any(grepl("EvalParams", d26_chain_classes(prob))))
  expect_false(is.null(prob@.cache$param_prog))
  expect_equal(v1, sum(pmin(value(q), 0)^2), tolerance = 1e-6)
  for (j in 2:3) {
    value(q) <- sin(seq_len(m)) * j
    expect_equal(psolve(prob, solver = "CLARABEL"), sum(pmin(value(q), 0)^2), tolerance = 1e-6)
    expect_equal(as.numeric(value(y)), as.numeric(pmax(value(q), 0)), tolerance = 1e-5)
  }
  ## The compiled program was reused, not rebuilt.
  expect_false(is.null(prob@.cache$param_prog))
})

## @cvxpy cvxpy/reductions/solvers/solving_chain.py:202-207
test_that("C0919: the same on the QP path (OSQP), with a parametric linear term too", {
  m <- 6
  y <- Variable(m); q <- Parameter(m); p <- Parameter(); value(q) <- sin(seq_len(m)); value(p) <- 1
  prob <- Problem(Minimize(sum_squares(y - q) + p * sum(y)), list(y >= 0, y <= 3))
  expect_true(is_dpp(prob))
  invisible(psolve(prob, solver = "OSQP"))
  expect_false(any(grepl("EvalParams", d26_chain_classes(prob))))
  expect_false(is.null(prob@.cache$param_prog))
  ## Closed form: y_i = clamp(q_i - p/2, 0, 3).
  closed <- function(qv, pv) {
    yv <- as.numeric(pmin(pmax(qv - pv / 2, 0), 3))
    list(y = yv, obj = sum((yv - qv)^2) + pv * sum(yv))
  }
  for (j in 1:3) {
    value(q) <- sin(seq_len(m)) * j; value(p) <- j
    cf <- closed(value(q), j)
    expect_equal(psolve(prob, solver = "OSQP"), cf$obj, tolerance = 1e-4)
    expect_equal(as.numeric(value(y)), cf$y, tolerance = 1e-3)
  }
})

## @cvxpy cvxpy/reductions/solvers/solving_chain.py:202-207 cvxpy/reductions/eliminate_pwl/canonicalizers/cumsum_canon.py:45-55
test_that("C0919: cumsum of a parameter inside a quadratic objective stays on the tensor path", {
  m <- 6
  y <- Variable(m); q <- Parameter(m); value(q) <- sin(seq_len(m))
  prob <- Problem(Minimize(sum_squares(y - cumsum(q))), list(y >= 0))
  expect_true(is_dpp(prob))
  v1 <- psolve(prob, solver = "CLARABEL")
  expect_false(any(grepl("EvalParams", d26_chain_classes(prob))))
  expect_false(is.null(prob@.cache$param_prog))
  expect_equal(v1, sum(pmin(cumsum(value(q)), 0)^2), tolerance = 1e-6)
  value(q) <- sin(seq_len(m)) * 3
  expect_equal(psolve(prob, solver = "CLARABEL"), sum(pmin(cumsum(value(q)), 0)^2), tolerance = 1e-6)
})

## @cvxpy cvxpy/utilities/coeff_extractor.py:244-262 test_coeff_extractor.py::test_issue_2402_vector
test_that("C0919: a parameter scaling an elementwise quadratic (diagonal path) keeps the tensor path", {
  ## minimize sum(pv_i w_i^2) + sum(cc_i w_i): w_i = -cc_i / (2 pv_i),
  ## objective -sum(cc_i^2 / (4 pv_i)). The linear term is what exposes a
  ## dropped or misplaced coefficient (test-quad-coeff-scaling.R).
  cc <- c(1, 2, 3, 4, 5)
  closed <- function(pvv) list(w = -cc / (2 * pvv), obj = -sum(cc^2 / (4 * pvv)))
  for (sv in c("CLARABEL", "OSQP")) {
    w <- Variable(5); pv <- Parameter(5, nonneg = TRUE); value(pv) <- 1:5
    prob <- Problem(Minimize(sum_entries(pv * square(w)) + t(cc) %*% w))
    expect_true(is_dpp(prob))
    cf <- closed(1:5)
    expect_equal(psolve(prob, solver = sv), cf$obj, tolerance = 1e-5)
    expect_equal(as.numeric(value(w)), cf$w, tolerance = 1e-4)
    expect_false(any(grepl("EvalParams", d26_chain_classes(prob))))
    expect_false(is.null(prob@.cache$param_prog))
    value(pv) <- 2 * (1:5)
    cf <- closed(2 * (1:5))
    expect_equal(psolve(prob, solver = sv), cf$obj, tolerance = 1e-5)
    expect_equal(as.numeric(value(w)), cf$w, tolerance = 1e-4)
  }
  ## The issue 2402 shape with a fixed w (CVXPY 1.9.3: 1700 for the pv * square(w) term
  ## alone, -13863.3412934 for the full objective with this Sigma).
  w <- Variable(5); p <- Parameter(value = 2, nonneg = TRUE)
  pv <- Parameter(5, value = 0:4, nonneg = TRUE)
  Sigma <- crossprod(matrix(sin(1:25), 5, 5)) / 100 + diag(5) * 1e-3
  r <- c(-0.48, 0.11, 0.09, -0.39, 0.03)
  fixed <- c(10, 11, 12, 13, 14)
  q2 <- Problem(Minimize(sum_entries(pv * square(w))), list(w == fixed))
  expect_equal(psolve(q2, solver = "CLARABEL"), 1700, tolerance = 1e-6)
  q4 <- Problem(Maximize(t(r) %*% w - p * quad_form(w, Sigma) -
                           sum_entries((pv * c(5, 6, 7, 8, 9)) * square(w))),
                list(w == fixed))
  expected <- sum(r * fixed) - 2 * drop(t(fixed) %*% Sigma %*% fixed) -
    sum((0:4) * c(5, 6, 7, 8, 9) * fixed^2)
  expect_equal(psolve(q4, solver = "CLARABEL"), expected, tolerance = 1e-6)
  expect_equal(expected, -13863.3412934, tolerance = 1e-9)
  expect_false(is.null(q4@.cache$param_prog))
})

## @cvxpy cvxpy/utilities/coeff_extractor.py:240-243 cvxpy/utilities/coeff_extractor.py:300-363
test_that("C0919: a parameter scaling an axis sum_squares (block path) keeps the tensor path", {
  ## minimize sum_j pv_j ||X[, j]||^2 + sum(C * X): X_ij = -C_ij / (2 pv_j).
  C <- matrix(1:6, 3, 2)
  closed <- function(pvv) {
    X <- -C / (2 * matrix(pvv, 3, 2, byrow = TRUE))
    list(X = X, obj = sum(rep(pvv, each = 3) * X^2) + sum(C * X))
  }
  ## sum_squares(X, axis = 2) is a 1 x 2 row (one entry per column), so the
  ## parameter must be a 1 x 2 row too, or `*` would broadcast to 2 x 2.
  X <- Variable(c(3, 2)); pv <- Parameter(c(1, 2), nonneg = TRUE); value(pv) <- matrix(c(1, 2), 1)
  prob <- Problem(Minimize(sum_entries(pv * sum_squares(X, axis = 2)) + sum_entries(C * X)))
  expect_true(is_dpp(prob))
  cf <- closed(c(1, 2))
  expect_equal(psolve(prob, solver = "CLARABEL"), cf$obj, tolerance = 1e-5)
  expect_equal(as.matrix(value(X)), cf$X, tolerance = 1e-4, ignore_attr = TRUE)
  expect_false(any(grepl("EvalParams", d26_chain_classes(prob))))
  expect_false(is.null(prob@.cache$param_prog))
  value(pv) <- matrix(c(3, 0.5), 1)
  cf <- closed(c(3, 0.5))
  expect_equal(psolve(prob, solver = "CLARABEL"), cf$obj, tolerance = 1e-5)
  expect_equal(as.matrix(value(X)), cf$X, tolerance = 1e-4, ignore_attr = TRUE)
})

## @cvxpy cvxpy/atoms/atom.py:300-308
test_that("C0428: is_dpp rejects an unsupported context, as CVXPY's ValueError", {
  x <- Variable()
  expect_error(is_dpp(x + 1, context = "bogus"), "Unsupported context")
  expect_true(is_dpp(x + 1, context = "DCP"))
  expect_true(is_dpp(x + 1))
})

## @cvxpy cvxpy/reductions/solvers/conic_solvers/cvxopt_conif.py:278-279
test_that("C0750: CVXOPT accepts max_iters as CVXPY's alias of maxiters", {
  require_solver("CVXOPT")
  make <- function() {
    w <- Variable(3)
    Problem(Minimize(sum_squares(w - c(1, 2, 3))), list(w >= 0, sum(w) <= 2))
  }
  ## One iteration cannot converge: the solver reports failure, as CVXPY does
  ## for `max_iters=1` (renamed to cvxopt's `maxiters`).
  expect_error(psolve(make(), solver = "CVXOPT", max_iters = 1L), class = "SolverError")
  expect_equal(psolve(make(), solver = "CVXOPT", max_iters = 50L), 5.5, tolerance = 1e-5)
})

## RECORDED DIVERGENCE (ADR D_PORT.14, claims audit C0155/C0199): upstream
## condition_number/gen_lambda_max `_grad` raise NotImplementedError; only
## dist_ratio returns None. CVXR returns NULL for all three.
## @cvxpy cvxpy/atoms/condition_number.py:45-56 cvxpy/atoms/gen_lambda_max.py:47-58 cvxpy/atoms/dist_ratio.py:83-84
test_that("C0155/C0199 (pin): DQCP atoms report an undefined gradient as NULL, like dist_ratio", {
  A <- Variable(c(2, 2), PSD = TRUE); value(A) <- diag(c(2, 1))
  B <- Variable(c(2, 2), PSD = TRUE); value(B) <- diag(2)
  g1 <- CVXR:::grad(condition_number(A))
  expect_length(g1, 1L); expect_true(all(vapply(g1, is.null, logical(1))))
  g2 <- CVXR:::grad(gen_lambda_max(A, B))
  expect_length(g2, 2L); expect_true(all(vapply(g2, is.null, logical(1))))
  a <- Variable(2); value(a) <- c(1, 0)
  g3 <- CVXR:::grad(dist_ratio(a, c(0, 0), c(3, 3)))
  expect_length(g3, 1L); expect_true(all(vapply(g3, is.null, logical(1))))
})

## @cvxpy cvxpy/problems/objective.py:61-95
test_that("C0511 (pin): an objective plus or minus the scalar 0 is the objective, on either side", {
  x <- Variable()
  expect_s3_class(Minimize(x) + 0, "CVXR::Minimize")
  expect_s3_class(0 + Minimize(x), "CVXR::Minimize")
  expect_s3_class(Minimize(x) - 0, "CVXR::Minimize")
  expect_s3_class(0 - Minimize(x), "CVXR::Maximize")
  expect_error(Minimize(x) + 1, "only be added")
  expect_error(Minimize(x) - 1, "only be subtracted")
})

## RECORDED DIVERGENCE (ADR D_PORT.14, claims audit C0809): upstream skips a
## row with no coefficients (scip_conif.py:421-432); CVXR keeps it as 0 == b / 0 <= b.
## @cvxpy cvxpy/reductions/solvers/conic_solvers/scip_conif.py:421-432
test_that("C0809 (pin): SCIP honors a constraint row with no variable coefficients", {
  require_solver("SCIP")
  t <- Variable()
  infeasible <- list(
    Problem(Minimize(t), list(t >= 0, Constant(1) <= 0)),
    Problem(Minimize(t), list(t >= 0, Constant(1) == 0)),
    Problem(Minimize(t), list(t >= 0, 0 * t == 1))
  )
  for (pr in infeasible) {
    expect_equal(psolve(pr, solver = "SCIP"), Inf)
    expect_equal(status(pr), "infeasible")
  }
  feasible <- Problem(Minimize(t), list(t >= 0, Constant(-1) <= 0))
  expect_equal(psolve(feasible, solver = "SCIP"), 0, tolerance = 1e-6)
  expect_equal(status(feasible), "optimal")
})
