## DGP CVXPY parity tests — closing critical gaps vs cvxpy/tests/test_dgp2dcp.py
## Tests ported from CVXPY commit 3b964472b (Release 1.8.1)
## Expected values verified against CVXPY using `uv run python`

## ══════════════════════════════════════════════════════════════════
## CRITICAL tests (10)
## ══════════════════════════════════════════════════════════════════

## ── test_unconstrained_monomial ──────────────────────────────────
## CVXPY: Minimize(x * y) with no constraints → unbounded (value ~ 0)
## In log-space, the DCP problem is unbounded below (→ -Inf),
## so exp(-Inf) = 0. Status: "unbounded".

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_unconstrained_monomial
test_that("CVXPY parity: unconstrained monomial minimize", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x * y), list())
  val <- psolve(prob, gp = TRUE, solver = "CLARABEL")
  ## Unbounded: x*y can be driven to 0+
  expect_equal(val, 0.0, tolerance = 1e-4)
  expect_equal(status(prob), "unbounded")
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_unconstrained_monomial
test_that("CVXPY parity: unconstrained monomial maximize", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Maximize(x * y), list())
  val <- psolve(prob, gp = TRUE, solver = "CLARABEL")
  ## Unbounded: x*y can be driven to +Inf
  expect_equal(val, Inf)
  expect_equal(status(prob), "unbounded")
})

## ── test_basic_equality_constraint ──────────────────────────────
## CVXPY: Minimize(x), x == 1.0 → value 1.0, x = 1.0

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_basic_equality_constraint
test_that("CVXPY parity: basic equality constraint", {
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x), list(x == 1.0))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
})

## ── test_add_canon ──────────────────────────────────────────────
## DGP add canonicalizer: X + Y in log-space becomes log(exp(X) + exp(Y))
## This tests the internal canonicalizer directly on constants.

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_add_canon
test_that("CVXPY parity: add_canon matrix", {
  X_val <- matrix(c(1, 4, 2, 5, 3, 6), nrow = 2, ncol = 3)
  Y_val <- matrix(c(2, 5, 3, 6, 4, 7), nrow = 2, ncol = 3)
  X <- Constant(X_val)
  Y <- Constant(Y_val)
  Z <- X + Y

  result <- CVXR:::.dgp_add_canon(Z, Z@args)
  canon_matrix <- result[[1L]]
  constraints <- result[[2L]]

  expect_equal(length(constraints), 0L)
  expect_equal(canon_matrix@shape, Z@shape)
  expected <- log(exp(X_val) + exp(Y_val))
  expect_equal(as.numeric(value(canon_matrix)), as.numeric(expected),
               tolerance = 1e-5)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_add_canon
test_that("CVXPY parity: add_canon with scalar promotion", {
  X_val <- matrix(c(1, 4, 2, 5, 3, 6), nrow = 2, ncol = 3)
  X <- Constant(X_val)
  y <- Constant(2.0)
  Z <- X + y

  result <- CVXR:::.dgp_add_canon(Z, Z@args)
  canon_matrix <- result[[1L]]
  constraints <- result[[2L]]

  expect_equal(length(constraints), 0L)
  expect_equal(canon_matrix@shape, Z@shape)
  expected <- log(exp(X_val) + exp(2.0))
  expect_equal(as.numeric(value(canon_matrix)), as.numeric(expected),
               tolerance = 1e-5)
})

## ── test_matmul_canon ───────────────────────────────────────────
## DGP matmul canonicalizer: (X @ Y)[i,j] = log(sum_k exp(X[i,k] + Y[k,j]))
## For constants X (2x3) and Y (3x1):
## Entry [1,1] = log(exp(1+1) + exp(2+2) + exp(3+3)) = log(exp(2)+exp(4)+exp(6))
## Entry [2,1] = log(exp(4+1) + exp(5+2) + exp(6+3)) = log(exp(5)+exp(7)+exp(9))

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_matmul_canon
test_that("CVXPY parity: matmul_canon", {
  X_val <- matrix(c(1, 4, 2, 5, 3, 6), nrow = 2, ncol = 3)
  Y_val <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
  X <- Constant(X_val)
  Y <- Constant(Y_val)
  Z <- X %*% Y

  result <- CVXR:::.dgp_mulexpression_canon(Z, Z@args)
  canon_matrix <- result[[1L]]
  constraints <- result[[2L]]

  expect_equal(length(constraints), 0L)
  expect_equal(canon_matrix@shape, c(2L, 1L))
  first_entry <- log(exp(2) + exp(4) + exp(6))
  second_entry <- log(exp(5) + exp(7) + exp(9))
  canon_val <- value(canon_matrix)
  expect_equal(as.numeric(canon_val[1, 1]), first_entry, tolerance = 1e-5)
  expect_equal(as.numeric(canon_val[2, 1]), second_entry, tolerance = 1e-5)
})

## ── test_trace_canon ────────────────────────────────────────────
## DGP trace canonicalizer: trace(X) → add_canon(diag(X))
## trace([[1, 5], [9, 14]]) in log-space = log(exp(1) + exp(14))

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_trace_canon
test_that("CVXPY parity: trace_canon", {
  X_val <- matrix(c(1, 9, 5, 14), nrow = 2, ncol = 2)
  X <- Constant(X_val)
  Y <- Trace(X)

  result <- CVXR:::.dgp_trace_canon(Y, Y@args)
  canon <- result[[1L]]
  constraints <- result[[2L]]

  expect_equal(length(constraints), 0L)
  expect_true(is_scalar(canon))
  expected <- log(exp(1) + exp(14))
  expect_equal(as.numeric(value(canon)), expected, tolerance = 1e-5)
})

## ── test_paper_example_sum_largest ──────────────────────────────
## CVXPY skips this test: "Enable test once sum_largest is implemented."
## We test whether sum_largest works in DGP context.
## Minimizing sum of 2 largest of [3*x0^0.5*x1^0.5, x0*x1+0.5*x1*x3^3, x2]
## subject to x0*x1 >= 16
## CVXPY expected: 28.0 (at x0=4, x1=4, x3~0, x2~0)

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_paper_example_sum_largest
test_that("CVXPY parity: paper example sum_largest (composition with posynomials)", {
  ## sum_largest with DGP is a feature that CVXPY itself skips.
  ## We test the composition variant if sum_largest works in DGP.
  skip_if_not(
    tryCatch({
      x_test <- Variable(4, pos = TRUE)
      e <- sum_largest(x_test, 2L)
      is_log_log_convex(e)
    }, error = function(e) FALSE),
    "sum_largest not supported in DGP context"
  )

  x <- Variable(4, pos = TRUE)
  obj <- Minimize(sum_largest(
    hstack(3 * power(x[1], 0.5) * power(x[2], 0.5),
           x[1] * x[2] + 0.5 * x[2] * power(x[4], 3),
           x[3]), 2L))
  constr <- list(x[1] * x[2] >= 16)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## Expected: 28.0 (3*sqrt(4)*sqrt(4) + 4*4 + 0.5*4*0 = 12+16 = 28)
  expect_equal(val, 28.0, tolerance = 0.5)
})

## ── test_paper_example_one_minus_pos ────────────────────────────
## CVXPY: Minimize(x*y), (y * one_minus_pos(x/y))^2 >= 1, x >= y/3
## This is a smoke test (CVXPY doesn't check exact values either).

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_paper_example_one_minus_pos
test_that("CVXPY parity: paper example one_minus_pos", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  obj <- Minimize(x * y)
  constr <- list(power(y * one_minus_pos(x / y), 2) >= 1, x >= y / 3)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## Smoke test: problem should solve successfully with finite value
  expect_true(is.finite(val))
  expect_true(val > 0)
  expect_equal(status(prob), "optimal")
})

## ── test_paper_example_eye_minus_inv ────────────────────────────
## CVXPY: Minimize(trace(eye_minus_inv(X))),
##   geo_mean(diag(X)) == 0.1, geo_mean(hstack(X[1,2], X[2,1])) == 0.1
## Expected: X = 0.1 * ones(2,2), value = 2.25
## (I - 0.1*ones)^{-1} = [[10/9, 1/9],[1/9, 10/9]]/... → trace ~ 2.25

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_paper_example_eye_minus_inv
test_that("CVXPY parity: paper example eye_minus_inv", {
  X <- Variable(c(2, 2), pos = TRUE)
  obj <- Minimize(matrix_trace(eye_minus_inv(X)))
  constr <- list(
    geo_mean(DiagMat(X)) == 0.1,
    geo_mean(hstack(X[1, 2], X[2, 1])) == 0.1
  )
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## Expected: trace((I - 0.1*ones(2,2))^{-1}) = 2.25
  expect_equal(val, 2.25, tolerance = 0.05)
  ## All entries of X should be ~0.1
  X_val <- as.numeric(value(X))
  expect_equal(X_val, rep(0.1, 4), tolerance = 0.01)
})

## ── test_paper_example_exp_log ──────────────────────────────────
## CVXPY: Minimize(x*y), exp(y/x) <= log(y)
## This is a smoke test in CVXPY.

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_paper_example_exp_log
test_that("CVXPY parity: paper example exp_log", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  obj <- Minimize(x * y)
  constr <- list(exp(y / x) <= log(y))
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## Smoke test: problem should solve successfully with finite value
  expect_true(is.finite(val))
  expect_true(val > 0)
  expect_equal(status(prob), "optimal")
})

## ── test_rank_one_nmf ───────────────────────────────────────────
## Rank-one non-negative matrix factorization via GP.
## CVXPY: minimize sum of element-wise max(X/xy, xy/X) for rank-1 approx xy^T
## with known entries and constraint x0*x1*x2 == 1.

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_rank_one_nmf
test_that("CVXPY parity: rank one NMF", {
  X <- Variable(c(3, 3), pos = TRUE)
  x <- Variable(3, pos = TRUE)
  y <- Variable(3, pos = TRUE)

  ## Outer product: xy[i,j] = x[i] * y[j]
  xy <- vstack(x[1] * y[1], x[1] * y[2], x[1] * y[3],
               x[2] * y[1], x[2] * y[2], x[2] * y[3],
               x[3] * y[1], x[3] * y[2], x[3] * y[3])
  xy <- reshape_expr(xy, c(3L, 3L))

  R <- max_elemwise(
    multiply(X, power(xy, -1)),
    multiply(power(X, -1), xy)
  )
  objective <- sum(R)
  constraints <- list(
    X[1, 1] == 1.0,
    X[1, 3] == 1.9,
    X[2, 2] == 0.8,
    X[3, 1] == 3.2,
    X[3, 2] == 5.9,
    x[1] * x[2] * x[3] == 1.0
  )
  prob <- Problem(Minimize(objective), constraints)
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## Smoke test: problem should solve with finite optimal value
  expect_true(is.finite(val))
  expect_true(val >= 9)  ## at least 9 (each max >= 1, 9 entries)
  expect_equal(status(prob), "optimal")
})

## ══════════════════════════════════════════════════════════════════
## MEDIUM tests (5)
## ══════════════════════════════════════════════════════════════════

## ── test_solving_non_dgp_problem_raises_error ───────────────────
## CVXPY: solve(gp=True) on non-DGP problem → DGPError with hint about DCP
## In CVXR: psolve(prob, gp = TRUE) → error "not DGP compliant"

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_solving_non_dgp_problem_raises_error
test_that("CVXPY parity: non-DGP with gp=TRUE raises error", {
  ## Minimize(-x) is not DGP (negative coefficient)
  x <- Variable()
  prob <- Problem(Minimize(-1.0 * x), list())
  expect_error(psolve(prob, gp = TRUE), "not DGP")
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_solving_non_dgp_problem_raises_error
test_that("CVXPY parity: non-DGP solves normally as DCP", {
  ## Same problem is DCP: Minimize(-x) → unbounded below
  x <- Variable()
  prob <- Problem(Minimize(-1.0 * x), list())
  val <- psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "unbounded")
  expect_equal(val, -Inf)
})

## ── test_solving_non_dcp_problem_raises_error ───────────────────
## CVXPY: solve() on DGP-only problem → DCPError with hint about DGP
## In CVXR: psolve(prob) → error suggesting gp = TRUE

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_solving_non_dcp_problem_raises_error
test_that("CVXPY parity: DGP-only problem without gp=TRUE suggests DGP", {
  ## x * y is not DCP (product of two variables) but is DGP
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x * y))
  expect_error(psolve(prob), "gp = TRUE")
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_solving_non_dcp_problem_raises_error
test_that("CVXPY parity: DGP-only problem solves with gp=TRUE", {
  ## Same problem solved with gp=TRUE → unbounded (value ~ 0)
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x * y))
  val <- psolve(prob, gp = TRUE, solver = "CLARABEL")
  expect_equal(val, 0.0, tolerance = 1e-4)
  expect_equal(status(prob), "unbounded")
})

## ── test_solving_non_dcp_problems_raises_detailed_error ─────────
## CVXPY: solve() on non-DCP problems → detailed DCPError messages

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_solving_non_dcp_problems_raises_detailed_error
test_that("CVXPY parity: non-DCP objective raises detailed error", {
  ## sum(x) - sum_squares(x) is not DCP (concave - convex = not DCP)
  x <- Variable(3)
  prob <- Problem(Minimize(sum(x) - sum_squares(x)))
  expect_error(psolve(prob), "not DCP")
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_solving_non_dcp_problems_raises_detailed_error
test_that("CVXPY parity: non-DCP constraint raises error", {
  ## x * x <= 5 is not DCP (quadratic in constraint LHS non-convex sense)
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x), list(x * x <= 5))
  expect_error(psolve(prob), "not DCP")
})

## ── test_simpler_eye_minus_inv ──────────────────────────────────
## CVXPY: Minimize(trace(eye_minus_inv(X))), diag(X) == 0.1,
##   hstack(X[0,1], X[1,0]) == 0.1
## Expected: X = 0.1*ones(2,2), value = 2.25

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_simpler_eye_minus_inv
test_that("CVXPY parity: simpler eye_minus_inv", {
  X <- Variable(c(2, 2), pos = TRUE)
  obj <- Minimize(matrix_trace(eye_minus_inv(X)))
  constr <- list(
    DiagMat(X) == matrix(c(0.1, 0.1), ncol = 1),
    hstack(X[1, 2], X[2, 1]) == matrix(c(0.1, 0.1), nrow = 1)
  )
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 2.25, tolerance = 0.05)
  X_val <- as.numeric(value(X))
  expect_equal(X_val, rep(0.1, 4), tolerance = 0.01)
})

## ── test_parameter_name ─────────────────────────────────────────
## CVXPY: Parameter with name "alpha" → after DGP reduction, the
## log-space parameter preserves the name.

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_parameter_name
test_that("CVXPY parity: parameter name preserved in DGP reduction", {
  param <- Parameter(pos = TRUE, name = "alpha")
  value(param) <- 1.0
  prob <- Problem(Minimize(param), list())

  ## Apply Dgp2Dcp reduction manually
  dgp2dcp <- Dgp2Dcp()
  result <- reduction_apply(dgp2dcp, prob)
  dcp_prob <- result[[1L]]

  ## The log-space parameter should preserve the name "alpha"
  dcp_params <- parameters(dcp_prob)
  expect_true(length(dcp_params) >= 1L)
  param_names <- vapply(dcp_params, function(p) p@.name, character(1L))
  expect_true("alpha" %in% param_names)
})

## ══════════════════════════════════════════════════════════════════
## Additional parity tests from CVXPY test_dgp2dcp.py
## ══════════════════════════════════════════════════════════════════

## ── test_prod (DGP curvature checks) ────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_prod
test_that("CVXPY parity: prod curvature in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  posy1 <- x * sqrt(y) + 3.0 * x * sqrt(y)
  posy2 <- x * sqrt(y) + 3.0 * power(x, 2) * sqrt(y)
  ## Product of two posynomials: log-log convex (posynomial)
  expect_true(is_log_log_convex(prod(posy1, posy2)))
  expect_false(is_log_log_concave(prod(posy1, posy2)))
  ## Product of posynomial and its inverse: not DGP
  expect_false(is_dgp(prod(posy1, 1 / posy1)))

  ## Product of two monomials: log-log affine
  m <- x * sqrt(y)
  expect_true(is_log_log_affine(prod(m, m)))
  ## Product of monomial and inverse posynomial: log-log concave
  expect_true(is_log_log_concave(prod(m, 1 / posy1)))
  expect_false(is_log_log_convex(prod(m, 1 / posy1)))
})

## ── test_geo_mean (unbounded) ───────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_geo_mean
test_that("CVXPY parity: geo_mean unbounded", {
  x <- Variable(3, pos = TRUE)
  p <- c(1, 2, 0.5)
  prob <- Problem(Minimize(geo_mean(x, p)), list())
  val <- psolve(prob, gp = TRUE, solver = "CLARABEL")
  ## geo_mean can be driven to 0 → unbounded
  expect_equal(val, 0.0, tolerance = 1e-4)
  expect_equal(status(prob), "unbounded")
})

## ── test_div ────────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_div
test_that("CVXPY parity: division constraint", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x * y), list(y / 3 <= x, y >= 1))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 1.0 / 3.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 1.0 / 3.0, tolerance = 1e-2)
})

## ── test_qp_solver_not_allowed ──────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_qp_solver_not_allowed
test_that("CVXPY parity: QP solver not allowed for DGP", {
  ## CVXPY raises an error when a QP solver (e.g., OSQP) is used with gp=TRUE.
  ## CVXR differs: DGP reduces to a DCP problem which may be solvable by QP solvers
  ## (the log-space problem can be an LP/QP). CVXR does not enforce solver type for DGP.
  skip("CVXR allows QP solvers for DGP (differs from CVXPY)")
})

## ── test_sum_scalar ─────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_sum_scalar
test_that("CVXPY parity: sum of scalar in DGP", {
  w <- Variable(pos = TRUE)
  h <- Variable(pos = TRUE)
  prob <- Problem(Minimize(h), list(w * h >= 10, sum(w) <= 5))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(h)), 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(w)), 5.0, tolerance = 1e-2)
})

## ── test_sum_squares_vector ─────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_sum_squares_vector
test_that("CVXPY parity: sum_squares vector in DGP", {
  w <- Variable(2, pos = TRUE)
  h <- Variable(2, pos = TRUE)
  prob <- Problem(Minimize(sum_squares(h)),
                  list(multiply(w, h) >= 10, sum(w) <= 10))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 8.0, tolerance = 0.1)
  expect_equal(as.numeric(value(h)), c(2, 2), tolerance = 0.1)
  expect_equal(as.numeric(value(w)), c(5, 5), tolerance = 0.1)
})

## ── test_trace (DGP solve) ──────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_trace
test_that("CVXPY parity: trace in DGP solve", {
  w <- Variable(c(1, 1), pos = TRUE)
  h <- Variable(pos = TRUE)
  prob <- Problem(Minimize(h), list(w * h >= 10, matrix_trace(w) <= 5))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(h)), 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(w)), 5.0, tolerance = 1e-2)
})

## ── test_parameter (DGP parameter re-solve) ─────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_parameter
test_that("CVXPY parity: parameter in DGP", {
  param <- Parameter(pos = TRUE)
  x <- Variable(pos = TRUE)

  value(param) <- 1.0
  prob <- Problem(Minimize(x), list(x == param))
  val1 <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val1, 1.0, tolerance = 1e-3)

  value(param) <- 2.0
  val2 <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val2, 2.0, tolerance = 1e-3)
})

## ── test_parameter DGP reduction transforms to log-space ────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_parameter
test_that("CVXPY parity: parameter value in log-space after DGP reduction", {
  param <- Parameter(pos = TRUE)
  value(param) <- 1.0
  prob <- Problem(Minimize(param), list())

  dgp2dcp <- Dgp2Dcp()
  result <- reduction_apply(dgp2dcp, prob)
  dcp_prob <- result[[1L]]

  ## The DCP problem's parameter should have value = log(1) = 0
  dcp_params <- parameters(dcp_prob)
  expect_true(length(dcp_params) >= 1L)
  expect_equal(as.numeric(value(dcp_params[[1L]])), log(1.0), tolerance = 1e-10)
})

## ── test_gmatmul ────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_gmatmul
test_that("CVXPY parity: gmatmul value evaluation", {
  x <- Variable(2, pos = TRUE)
  A <- matrix(c(-5, 1, 2, -3), 2, 2)  ## col-major: [[-5,2],[1,-3]]
  b <- c(3, 2)
  expr <- gmatmul(A, x)
  value(x) <- b
  ## gmatmul(A, x)[i] = prod_j x[j]^A[i,j]
  ## entry 1: 3^(-5) * 2^(2) = 4/243 ~ 0.01646
  ## entry 2: 3^(1) * 2^(-3) = 3/8 = 0.375
  expected <- c(3^(-5) * 2^(2), 3^(1) * 2^(-3))
  expect_equal(as.numeric(value(expr)), expected, tolerance = 1e-5)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_gmatmul
test_that("CVXPY parity: gmatmul solve", {
  x <- Variable(2, pos = TRUE)
  A <- matrix(c(-5, 1, 2, -3), 2, 2)  ## col-major: [[-5,2],[1,-3]]
  b <- c(3, 2)
  prob <- Problem(Minimize(1.0), list(gmatmul(A, x) == b))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## Solution: x = exp(solve(A, log(b)))
  expected_x <- exp(solve(A, log(b)))
  expect_equal(as.numeric(value(x)), as.numeric(expected_x), tolerance = 1e-2)
})

## ── test_xexp ───────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_xexp
test_that("CVXPY parity: xexp value evaluation", {
  x <- Variable(2, pos = TRUE)
  b <- c(1, 0.5)
  expr <- xexp(x)
  value(x) <- b
  ## xexp(x) = x * exp(x)
  expected <- c(1 * exp(1), 0.5 * exp(0.5))
  expect_equal(as.numeric(value(expr)), expected, tolerance = 1e-5)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_xexp
test_that("CVXPY parity: xexp solve", {
  x <- Variable(2, pos = TRUE)
  b <- c(1, 0.5)
  prob <- Problem(Minimize(prod(xexp(x))), list(x >= b))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## At lower bounds: xexp([1, 0.5]) = [e, 0.5*sqrt(e)]
  ## product = e * 0.5 * sqrt(e) = 0.5 * e^1.5
  expected <- 0.5 * exp(1)^1.5
  expect_equal(val, expected, tolerance = 1e-2)
})

## ── test_pnorm ──────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_pnorm
test_that("CVXPY parity: pnorm p=2 in DGP", {
  x <- Variable(pos = TRUE)
  arr <- c(3, 4)
  prob <- Problem(Minimize(p_norm(Constant(arr), 2) * power(x, 2)),
                  list(x >= 1))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 5.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_pnorm
test_that("CVXPY parity: pnorm p=3 in DGP", {
  x <- Variable(pos = TRUE)
  arr <- c(1.5, 3, 2)
  l3_norm <- (sum(arr^3))^(1/3)
  prob <- Problem(Minimize(p_norm(Constant(arr), 3) * power(x, 2)),
                  list(x >= 1))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, l3_norm, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
})

## ── test_maximum (DGP solve) ────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_maximum
test_that("CVXPY parity: maximum in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  obj <- Minimize(max_elemwise(x * sqrt(y), 3 * x * sqrt(y)))
  prob <- Problem(obj, list(x == 1, y == 4))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## max(1*2, 3*1*2) = max(2, 6) = 6
  expect_equal(val, 6.0, tolerance = 0.1)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 4.0, tolerance = 1e-2)
})

## ── test_minimum (DGP solve) ────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_minimum
test_that("CVXPY parity: minimum in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prod1 <- x * sqrt(y)
  prod2 <- 3 * x * sqrt(y)
  posy <- prod1 + prod2
  obj <- Maximize(min_elemwise(prod1, prod2, 1 / posy))
  prob <- Problem(obj, list(x == 1, y == 4))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## min(2, 6, 1/(2+6)) = min(2, 6, 0.125) = 0.125
  expect_equal(val, 1.0 / (2.0 + 6.0), tolerance = 1e-2)
})

## ── test_documentation_prob ─────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_documentation_prob
test_that("CVXPY parity: documentation example", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  prob <- Problem(Maximize(x * y * z),
                  list(4 * x * y * z + 2 * x * z <= 10,
                       x <= 2 * y, y <= 2 * x, z >= 1))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## Smoke test: should solve successfully
  expect_true(is.finite(val))
  expect_true(val > 0)
  expect_equal(status(prob), "optimal")
})

## ── test_solver_error ───────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_solver_error
test_that("CVXPY parity: solver error propagation in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x * y), list())

  dgp2dcp <- Dgp2Dcp()
  result <- reduction_apply(dgp2dcp, prob)
  dcp_prob <- result[[1L]]
  inverse_data <- result[[2L]]

  ## Simulate a solver error solution
  soln <- Solution(SOLVER_ERROR, NA_real_, list(), list(), list())
  dgp_soln <- reduction_invert(dgp2dcp, soln, inverse_data)
  expect_equal(dgp_soln@status, SOLVER_ERROR)
})

## ── test_one_minus_pos (direct) ─────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_one_minus_pos
test_that("CVXPY parity: one_minus_pos direct", {
  x <- Variable(pos = TRUE)
  prob <- Problem(Maximize(x), list(one_minus_pos(x) >= 0.4))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## max x s.t. 1-x >= 0.4 → x <= 0.6 → optimal x = 0.6
  expect_equal(val, 0.6, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 0.6, tolerance = 1e-2)
})

## ── test_pf_matrix_completion ───────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_pf_matrix_completion
test_that("CVXPY parity: PF matrix completion", {
  X <- Variable(c(3, 3), pos = TRUE)
  obj <- Minimize(pf_eigenvalue(X))
  constr <- list(
    X[1, 1] == 1.0, X[1, 3] == 1.9,
    X[2, 2] == 0.8,
    X[3, 1] == 3.2, X[3, 2] == 5.9,
    X[1, 2] * X[2, 1] * X[2, 3] * X[3, 3] == 1.0
  )
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## Smoke test: should produce finite optimal
  expect_true(is.finite(val))
  expect_true(val > 0)
  expect_equal(status(prob), "optimal")
})

## ── Unconstrained monomial: reduction-level test ────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_unconstrained_monomial
test_that("CVXPY parity: unconstrained monomial reduction internals", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prod_expr <- x * y
  prob <- Problem(Minimize(prod_expr), list())

  dgp2dcp <- Dgp2Dcp()
  result <- reduction_apply(dgp2dcp, prob)
  dcp_prob <- result[[1L]]

  ## In log-space, the DCP problem objective should be addition of two variables
  ## (log(x*y) = log(x) + log(y))
  obj_expr <- dcp_prob@objective@args[[1L]]
  ## It should be an AddExpression with 2 args (both Variables)
  expect_true(S7_inherits(obj_expr, AddExpression))
  expect_equal(length(obj_expr@args), 2L)
  expect_true(S7_inherits(obj_expr@args[[1L]], Variable))
  expect_true(S7_inherits(obj_expr@args[[2L]], Variable))
})

## ── Basic equality: reduction-level test ────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_basic_equality_constraint
test_that("CVXPY parity: basic equality reduction internals", {
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x), list(x == 1.0))

  dgp2dcp <- Dgp2Dcp()
  result <- reduction_apply(dgp2dcp, prob)
  dcp_prob <- result[[1L]]

  ## In log-space, the objective should be a single Variable
  obj_expr <- dcp_prob@objective@args[[1L]]
  expect_true(S7_inherits(obj_expr, Variable))
})

## ══════════════════════════════════════════════════════════════════
## Additional CVXPY parity tests: test_max, test_min, test_sum_largest
## ══════════════════════════════════════════════════════════════════

## ── test_max (DGP: max of hstack of posynomials) ─────────────────
## CVXPY: Minimize(max(hstack([prod1, prod2]))), x == 1, y == 4
## prod1 = x * y^0.5 = 1 * 2 = 2, prod2 = 3*x*y^0.5 = 6
## max(2, 6) = 6

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_max
test_that("CVXPY parity: test_max — max of hstack of posynomials in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prod1 <- x * power(y, 0.5)
  prod2 <- 3.0 * x * power(y, 0.5)
  obj <- Minimize(max_entries(hstack(prod1, prod2)))
  constr <- list(x == 1.0, y == 4.0)
  dgp <- Problem(obj, constr)
  val <- psolve(dgp, solver = "CLARABEL", gp = TRUE)
  expect_equal(val, 6.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 4.0, tolerance = 1e-2)
})

## ── test_min (DGP: max of min of monomial, posynomial inverse) ───
## CVXPY: Maximize(min(hstack([prod1, prod2, 1/posy]))), x == 1, y == 4
## prod1 = 2, prod2 = 6, posy = prod1 + prod2 = 8, 1/posy = 0.125
## min(2, 6, 0.125) = 0.125

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_min
test_that("CVXPY parity: test_min — min of hstack in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prod1 <- x * power(y, 0.5)
  prod2 <- 3.0 * x * power(y, 0.5)
  posy <- prod1 + prod2
  obj <- Maximize(min_entries(hstack(prod1, prod2, 1 / posy)))
  constr <- list(x == 1.0, y == 4.0)
  dgp <- Problem(obj, constr)
  val <- psolve(dgp, solver = "CLARABEL", gp = TRUE)
  ## min(2, 6, 1/8) = 0.125
  expect_equal(val, 1.0 / (2.0 + 6.0), tolerance = 1e-2)
})

## ── test_sum_largest (DGP: skipped in CVXPY) ─────────────────────
## CVXPY: self.skipTest("Enable test once sum_largest is implemented.")

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_sum_largest
test_that("CVXPY parity: test_sum_largest — skipped (not implemented in DGP)", {
  skip("CVXPY skips: sum_largest not implemented in DGP")
})

## ══════════════════════════════════════════════════════════════════
## DGP curvature tests from CVXPY test_dgp.py::TestDgp
## Ported from CVXPY commit 3b964472b (Release 1.8.1)
## ══════════════════════════════════════════════════════════════════

## ── test_product ──────────────────────────────────────────────────
## CVXPY: product of positive variables is DGP, log-log convex and concave.
## Multiplying by negative constant breaks DGP.

## @cvxpy test_dgp.py::TestDgp::test_product
test_that("CVXPY parity: product of positive variables — DGP curvature", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## x * y is DGP, log-log affine

  prod_xy <- x * y
  expect_true(is_dgp(prod_xy))
  expect_true(is_log_log_convex(prod_xy))
  expect_true(is_log_log_concave(prod_xy))

  ## (x*y)^2 is still DGP and log-log affine (monomial)
  prod2 <- prod_xy * prod_xy
  expect_true(is_dgp(prod2))
  expect_true(is_log_log_convex(prod2))
  expect_true(is_log_log_concave(prod2))

  ## 5.0 * (x*y)^2 is still DGP (positive constant times monomial)
  prod3 <- 5.0 * prod2
  expect_true(is_dgp(prod3))
  expect_true(is_log_log_convex(prod3))
  expect_true(is_log_log_concave(prod3))

  ## -5.0 * (x*y)^2 breaks DGP (negative constant)
  prod4 <- -5.0 * prod2
  expect_false(is_dgp(prod4))
  expect_false(is_log_log_convex(prod4))
  expect_false(is_log_log_concave(prod4))
})

## ── test_product_with_unconstrained_variables_is_not_dgp ──────────

## @cvxpy test_dgp.py::TestDgp::test_product_with_unconstrained_variables_is_not_dgp
test_that("CVXPY parity: product of unconstrained variables is not DGP", {
  x <- Variable()
  y <- Variable()

  ## Unconstrained variables: x * y is not DGP
  prod_xy <- x * y
  expect_false(is_dgp(prod_xy))
  expect_false(is_log_log_convex(prod_xy))
  expect_false(is_log_log_concave(prod_xy))

  ## One positive, one unconstrained: still not DGP
  z <- Variable(pos = TRUE)
  prod_xz <- x * z
  expect_false(is_dgp(prod_xz))
  expect_false(is_log_log_convex(prod_xz))
  expect_false(is_log_log_concave(prod_xz))
})

## ── test_division ─────────────────────────────────────────────────

## @cvxpy test_dgp.py::TestDgp::test_division
test_that("CVXPY parity: division curvature in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## x / y is log-log affine (monomial)
  div1 <- x / y
  expect_true(is_log_log_affine(div1))

  ## posynomial / monomial is log-log convex
  posynomial <- 5.0 * x * y + 1.2 * y * y
  div2 <- x / y
  expect_true(is_log_log_affine(div2))

  ## posynomial / monomial = log-log convex (posynomial numerator)
  div3 <- posynomial / (3.0 * x * power(y, -0.1))
  expect_true(is_log_log_convex(div3))
  expect_false(is_log_log_concave(div3))
  expect_true(is_dgp(div3))

  ## posynomial / posynomial is NOT DGP
  div4 <- posynomial / (3.0 * x + y)
  expect_false(is_log_log_convex(div4))
  expect_false(is_log_log_concave(div4))
  expect_false(is_dgp(div4))
})

## ── test_add ──────────────────────────────────────────────────────

## @cvxpy test_dgp.py::TestDgp::test_add
test_that("CVXPY parity: addition curvature in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## x + y is DGP, log-log convex (posynomial), NOT log-log concave
  expr <- x + y
  expect_true(is_dgp(expr))
  expect_true(is_log_log_convex(expr))
  expect_false(is_log_log_concave(expr))

  ## More complex posynomial
  posynomial <- 5.0 * x * y + 1.2 * y * y
  expect_true(is_dgp(posynomial))
  expect_true(is_log_log_convex(posynomial))
})

## ── test_add_with_unconstrained_variables_is_not_dgp ──────────────

## @cvxpy test_dgp.py::TestDgp::test_add_with_unconstrained_variables_is_not_dgp
test_that("CVXPY parity: addition with unconstrained variables is not DGP", {
  x <- Variable()
  y <- Variable(pos = TRUE)

  ## x + y where x is unconstrained: not DGP
  expr <- x + y
  expect_false(is_dgp(expr))
  expect_false(is_log_log_convex(expr))
  expect_false(is_log_log_concave(expr))

  ## Posynomial with unconstrained variable: not DGP
  posynomial <- 5.0 * x * y + 1.2 * y * y
  expect_false(is_dgp(posynomial))
  expect_false(is_log_log_convex(posynomial))
  expect_false(is_log_log_concave(posynomial))
})

## ── test_monomials ────────────────────────────────────────────────

## @cvxpy test_dgp.py::TestDgp::test_monomials
test_that("CVXPY parity: monomial curvature in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)

  ## 5 * x^0.1 * y^(-0.1) * z^3 is a monomial (log-log affine)
  monomial <- 5.0 * power(x, 0.1) * power(y, -0.1) * power(z, 3)
  expect_true(is_dgp(monomial))
  expect_true(is_log_log_convex(monomial))
  expect_true(is_log_log_concave(monomial))

  ## Multiplying by -1 breaks DGP
  neg_monomial <- monomial * (-1.0)
  expect_false(is_dgp(neg_monomial))
  expect_false(is_log_log_convex(neg_monomial))
  expect_false(is_log_log_concave(neg_monomial))
})

## ── test_maximum ──────────────────────────────────────────────────
## CVXPY: maximum of monomial, posynomial, posynomial^2 is DGP, log-log convex.

## @cvxpy test_dgp.py::TestDgp::test_maximum
test_that("CVXPY parity: maximum curvature in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)

  monomial <- 5.0 * power(x, 0.1) * power(y, -0.1) * power(z, 3)
  posynomial <- 5.0 * x * y + 1.2 * y * y
  another_posynomial <- posynomial * posynomial

  ## max(monomial, posynomial, posynomial^2) is DGP, log-log convex
  expr <- max_elemwise(monomial, posynomial, another_posynomial)
  expect_true(is_dgp(expr))
  expect_true(is_log_log_convex(expr))
  expect_false(is_log_log_concave(expr))

  ## posynomial * max(...) is still log-log convex
  expr2 <- posynomial * expr
  expect_true(is_dgp(expr2))
  expect_true(is_log_log_convex(expr2))
  expect_false(is_log_log_concave(expr2))

  ## posynomial * max(...) + max(...) is still log-log convex
  expr3 <- posynomial * expr2 + expr2
  expect_true(is_dgp(expr3))
  expect_true(is_log_log_convex(expr3))
})

## ── test_minimum ──────────────────────────────────────────────────
## CVXPY: minimum of monomial, 1/posynomial, 1/posynomial^2 is DGP, log-log concave.

## @cvxpy test_dgp.py::TestDgp::test_minimum
test_that("CVXPY parity: minimum curvature in DGP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)

  monomial <- 5.0 * power(x, 0.1) * power(y, -0.1) * power(z, 3)
  posynomial <- 5.0 * x * y + 1.2 * y * y
  another_posynomial <- posynomial * posynomial

  ## min(monomial, 1/posynomial, 1/another_posynomial) is DGP, log-log concave
  expr <- min_elemwise(monomial, 1 / posynomial, 1 / another_posynomial)
  expect_true(is_dgp(expr))
  expect_false(is_log_log_convex(expr))
  expect_true(is_log_log_concave(expr))

  ## (1/posynomial) * min(...) is still log-log concave
  expr2 <- (1 / posynomial) * expr
  expect_true(is_dgp(expr2))
  expect_false(is_log_log_convex(expr2))
  expect_true(is_log_log_concave(expr2))

  ## (expr)^2 with concave base: log-log concave (even power of concave monomial)
  expr3 <- expr2^2
  expect_true(is_dgp(expr3))
  expect_false(is_log_log_convex(expr3))
  expect_true(is_log_log_concave(expr3))
})

## ── test_constant ─────────────────────────────────────────────────

## @cvxpy test_dgp.py::TestDgp::test_constant
test_that("CVXPY parity: constant DGP status", {
  x <- Constant(1.0)
  expect_true(is_dgp(x))

  ## Negative constant is not DGP
  expect_false(is_dgp(-1.0 * x))
})

## ── test_geo_mean ─────────────────────────────────────────────────

## @cvxpy test_dgp.py::TestDgp::test_geo_mean
test_that("CVXPY parity: geo_mean curvature in DGP", {
  x <- Variable(3, pos = TRUE)
  p <- c(1, 2, 0.5)
  gm <- geo_mean(x, p)

  expect_true(is_dgp(gm))
  expect_true(is_log_log_affine(gm))
  expect_true(is_log_log_convex(gm))
  expect_true(is_log_log_concave(gm))
})

## ── test_geo_mean_scalar1 ─────────────────────────────────────────

## @cvxpy test_dgp.py::TestDgp::test_geo_mean_scalar1
test_that("CVXPY parity: geo_mean scalar1 — solve", {
  x <- Variable(1, pos = TRUE)
  p <- c(2)
  gm <- geo_mean(x, p)
  expect_true(is_dgp(gm))

  prob <- Problem(Maximize(gm), list(x == 2))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## geo_mean(x, [2]) with scalar x = x (weights normalize to 1)
  ## CVXPY: value = 2
  expect_equal(val, 2, tolerance = 1e-2)
})

## ── test_geo_mean_scalar2 ─────────────────────────────────────────

## @cvxpy test_dgp.py::TestDgp::test_geo_mean_scalar2
test_that("CVXPY parity: geo_mean scalar2 — scalar variable", {
  x <- Variable(pos = TRUE)
  p <- c(2)
  gm <- geo_mean(x, p)
  expect_true(is_dgp(gm))
})

## ── test_inv_prod ─────────────────────────────────────────────────
## CVXPY: inv_prod(x[0]) + inv_prod(x[:2]) with sum(x)==2, then
## compare inv_prod(x[:1]) vs inv_pos(x[0]).
## inv_prod of a scalar = inv_pos of that scalar.

## @cvxpy test_dgp.py::TestDgp::test_inv_prod
test_that("CVXPY parity: inv_prod consistency with inv_pos", {
  x <- Variable(2)

  ## Problem 1: inv_prod(x[1]) + inv_prod(x[1:2]) s.t. sum(x)==2
  ## In CVXPY: x[0] is scalar, x[:2] is the full vector
  prob1 <- Problem(
    Minimize(inv_prod(x[1, 1]) + inv_prod(x)),
    list(sum(x) == 2)
  )
  psolve(prob1, solver = "CLARABEL")

  ## Problem 2: inv_prod(x[1:1]) + inv_prod(x[1:2]) s.t. sum(x)==2
  ## In CVXPY: x[:1] is a 1-element slice
  prob2 <- Problem(
    Minimize(inv_prod(x[1, 1, drop = FALSE]) + inv_prod(x)),
    list(sum(x) == 2)
  )
  psolve(prob2, solver = "CLARABEL")

  ## Problem 3: inv_pos(x[1]) + inv_prod(x[1:2]) s.t. sum(x)==2
  prob3 <- Problem(
    Minimize(inv_pos(x[1, 1]) + inv_prod(x)),
    list(sum(x) == 2)
  )
  psolve(prob3, solver = "CLARABEL")

  ## Key assertion: prob2 and prob3 should have same value
  ## (inv_prod of a 1-element vector == inv_pos of that scalar)
  expect_equal(value(prob2), value(prob3), tolerance = 1e-3)
})

## ── test_builtin_sum ──────────────────────────────────────────────
## CVXPY: sum(x) for positive variables is log-log convex

## @cvxpy test_dgp.py::TestDgp::test_builtin_sum
test_that("CVXPY parity: builtin sum of positive variables is log-log convex", {
  x <- Variable(2, pos = TRUE)
  ## In R, sum(x) dispatches to SumEntries via Summary group
  expect_true(is_log_log_convex(sum(x)))
})

## ── test_gmatmul ──────────────────────────────────────────────────
## CVXPY: gmatmul errors on non-constant A and non-positive x.
## gmatmul(A, x) with constant A and positive x is DGP, log-log affine.

## @cvxpy test_dgp.py::TestDgp::test_gmatmul
test_that("CVXPY parity: gmatmul errors and curvature", {
  ## Error: A must be constant
  x <- Variable(2, pos = TRUE)
  A_var <- Variable(c(2, 2))
  expect_error(gmatmul(A_var, x), ".*constant.*")

  ## Error: x must be positive
  x_unc <- Variable(2)
  A <- matrix(1, 4, 2)
  expect_error(gmatmul(A, x_unc), ".*positive.*")

  ## Valid: constant A, positive x
  x3 <- Variable(3, pos = TRUE)
  A3 <- matrix(1, 4, 3)
  gm <- gmatmul(A3, x3)
  expect_true(is_dgp(gm))
  expect_true(is_log_log_affine(gm))
  expect_true(is_log_log_convex(gm))
  expect_true(is_log_log_concave(gm))
  expect_true(is_nonneg(gm))
  ## is_incr for first arg (idx=1 in R convention, CVXPY idx=0)
  expect_true(CVXR:::is_incr(gm, 1))
  expect_true(CVXR:::is_decr(gmatmul(-A3, x3), 1))

  ## Matrix x: (2,3) positive variable
  x_mat <- Variable(c(2, 3), pos = TRUE)
  A_mat <- matrix(c(2, 0, -1, 3), 2, 2)
  gm_mat <- gmatmul(A_mat, x_mat)
  expect_true(is_dgp(gm_mat))
  expect_true(is_log_log_affine(gm_mat))
  expect_true(is_log_log_convex(gm_mat))
  expect_true(is_log_log_concave(gm_mat))
  ## Mixed-sign A: not monotone
  expect_false(CVXR:::is_incr(gm_mat, 1))
  expect_false(CVXR:::is_decr(gm_mat, 1))
})

## ── test_power_sign ───────────────────────────────────────────────
## CVXPY: (positive variable)^1 is nonneg, not nonpos.

## @cvxpy test_dgp.py::TestDgp::test_power_sign
test_that("CVXPY parity: power sign for positive variable", {
  x <- Variable(pos = TRUE)
  expr <- x^1
  expect_true(is_nonneg(expr))
  expect_false(is_nonpos(expr))
})

## ── test_sparse_constant_not_allowed ──────────────────────────────
## CVXPY: sparse matrix constants are not log-log constant.

## @cvxpy test_dgp.py::TestDgp::test_sparse_constant_not_allowed
test_that("CVXPY parity: sparse constant is not log-log constant", {
  sparse_mat <- Constant(Matrix::sparseMatrix(
    i = c(1, 1), j = 1:2, x = c(1.0, 2.0), dims = c(1, 2)
  ))
  expect_false(CVXR:::is_log_log_constant(sparse_mat))
})
