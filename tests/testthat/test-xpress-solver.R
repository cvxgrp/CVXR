## Tests for XPRESS solver integration
## Mirrors CVXPY test_conic_solvers.py::TestXPRESS
## All tests skip cleanly when xpress is not available.

# ══════════════════════════════════════════════════════════════════
# Setup: skip entire file if XPRESS not available
# ══════════════════════════════════════════════════════════════════
require_solver("XPRESS")

# ══════════════════════════════════════════════════════════════════
# Unit tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("XPRESS_SOLVER constant is exported", {
  expect_equal(CVXR::XPRESS_SOLVER, "XPRESS")
})

## @cvxpy NONE
test_that("XPRESS appears in installed_solvers when available", {
  expect_true("XPRESS" %in% installed_solvers())
})

# ══════════════════════════════════════════════════════════════════
# Standard LP tests (lp_0 through lp_4)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_lp_0
## lp_0: min ||x||_1 + 1 s.t. x == 0 => obj=1, x=[0,0]
test_that("XPRESS lp_0: norm1 with equality", {
  x <- Variable(2)
  prob <- Problem(Minimize(cvxr_norm(x, 1) + 1.0), list(x == 0))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, 1.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0, 0), ncol = 1), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_lp_1
## lp_1: min -4x0 - 5x1 s.t. 2x0+x1<=3, x0+2x1<=3, x>=0 => obj=-9, x=[1,1]
test_that("XPRESS lp_1: LP with dual values", {
  x <- Variable(2)
  c1 <- (2 * x[1] + x[2] <= 3)
  c2 <- (x[1] + 2 * x[2] <= 3)
  c3 <- (x[1] >= 0)
  c4 <- (x[2] >= 0)
  prob <- Problem(Minimize(-4 * x[1] - 5 * x[2]), list(c1, c2, c3, c4))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, -9.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(1, 1), ncol = 1), tolerance = 1e-3)
  ## Dual values: c1->1, c2->2, c3->0, c4->0
  expect_equal(as.numeric(dual_value(c1)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2)), 2.0, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_lp_2
## lp_2: min x0 + 0.5*x1 s.t. x0>=-100, x0<=-10, x1==1 => obj=-99.5
test_that("XPRESS lp_2: LP with equality and dual values", {
  x <- Variable(2)
  c1 <- (x[1] >= -100)
  c2 <- (x[1] <= -10)
  c3 <- (x[2] == 1)
  prob <- Problem(Minimize(x[1] + 0.5 * x[2]), list(c1, c2, c3))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, -99.5, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(-100, 1), ncol = 1), tolerance = 1e-3)
  ## Dual values: c1->1, c2->0, c3->-0.5
  expect_equal(as.numeric(dual_value(c1)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c3)), -0.5, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_lp_3
## lp_3: unbounded — min sum(x) s.t. x <= 1 => obj=-Inf
test_that("XPRESS lp_3: unbounded LP", {
  x <- Variable(5)
  prob <- Problem(Minimize(sum_entries(x)), list(x <= 1))
  val <- psolve(prob, solver = "XPRESS")
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate"))
  expect_true(is.infinite(val) && val < 0)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_lp_4
## lp_4: infeasible — min sum(x) s.t. x <= 0, x >= 1 => infeasible
test_that("XPRESS lp_4: infeasible LP", {
  x <- Variable(5)
  prob <- Problem(Minimize(sum_entries(x)), list(x <= 0, x >= 1))
  val <- psolve(prob, solver = "XPRESS")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
  expect_true(is.infinite(val) && val > 0)
})

# ══════════════════════════════════════════════════════════════════
# Standard SOCP tests (socp_0 through socp_2)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_socp_0
## socp_0: min ||x||_2 + 1 s.t. x == 0 => obj=1, x=[0,0]
test_that("XPRESS socp_0: norm2 with equality", {
  x <- Variable(2)
  prob <- Problem(Minimize(cvxr_norm(x, 2) + 1), list(x == 0))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, 1.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0, 0), ncol = 1), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_socp_1
## socp_1: min 3x0+2x1+x2 s.t. ||x||<=y, x0+x1+3x2>=1, y<=5
test_that("XPRESS socp_1: SOCP with SOC constraint", {
  x <- Variable(3)
  y <- Variable()
  prob <- Problem(Minimize(3 * x[1] + 2 * x[2] + x[3]),
                  list(cvxr_norm(x, 2) <= y,
                       x[1] + x[2] + 3 * x[3] >= 1.0,
                       y <= 5))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, -13.5486, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 5.0, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_socp_2
## socp_2: SOCP reformulation of lp_1
test_that("XPRESS socp_2: SOCP reformulation of LP", {
  x <- Variable(2)
  prob <- Problem(Minimize(-4 * x[1] - 5 * x[2]),
                  list(2 * x[1] + x[2] <= 3,
                       cvxr_norm(Reshape(x[1] + 2 * x[2], c(1, 1)), 2) <= 3,
                       x[1] >= 0,
                       x[2] >= 0))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, -9.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(1, 1), ncol = 1), tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# Standard MI-LP tests (mi_lp_0 through mi_lp_5)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_mi_lp_0
## mi_lp_0: min ||x||_1 + 1 s.t. x == bool_var, bool_var == 0
test_that("XPRESS mi_lp_0: boolean with norm1", {
  x <- Variable(2)
  bv <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(cvxr_norm(x, 1) + 1.0),
                  list(x == bv, bv == 0))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, 1.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0, 0), ncol = 1), tolerance = 1e-3)
  expect_equal(as.numeric(value(bv)), 0, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_mi_lp_1
## mi_lp_1: min -4x0-5x1 s.t. 2x0+x1<=intvar, x0+2x1<=3*boolvar, ...
test_that("XPRESS mi_lp_1: mixed boolean + integer LP", {
  x <- Variable(2)
  boolvar <- Variable(boolean = TRUE)
  intvar <- Variable(integer = TRUE)
  prob <- Problem(Minimize(-4 * x[1] - 5 * x[2]),
                  list(2 * x[1] + x[2] <= intvar,
                       x[1] + 2 * x[2] <= 3 * boolvar,
                       x >= 0,
                       intvar == 3 * boolvar,
                       intvar == 3))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, -9.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(1, 1), ncol = 1), tolerance = 1e-3)
  expect_equal(as.numeric(value(boolvar)), 1, tolerance = 1e-3)
  expect_equal(as.numeric(value(intvar)), 3, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_mi_lp_2
## mi_lp_2: 50-item knapsack => obj=8373
test_that("XPRESS mi_lp_2: knapsack problem", {
  coeffs <- matrix(c(
     1, 94, 485,  2, 506, 326,  3, 416, 248,  4, 992, 421,  5, 649, 322,
     6, 237, 795,  7, 457, 43,   8, 815, 845,  9, 446, 955, 10, 422, 252,
    11, 791, 9,   12, 359, 901, 13, 667, 122, 14, 598, 94,  15, 7,   738,
    16, 544, 574, 17, 334, 715, 18, 766, 882, 19, 994, 367, 20, 893, 984,
    21, 633, 299, 22, 131, 433, 23, 428, 682, 24, 700, 72,  25, 617, 874,
    26, 874, 138, 27, 720, 856, 28, 419, 145, 29, 794, 995, 30, 196, 529,
    31, 997, 199, 32, 116, 277, 33, 908, 97,  34, 539, 719, 35, 707, 242,
    36, 569, 107, 37, 537, 122, 38, 931, 70,  39, 726, 98,  40, 487, 600,
    41, 772, 645, 42, 513, 267, 43, 81,  972, 44, 943, 895, 45, 58,  213,
    46, 303, 748, 47, 764, 487, 48, 536, 923, 49, 724, 29,  50, 789, 674
  ), ncol = 3, byrow = TRUE)
  n <- 50L
  profits <- coeffs[, 2]
  weights <- coeffs[, 3]
  X <- Variable(n, boolean = TRUE)
  prob <- Problem(Maximize(sum_entries(Multiply(profits, X))),
                  list(sum_entries(Multiply(weights, X)) <= 995))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, 8373, tolerance = 1)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_mi_lp_3
## mi_lp_3: infeasible boolean => obj=+Inf
test_that("XPRESS mi_lp_3: infeasible boolean LP", {
  x <- Variable(4, boolean = TRUE)
  prob <- Problem(Maximize(Constant(1)),
                  list(x[1] + x[2] + x[3] + x[4] <= 2,
                       x[1] + x[2] + x[3] + x[4] >= 2,
                       x[1] + x[2] <= 1,
                       x[1] + x[3] <= 1,
                       x[1] + x[4] <= 1,
                       x[3] + x[4] <= 1,
                       x[2] + x[4] <= 1,
                       x[2] + x[3] <= 1))
  val <- psolve(prob, solver = "XPRESS")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_mi_lp_5
## mi_lp_5: infeasible boolean — Sagemath ticket #31962
test_that("XPRESS mi_lp_5: infeasible boolean LP (Sagemath)", {
  z <- Variable(11, boolean = TRUE)
  prob <- Problem(Minimize(Constant(0)),
                  list(z[3] + z[2] == 1,
                       z[5] + z[4] == 1,
                       z[7] + z[6] == 1,
                       z[9] + z[8] == 1,
                       z[11] + z[10] == 1,
                       z[5] + z[2] <= 1,
                       z[3] + z[4] <= 1,
                       z[7] + z[3] <= 1,
                       z[2] + z[6] <= 1,
                       z[9] + z[7] <= 1,
                       z[6] + z[8] <= 1,
                       z[11] + z[9] <= 1,
                       z[8] + z[10] <= 1,
                       z[10] + z[5] <= 1,
                       z[4] + z[11] <= 1))
  val <- psolve(prob, solver = "XPRESS")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

# ══════════════════════════════════════════════════════════════════
# Standard MI-SOCP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_mi_socp_1
## mi_socp_1: min 3x0+2x1+x2+y0+2y1 s.t. ||x||<=y0, ||x||<=y1, ...
test_that("XPRESS mi_socp_1: mixed-integer SOCP", {
  x <- Variable(3)
  y <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(3 * x[1] + 2 * x[2] + x[3] + y[1] + 2 * y[2]),
                  list(cvxr_norm(x, 2) <= y[1],
                       cvxr_norm(x, 2) <= y[2],
                       x[1] + x[2] + 3 * x[3] >= 0.1,
                       y <= 5))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, 0.2136, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), c(1, 1), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_mi_socp_2
## mi_socp_2: SOCP reformulation of mi_lp_1
test_that("XPRESS mi_socp_2: MI-SOCP reformulation of MI-LP", {
  x <- Variable(2)
  bool_var <- Variable(boolean = TRUE)
  int_var <- Variable(integer = TRUE)
  prob <- Problem(Minimize(-4 * x[1] - 5 * x[2]),
                  list(2 * x[1] + x[2] <= int_var,
                       (x[1] + 2 * x[2])^2 <= 9 * bool_var,
                       x >= 0,
                       int_var == 3 * bool_var,
                       int_var == 3))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, -9.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(1, 1), ncol = 1), tolerance = 1e-3)
  expect_equal(as.numeric(value(bool_var)), 1, tolerance = 1e-3)
  expect_equal(as.numeric(value(int_var)), 3, tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# Warm-start and solver params
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_warm_start
test_that("XPRESS warm-start with parameter changes", {
  x <- Variable(2)

  A <- Parameter(c(2, 2))
  b <- Parameter(2)
  h <- Parameter(2)
  cc <- Parameter(2)

  value(A) <- matrix(c(1, 0, 0, 0), 2, 2)
  value(b) <- c(1, 0)
  value(h) <- c(2, 2)
  value(cc) <- c(1, 1)

  prob <- Problem(Maximize(cc[1] * x[1] + cc[2] * x[2]),
                  list(x[1] <= h[1], x[2] <= h[2], A %*% x == b))
  val <- psolve(prob, solver = "XPRESS", warm_start = TRUE)
  expect_equal(val, 3.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(1, 2), ncol = 1), tolerance = 1e-3)

  ## Change A and b
  value(A) <- matrix(c(0, 0, 0, 1), 2, 2)
  value(b) <- c(0, 1)
  x2 <- Variable(2)
  prob2 <- Problem(Maximize(cc[1] * x2[1] + cc[2] * x2[2]),
                   list(x2[1] <= h[1], x2[2] <= h[2], A %*% x2 == b))
  val2 <- psolve(prob2, solver = "XPRESS", warm_start = TRUE)
  expect_equal(val2, 3.0, tolerance = 1e-4)
  expect_equal(value(x2), matrix(c(2, 1), ncol = 1), tolerance = 1e-3)

  ## Change h
  value(A) <- matrix(c(1, 0, 0, 0), 2, 2)
  value(b) <- c(1, 0)
  value(h) <- c(1, 1)
  x3 <- Variable(2)
  prob3 <- Problem(Maximize(cc[1] * x3[1] + cc[2] * x3[2]),
                   list(x3[1] <= h[1], x3[2] <= h[2], A %*% x3 == b))
  val3 <- psolve(prob3, solver = "XPRESS", warm_start = TRUE)
  expect_equal(val3, 2.0, tolerance = 1e-4)
  expect_equal(value(x3), matrix(c(1, 1), ncol = 1), tolerance = 1e-3)

  ## Change c
  value(h) <- c(2, 2)
  value(cc) <- c(2, 1)
  x4 <- Variable(2)
  prob4 <- Problem(Maximize(cc[1] * x4[1] + cc[2] * x4[2]),
                   list(x4[1] <= h[1], x4[2] <= h[2], A %*% x4 == b))
  val4 <- psolve(prob4, solver = "XPRESS", warm_start = TRUE)
  expect_equal(val4, 4.0, tolerance = 1e-4)
  expect_equal(value(x4), matrix(c(1, 2), ncol = 1), tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_params
test_that("XPRESS solver params passthrough", {
  set.seed(0)
  n <- 10L; m <- 4L
  A_mat <- matrix(rnorm(m * n), m, n)
  x0 <- rnorm(n)
  y_vec <- A_mat %*% x0

  z <- Variable(n)
  prob <- Problem(Minimize(cvxr_norm(z, 1)),
                  list(A_mat %*% z == y_vec))
  ## Pass solver-specific controls
  val <- psolve(prob, solver = "XPRESS", LPITERLIMIT = 1000L, MAXTIME = 1000L)
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_iis_none
## DEFERRED: IIS not implemented in R XPRESS interface

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_iis_full
## DEFERRED: IIS not implemented in R XPRESS interface

## @cvxpy test_conic_solvers.py::TestXPRESS::test_xpress_lp_bound_attr
## N/A: BOUNDED_VARIABLES is deferred in CVXR (Tier 2)

# ══════════════════════════════════════════════════════════════════
# QP-path tests (no CVXPY counterpart — QP path is R-specific)
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("XPRESS QP via QP path", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + x[1]),
                  list(x[1] + x[2] == 1, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, 0.875, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0.25, 0.75), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("XPRESS MIQP (integer + quadratic)", {
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(sum_squares(x) + x[1]),
                  list(x[1] + x[2] == 2, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "XPRESS")
  ## (0,2)->4, (1,1)->3, (2,0)->6
  expect_equal(val, 3.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(1, 1), ncol = 1), tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# Dual value tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("XPRESS LP duals", {
  x <- Variable(2)
  c1 <- (x[1] + x[2] >= 3)
  c2 <- (x[1] >= 0)
  c3 <- (x[2] >= 0)
  prob <- Problem(Minimize(x[1] + 2 * x[2]), list(c1, c2, c3))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, 3.0, tolerance = 1e-4)
  d1 <- dual_value(c1)
  expect_true(!is.null(d1))
  expect_equal(as.numeric(d1), 1.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("XPRESS QP duals", {
  x <- Variable(2)
  c1 <- (x[1] + x[2] == 1)
  c2 <- (x[1] >= 0)
  c3 <- (x[2] >= 0)
  prob <- Problem(Minimize(sum_squares(x) + x[1]), list(c1, c2, c3))
  val <- psolve(prob, solver = "XPRESS")
  expect_equal(val, 0.875, tolerance = 1e-4)
  d1 <- dual_value(c1)
  expect_true(!is.null(d1))
  expect_true(abs(as.numeric(d1)) > 0.01)
})

# ══════════════════════════════════════════════════════════════════
# Cross-validation: XPRESS vs Clarabel
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("XPRESS LP matches Clarabel", {
  x <- Variable(2)
  prob_x <- Problem(Minimize(x[1] + 2 * x[2]),
                    list(x[1] + x[2] >= 3, x[1] >= 0, x[2] >= 0))
  val_x <- psolve(prob_x, solver = "XPRESS")

  y <- Variable(2)
  prob_c <- Problem(Minimize(y[1] + 2 * y[2]),
                    list(y[1] + y[2] >= 3, y[1] >= 0, y[2] >= 0))
  val_c <- psolve(prob_c, solver = "CLARABEL")

  expect_equal(val_x, val_c, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("XPRESS QP matches Clarabel", {
  x <- Variable(2)
  prob_x <- Problem(Minimize(sum_squares(x) + x[1]),
                    list(x[1] + x[2] == 1, x[1] >= 0, x[2] >= 0))
  val_x <- psolve(prob_x, solver = "XPRESS")

  y <- Variable(2)
  prob_c <- Problem(Minimize(sum_squares(y) + y[1]),
                    list(y[1] + y[2] == 1, y[1] >= 0, y[2] >= 0))
  val_c <- psolve(prob_c, solver = "CLARABEL")

  expect_equal(val_x, val_c, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("XPRESS SOCP matches Clarabel", {
  x <- Variable(2)
  prob_x <- Problem(Minimize(x[1] + x[2]),
                    list(cvxr_norm(x, 2) <= 1))
  val_x <- psolve(prob_x, solver = "XPRESS")

  y <- Variable(2)
  prob_c <- Problem(Minimize(y[1] + y[2]),
                    list(cvxr_norm(y, 2) <= 1))
  val_c <- psolve(prob_c, solver = "CLARABEL")

  expect_equal(val_x, val_c, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("XPRESS LP dual signs match Clarabel", {
  x <- Variable(2)
  c1 <- (x[1] + x[2] >= 3)
  c2 <- (x[1] >= 0)
  c3 <- (x[2] >= 0)
  prob_x <- Problem(Minimize(x[1] + 2 * x[2]), list(c1, c2, c3))
  psolve(prob_x, solver = "XPRESS")
  d1_x <- dual_value(c1)

  y <- Variable(2)
  c1c <- (y[1] + y[2] >= 3)
  c2c <- (y[1] >= 0)
  c3c <- (y[2] >= 0)
  prob_c <- Problem(Minimize(y[1] + 2 * y[2]), list(c1c, c2c, c3c))
  psolve(prob_c, solver = "CLARABEL")
  d1_c <- dual_value(c1c)

  expect_equal(as.numeric(d1_x), as.numeric(d1_c), tolerance = 1e-3)
})
