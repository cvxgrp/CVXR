## CVXPY cone2cone (test_cone2cone.py) parity gap tests
## Tests ported from CVXPY commit 3b964472b (Release 1.8.1)
## Expected values verified against CVXPY source and sth_* helpers.
##
## Note: CVXR's cone2cone/affine2direct module is NOT IMPLEMENTED.
## The CVXPY tests (TestDualize, TestSlacks) verify the Dualize and Slacks
## internal reductions by manually building and solving the dual/slack
## reformulations. Since these reductions don't exist in CVXR, we test
## that the same underlying problems solve correctly through CVXR's
## standard solver chain, verifying identical optimal values and primals.
##
## TestPowND tests use PowConeND which IS implemented in CVXR.

# ======================================================================
# TestDualize — verify underlying problems solve correctly
# ======================================================================

## @cvxpy test_cone2cone.py::TestDualize::test_lp_1
test_that("cone2cone gap: Dualize LP 1 — typical LP", {
  skip_if_not_installed("clarabel")
  sth <- sth_lp_1()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestDualize::test_lp_2
test_that("cone2cone gap: Dualize LP 2 — typical LP", {
  skip_if_not_installed("clarabel")
  sth <- sth_lp_2()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestDualize::test_lp_3
test_that("cone2cone gap: Dualize LP 3 — unbounded LP", {
  skip_if_not_installed("clarabel")
  sth <- sth_lp_3()
  psolve(sth$prob, solver = "CLARABEL")
  ## Unbounded: status should be unbounded or value -Inf
  expect_true(status(sth$prob) %in% c("unbounded", "unbounded_inaccurate") ||
                value(sth$prob) == -Inf)
})

## @cvxpy test_cone2cone.py::TestDualize::test_lp_4
test_that("cone2cone gap: Dualize LP 4 — infeasible LP", {
  skip_if_not_installed("clarabel")
  sth <- sth_lp_4()
  psolve(sth$prob, solver = "CLARABEL")
  ## Infeasible: status should be infeasible or value +Inf
  expect_true(status(sth$prob) %in% c("infeasible", "infeasible_inaccurate") ||
                value(sth$prob) == Inf)
})

## @cvxpy test_cone2cone.py::TestDualize::test_lp_5
test_that("cone2cone gap: Dualize LP 5 — LP with redundant constraints", {
  skip_if_not_installed("clarabel")
  sth <- sth_lp_5()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestDualize::test_socp_0
test_that("cone2cone gap: Dualize SOCP 0", {
  skip_if_not_installed("clarabel")
  sth <- sth_socp_0()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestDualize::test_socp_1
test_that("cone2cone gap: Dualize SOCP 1", {
  skip_if_not_installed("clarabel")
  sth <- sth_socp_1()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestDualize::test_socp_2
test_that("cone2cone gap: Dualize SOCP 2", {
  skip_if_not_installed("clarabel")
  sth <- sth_socp_2()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestDualize::test_socp_3_axis_0
test_that("cone2cone gap: Dualize SOCP 3 axis=0", {
  skip_if_not_installed("clarabel")
  ## CVXPY axis=0 -> R axis=2 (column-wise)
  sth <- sth_socp_3_ax0()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestDualize::test_socp_3_axis_1
test_that("cone2cone gap: Dualize SOCP 3 axis=1", {
  skip_if_not_installed("clarabel")
  ## CVXPY axis=1 -> R axis=1 (row-wise)
  sth <- sth_socp_3_ax1()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestDualize::test_expcone_1
test_that("cone2cone gap: Dualize ExpCone 1", {
  skip_if_not_installed("clarabel")
  sth <- sth_expcone_1()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestDualize::test_expcone_socp_1
test_that("cone2cone gap: Dualize ExpCone+SOCP 1", {
  skip_if_not_installed("scs")
  sth <- sth_expcone_socp_1()
  psolve(sth$prob, solver = "SCS")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-2)
})

## @cvxpy test_cone2cone.py::TestDualize::test_pcp_2
test_that("cone2cone gap: Dualize PowCone 2", {
  skip_if_not_installed("scs")
  sth <- sth_pcp_2()
  psolve(sth$prob, solver = "SCS")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-2)
})

# ======================================================================
# TestSlacks — verify underlying problems solve correctly
# ======================================================================

## @cvxpy test_cone2cone.py::TestSlacks::test_lp_2
test_that("cone2cone gap: Slacks LP 2 — typical LP", {
  skip_if_not_installed("clarabel")
  sth <- sth_lp_2()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
  ## Verify primal
  x_val <- as.numeric(value(variables(sth$prob)[[1]]))
  expect_equal(x_val, sth$expect_x, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_lp_3
test_that("cone2cone gap: Slacks LP 3 — unbounded LP", {
  skip_if_not_installed("clarabel")
  sth <- sth_lp_3()
  psolve(sth$prob, solver = "CLARABEL")
  expect_true(status(sth$prob) %in% c("unbounded", "unbounded_inaccurate") ||
                value(sth$prob) == -Inf)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_lp_4
test_that("cone2cone gap: Slacks LP 4 — infeasible LP", {
  skip_if_not_installed("clarabel")
  sth <- sth_lp_4()
  psolve(sth$prob, solver = "CLARABEL")
  expect_true(status(sth$prob) %in% c("infeasible", "infeasible_inaccurate") ||
                value(sth$prob) == Inf)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_socp_2
test_that("cone2cone gap: Slacks SOCP 2", {
  skip_if_not_installed("clarabel")
  sth <- sth_socp_2()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
  x_val <- as.numeric(value(variables(sth$prob)[[1]]))
  expect_equal(x_val, sth$expect_x, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_socp_3
test_that("cone2cone gap: Slacks SOCP 3 (both axes)", {
  skip_if_not_installed("clarabel")

  ## axis=0 (R: sth_socp_3_ax0)
  sth0 <- sth_socp_3_ax0()
  psolve(sth0$prob, solver = "CLARABEL")
  expect_equal(value(sth0$prob), sth0$expect_obj, tolerance = 1e-3)
  x_val0 <- as.numeric(value(variables(sth0$prob)[[1]]))
  expect_equal(x_val0, sth0$expect_x, tolerance = 1e-3)

  ## axis=1 (R: sth_socp_3_ax1)
  sth1 <- sth_socp_3_ax1()
  psolve(sth1$prob, solver = "CLARABEL")
  expect_equal(value(sth1$prob), sth1$expect_obj, tolerance = 1e-3)
  x_val1 <- as.numeric(value(variables(sth1$prob)[[1]]))
  expect_equal(x_val1, sth1$expect_x, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_expcone_1
test_that("cone2cone gap: Slacks ExpCone 1", {
  skip_if_not_installed("clarabel")
  sth <- sth_expcone_1()
  psolve(sth$prob, solver = "CLARABEL")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_expcone_socp_1
test_that("cone2cone gap: Slacks ExpCone+SOCP 1", {
  skip_if_not_installed("scs")
  sth <- sth_expcone_socp_1()
  psolve(sth$prob, solver = "SCS")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-2)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_pcp_1
test_that("cone2cone gap: Slacks PowCone 1", {
  skip_if_not_installed("scs")
  sth <- sth_pcp_1()
  psolve(sth$prob, solver = "SCS")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-2)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_pcp_2
test_that("cone2cone gap: Slacks PowCone 2", {
  skip_if_not_installed("scs")
  sth <- sth_pcp_2()
  psolve(sth$prob, solver = "SCS")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-2)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_mi_lp_1
test_that("cone2cone gap: Slacks MI LP 1", {
  skip_if_not_installed("highs")
  sth <- sth_mi_lp_1()
  psolve(sth$prob, solver = "HIGHS")
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_mi_socp_1
test_that("cone2cone gap: Slacks MI SOCP 1", {
  skip("Known bug in ECOS BB — same as CVXPY skip")
  sth <- sth_mi_socp_1()
  psolve(sth$prob)
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_cone2cone.py::TestSlacks::test_mi_socp_2
test_that("cone2cone gap: Slacks MI SOCP 2", {
  ## Need a MIP-capable SOCP solver
  has_mi <- FALSE
  if (requireNamespace("gurobi", quietly = TRUE)) has_mi <- TRUE
  if (requireNamespace("Rmosek", quietly = TRUE)) has_mi <- TRUE
  skip_if(!has_mi, "No mixed-integer SOCP solver installed")

  sth <- sth_mi_socp_2()
  psolve(sth$prob)
  expect_equal(value(sth$prob), sth$expect_obj, tolerance = 1e-3)
})

# ======================================================================
# TestPowND — PowConeND problems
# ======================================================================

## @cvxpy test_cone2cone.py::TestPowND::test_pcp_3a
test_that("cone2cone gap: PowND pcp_3 axis=0", {
  skip_if_not_installed("clarabel")

  ## Reformulation of pcp_2 using PowConeND
  ## max  x3 + x4 - x0
  ## s.t. x0 + x1 + x2/2 == 2,
  ##      (W, z) in PowND(alpha, axis=0)
  ## CVXPY axis=0 -> R axis=2 (column-wise)
  x <- Variable(3)
  expect_x <- c(0.06393515, 0.78320961, 2.30571048)
  hypos <- Variable(2)
  objective <- Maximize(sum_entries(hypos) - x[1])

  W <- bmat(list(list(x[1], x[3]),
                 list(x[2], 1.0)))
  alpha <- matrix(c(0.2, 0.4, 0.8, 0.6), 2, 2, byrow = TRUE)

  constraints <- list(
    x[1] + x[2] + 0.5 * x[3] == 2,
    PowConeND(W, hypos, alpha, axis = 2L)
  )
  prob <- Problem(objective, constraints)
  psolve(prob, solver = "CLARABEL")
  expect_equal(value(prob), 1.8073, tolerance = 1e-2)
  x_val <- as.numeric(value(x))
  expect_equal(x_val, expect_x, tolerance = 1e-2)
})

## @cvxpy test_cone2cone.py::TestPowND::test_pcp_3b
test_that("cone2cone gap: PowND pcp_3 axis=1", {
  skip_if_not_installed("clarabel")

  ## Same problem but transposed: axis=1 (CVXPY axis=1 -> R axis=1)
  ## ConeMatrixStuffing now transposes PowConeND(axis=1) -> PowConeND(axis=2)
  ## before extracting alpha columns (matching CVXPY lines 382-388).
  x <- Variable(3)
  expect_x <- c(0.06393515, 0.78320961, 2.30571048)
  hypos <- Variable(2)
  objective <- Maximize(sum_entries(hypos) - x[1])

  ## Transpose W and alpha for axis=1
  W <- bmat(list(list(x[1], x[2]),
                 list(x[3], 1.0)))
  alpha <- matrix(c(0.2, 0.8, 0.4, 0.6), 2, 2, byrow = TRUE)

  constraints <- list(
    x[1] + x[2] + 0.5 * x[3] == 2,
    PowConeND(W, hypos, alpha, axis = 1L)
  )
  prob <- Problem(objective, constraints)
  psolve(prob, solver = "CLARABEL")
  expect_equal(value(prob), 1.8073, tolerance = 1e-2)
  x_val <- as.numeric(value(x))
  expect_equal(x_val, expect_x, tolerance = 1e-2)
})

## @cvxpy test_cone2cone.py::TestPowND::test_pcp_4a
test_that("cone2cone gap: PowND pcp_4a (Fisher market CEEI)", {
  skip_if_not_installed("clarabel")

  ## Fisher market equilibrium: competitive equilibrium from equal incomes
  set.seed(0)
  n_buyer <- 4L
  n_items <- 6L
  V <- matrix(runif(n_buyer * n_items), n_buyer, n_items)
  X <- Variable(c(n_buyer, n_items), nonneg = TRUE)
  ## u = sum(V * X, axis=1) => R: apply(V*X, 1, sum) => sum over columns per row
  u <- sum_entries(V * X, axis = 1L)

  b <- rep(1, n_buyer) / n_buyer  # CEEI: equal budgets

  ## First solve as log formulation to get reference
  log_objective <- Maximize(sum_entries(b * log(u)))
  log_cons <- list(sum_entries(X, axis = 2L) <= 1)
  log_prob <- Problem(log_objective, log_cons)
  psolve(log_prob, solver = "CLARABEL")
  expect_X <- value(X)
  log_opt <- value(log_prob)

  ## Power cone formulation
  z <- Variable()
  pow_objective <- Maximize(z)
  pow_cons <- list(
    sum_entries(X, axis = 2L) <= 1,
    PowConeND(W = u, z = z, alpha = b)
  )
  pow_prob <- Problem(pow_objective, pow_cons)
  psolve(pow_prob, solver = "CLARABEL")

  ## PowCone objective should equal exp(log_opt)
  expect_equal(value(pow_prob), exp(log_opt), tolerance = 1e-2)
})

## @cvxpy test_cone2cone.py::TestPowND::test_pcp_4b
test_that("cone2cone gap: PowND pcp_4b (Fisher market non-CEEI)", {
  skip_if_not_installed("clarabel")

  ## Fisher market with non-equal budgets
  set.seed(0)
  n_buyer <- 4L
  n_items <- 6L
  V <- matrix(runif(n_buyer * n_items), n_buyer, n_items)
  X <- Variable(c(n_buyer, n_items), nonneg = TRUE)
  u <- sum_entries(V * X, axis = 1L)

  b <- c(0.3, 0.15, 0.2, 0.35)

  ## Reference via log formulation
  log_objective <- Maximize(sum_entries(b * log(u)))
  log_cons <- list(sum_entries(X, axis = 2L) <= 1)
  log_prob <- Problem(log_objective, log_cons)
  psolve(log_prob, solver = "CLARABEL")
  log_opt <- value(log_prob)

  ## Power cone formulation
  z <- Variable()
  pow_objective <- Maximize(z)
  pow_cons <- list(
    sum_entries(X, axis = 2L) <= 1,
    PowConeND(W = u, z = z, alpha = b)
  )
  pow_prob <- Problem(pow_objective, pow_cons)
  psolve(pow_prob, solver = "CLARABEL")

  expect_equal(value(pow_prob), exp(log_opt), tolerance = 1e-2)
})
