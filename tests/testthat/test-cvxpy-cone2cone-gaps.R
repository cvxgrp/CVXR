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

# ======================================================================
# TestPSDUtils — round-trip of the shared svec packer/unpacker
# (utilities/psd_utils.R; ported at 1.9.2 under ADR D_19.6)
# ======================================================================

.random_psd <- function(n, seed) {
  set.seed(seed)
  A <- matrix(rnorm(n * n), n, n)
  A %*% t(A)
}

## @cvxpy test_cone2cone.py::TestPSDUtils::test_tri_to_full_round_trip
test_that("psd_utils: psd_format_mat -> tri_to_full recovers the symmetrized matrix", {
  n <- 4L
  X <- Variable(c(n, n), symmetric = TRUE)
  con <- PSD(X)
  M_val <- .random_psd(n, 42L)

  for (tri_kind in c(TriangleKind$LOWER, TriangleKind$UPPER)) {
    for (sqrt2 in c(TRUE, FALSE)) {
      M <- psd_format_mat(con, tri_kind, sqrt2)
      svec <- as.vector(M %*% as.vector(M_val))   # as.vector() is F-order in R
      recovered <- tri_to_full(svec, n, tri_kind, sqrt2)
      expect_equal(recovered, M_val, tolerance = 1e-12)
    }
  }
})

## @cvxpy test_cone2cone.py::TestPSDUtils::test_tri_to_full_n1
test_that("psd_utils: 1x1 matrix round-trips (edge case: no off-diagonals)", {
  X <- Variable(c(1L, 1L), symmetric = TRUE)
  con <- PSD(X)
  M_val <- matrix(5.0, 1L, 1L)

  for (tri_kind in c(TriangleKind$LOWER, TriangleKind$UPPER)) {
    for (sqrt2 in c(TRUE, FALSE)) {
      M <- psd_format_mat(con, tri_kind, sqrt2)
      svec <- as.vector(M %*% as.vector(M_val))
      recovered <- tri_to_full(svec, 1L, tri_kind, sqrt2)
      expect_equal(recovered, M_val, tolerance = 1e-12)
    }
  }
})

## Guards the CRAN 1.9.1 bug class (PSD duals scrambled for n >= 3, 893ca4b):
## the triangle kinds are NOT interchangeable.  The upstream round-trip tests
## above cannot catch this on their own -- pack and unpack with the SAME kind
## and any consistent mis-enumeration cancels out.  So pack the lower triangle
## and unpack it as the upper one: the entries must land in different places.
## (Not a transpose of the recovered matrix: `right` is symmetric, so the wrong
## kind permutes the off-diagonals rather than transposing the result.)
## @cvxpy NONE
test_that("psd_utils: the two triangle kinds are not interchangeable for n >= 3", {
  n <- 3L
  X <- Variable(c(n, n), symmetric = TRUE)
  con <- PSD(X)
  M_val <- .random_psd(n, 11L)

  M_lower <- psd_format_mat(con, TriangleKind$LOWER, TRUE)
  svec <- as.vector(M_lower %*% as.vector(M_val))

  right <- tri_to_full(svec, n, TriangleKind$LOWER, TRUE)
  wrong <- tri_to_full(svec, n, TriangleKind$UPPER, TRUE)
  expect_equal(right, M_val, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(wrong, M_val)))
})

# ======================================================================
# ExactCone2Cone / PSD -> SvecPSD  (ADR D_19.6 step B)
# ======================================================================

## @cvxpy test_cone2cone.py::TestExactApproxCone2Cone::test_exact_cone_conversions_map
test_that("EXACT_CONE_CONVERSIONS maps PSD to SvecPSD", {
  srcs <- lapply(EXACT_CONE_CONVERSIONS, `[[`, "source")
  expect_true(any(vapply(srcs, identical, logical(1L), PSD)))
  entry <- Filter(function(e) identical(e$source, PSD), EXACT_CONE_CONVERSIONS)[[1L]]
  expect_true(any(vapply(entry$targets, identical, logical(1L), SvecPSD)))
  ## PowCone3D is not a conversion source (it is a target of PowConeND upstream).
  expect_false(any(vapply(srcs, identical, logical(1L), PowCone3D)))
})

## @cvxpy test_cone2cone.py::TestExactApproxCone2Cone::test_exact_cone_conversions_is_dag
test_that("EXACT_CONE_CONVERSIONS is a DAG (no cone reaches itself)", {
  reach <- lapply(EXACT_CONE_CONVERSIONS, function(e) e$targets)
  srcs  <- lapply(EXACT_CONE_CONVERSIONS, `[[`, "source")
  repeat {
    changed <- FALSE
    for (i in seq_along(srcs)) {
      for (t in reach[[i]]) {
        j <- which(vapply(srcs, identical, logical(1L), t))
        if (length(j) == 0L) next
        for (tt in reach[[j]]) {
          if (!any(vapply(reach[[i]], identical, logical(1L), tt))) {
            reach[[i]] <- c(reach[[i]], list(tt))
            changed <- TRUE
          }
        }
      }
    }
    if (!changed) break
  }
  for (i in seq_along(srcs)) {
    expect_false(any(vapply(reach[[i]], identical, logical(1L), srcs[[i]])),
                 info = "cycle detected in EXACT_CONE_CONVERSIONS")
  }
})

## @cvxpy test_cone2cone.py::TestExactApproxCone2Cone::test_exact_cone2cone_target_cones_filtering
test_that("ExactCone2Cone(target_cones=) converts only the cones asked for", {
  X <- Variable(c(3L, 3L), symmetric = TRUE)
  con <- PSD(X)
  ctx <- SolverInfo(solver_name = "CLARABEL",
                    psd_triangle_kind = TriangleKind$UPPER,
                    psd_sqrt2_scaling = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)), list(con))

  ## Asked for PSD: converted.
  red <- ExactCone2Cone(target_cones = list(PSD), solver_context = ctx)
  out <- reduction_apply(red, prob)
  expect_true(.s7_is(out[[1L]]@constraints[[1L]], SvecPSD))

  ## Asked for something else: left alone.
  red2 <- ExactCone2Cone(target_cones = list(SOC), solver_context = ctx)
  out2 <- reduction_apply(red2, prob)
  expect_true(.s7_is(out2[[1L]]@constraints[[1L]], PSD))
})

## The svec length is n(n+1)/2 per cone, NOT n^2 -- this is what the solver's
## cone dims are checked against, and getting it wrong is a hard solver failure
## ("cone dimensions ... not equal to num rows in A").
## @cvxpy NONE
test_that("SvecPSD reports the packed size, and the chain packs PSD for svec solvers", {
  n <- 4L
  X <- Variable(c(n, n), symmetric = TRUE)
  sv <- SvecPSD(Variable(c((n * (n + 1L)) %/% 2L, 1L)), n = n)
  expect_equal(constr_size(sv), (n * (n + 1L)) %/% 2L)
  expect_equal(num_cones(sv), 1L)
  expect_equal(cone_sizes(sv), n)

  skip_if_not_installed("clarabel")
  prob <- Problem(Minimize(matrix_trace(X)), list(X %>>% diag(n)))
  d <- problem_data(prob, solver = "CLARABEL")
  expect_equal(nrow(d[[1L]]$A), (n * (n + 1L)) %/% 2L)
  expect_equal(d[[1L]]$dims@psd, n)
})

## A solver that takes FULL PSD matrices (CVXOPT/cccp `psdc`) must NOT be
## converted: its PSD_TRIANGLE_KIND is NA, so expand_cones leaves PSD alone.
## @cvxpy NONE
test_that("a full-matrix PSD solver gets no svec conversion", {
  skip_if_not_installed("cccp")
  n <- 3L
  X <- Variable(c(n, n), symmetric = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)), list(X %>>% diag(n)))
  d <- problem_data(prob, solver = "CVXOPT")
  expect_equal(nrow(d[[1L]]$A), n * n)
})
