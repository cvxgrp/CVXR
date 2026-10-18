## Tests for MOSEK solver integration
## All tests are wrapped in skip_if_not_installed("Rmosek") so they
## pass cleanly when Rmosek is not available.

## @cvxpy NONE
test_that("MOSEK LP", {
  skip_if_not_installed("Rmosek")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, 3.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(3.0, 0.0), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK QP", {
  skip_if_not_installed("Rmosek")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + x[1]),
                  list(x[1] + x[2] == 1, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, 0.875, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0.25, 0.75), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK SOCP", {
  skip_if_not_installed("Rmosek")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(cvxr_norm(x, 2) <= 1))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, -sqrt(2), tolerance = 1e-4)
  expect_equal(value(x), matrix(c(-1/sqrt(2), -1/sqrt(2)), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK ExpCone via entr", {
  skip_if_not_installed("Rmosek")
  x <- Variable(2)
  prob <- Problem(Maximize(entr(x[1]) + entr(x[2])),
                  list(x[1] + x[2] == 1))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, log(2), tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0.5, 0.5), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK PowCone via pnorm", {
  skip_if_not_installed("Rmosek")
  x <- Variable(3)
  prob <- Problem(Minimize(cvxr_norm(x, 1.5)),
                  list(sum(x) == 1, x >= 0))
  val <- psolve(prob, solver = "MOSEK")
  ## Reference: CVXPY Clarabel gives ~0.6934
  expect_equal(val, 0.6934, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), rep(1/3, 3), tolerance = 1e-2)
})

## @cvxpy NONE
test_that("MOSEK SDP via lambda_max", {
  skip_if_not_installed("Rmosek")
  X <- Variable(shape = c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)),
                  list(X[1, 1] == 1, X[2, 2] == 2, X[1, 2] >= 0.5))
  val <- psolve(prob, solver = "MOSEK")
  expected <- 1.5 + sqrt(2) / 2  # (3 + sqrt(2)) / 2
  expect_equal(val, expected, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK infeasible problem", {
  skip_if_not_installed("Rmosek")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 1, x <= -1))
  val <- psolve(prob, solver = "MOSEK")
  expect_true(is.infinite(val) && val > 0)
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy NONE
test_that("MOSEK unbounded problem", {
  skip_if_not_installed("Rmosek")
  x <- Variable()
  prob <- Problem(Minimize(x))
  val <- psolve(prob, solver = "MOSEK")
  expect_true(is.infinite(val) && val < 0)
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate"))
})

## @cvxpy NONE
test_that("MOSEK LP matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x[1] >= 0, x[2] >= 0))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK QP matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + x[1]),
                  list(x[1] + x[2] == 1, x[1] >= 0, x[2] >= 0))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK SOCP matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(cvxr_norm(x, 2) <= 1))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK ExpCone matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Maximize(entr(x[1]) + entr(x[2])),
                  list(x[1] + x[2] == 1))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK SDP matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  X <- Variable(shape = c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)),
                  list(X[1, 1] == 1, X[2, 2] == 2, X[1, 2] >= 0.5))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK log_det problem", {
  skip_if_not_installed("Rmosek")
  X <- Variable(shape = c(2, 2), symmetric = TRUE)
  prob <- Problem(Maximize(log_det(X)),
                  list(X[1, 1] <= 2, X[2, 2] <= 3, X == t(X)))
  val <- psolve(prob, solver = "MOSEK")
  ## log(det(diag(2,3))) = log(6)
  expect_equal(val, log(6), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK norm_nuc problem", {
  skip_if_not_installed("Rmosek")
  X <- Variable(shape = c(2, 2))
  prob <- Problem(Minimize(norm_nuc(X)),
                  list(X[1, 1] == 1, X[2, 2] == 1, X[1, 2] == 0, X[2, 1] == 0))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, 2.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK solver name constant is exported", {
  expect_equal(CVXR::MOSEK_SOLVER, "MOSEK")
})

# ── Simplex, bfs, save_file (claims-audit rows C0803 / C0796) ────────────────
## Expected values below were run on CVXPY 1.9.3 + MOSEK 11 (2026-10-09).

## @cvxpy mosek_conif.py::MOSEK::invert
test_that("MOSEK simplex optimizers return the basic solution", {
  skip_if_not_installed("Rmosek")
  ## Rmosek returns ONLY sol$bas for a simplex optimizer. Until 1.9.2.9146
  ## the interface always read sol$itr and so reported "Solver MOSEK failed"
  ## for every simplex solve; mosek_conif.py:463-476 reads the basic solution
  ## for primal_simplex / dual_simplex. All three parameter carriers, both
  ## the string names and MOSEK's optimizertype codes (dual_simplex = 1,
  ## primal_simplex = 6). CVXPY: optimal, 2.5, x = (1, 1.5), duals (0.5, 0)
  ## and 0.5.
  carriers <- list(
    list(MSK_IPAR_OPTIMIZER = "MSK_OPTIMIZER_PRIMAL_SIMPLEX"),
    list(iparam = list(OPTIMIZER = "MSK_OPTIMIZER_DUAL_SIMPLEX")),
    list(mosek_params = list(MSK_IPAR_OPTIMIZER = 6L)),
    list(iparam = list(OPTIMIZER = 1L)))
  for (opt in carriers) {
    x <- Variable(2)
    c1 <- x >= 1
    c2 <- x[1] + 2 * x[2] >= 4
    prob <- Problem(Minimize(sum(x)), list(c1, c2))
    val <- do.call(psolve, c(list(prob, solver = "MOSEK"), opt))
    lab <- paste(deparse(opt), collapse = "")
    expect_equal(status(prob), "optimal", label = lab)
    expect_equal(val, 2.5, tolerance = 1e-6, label = lab)
    expect_equal(as.numeric(value(x)), c(1, 1.5), tolerance = 1e-6, label = lab)
    expect_equal(as.numeric(dual_value(c1)), c(0.5, 0), tolerance = 1e-6, label = lab)
    expect_equal(as.numeric(dual_value(c2)), 0.5, tolerance = 1e-6, label = lab)
  }
})

## @cvxpy mosek_conif.py::MOSEK::invert
test_that("MOSEK simplex certificates are read from the basic solution", {
  skip_if_not_installed("Rmosek")
  z <- Variable(2)
  prob <- Problem(Minimize(sum(z)), list(z <= 1, sum(z) >= 10))
  val <- psolve(prob, solver = "MOSEK",
                MSK_IPAR_OPTIMIZER = "MSK_OPTIMIZER_PRIMAL_SIMPLEX")
  expect_equal(status(prob), "infeasible")
  expect_equal(val, Inf)

  z <- Variable(2)
  prob <- Problem(Minimize(-sum(z)), list(z >= 0, z[1] == z[2]))
  val <- psolve(prob, solver = "MOSEK",
                MSK_IPAR_OPTIMIZER = "MSK_OPTIMIZER_PRIMAL_SIMPLEX")
  expect_equal(status(prob), "unbounded")
  expect_equal(val, -Inf)
})

## @cvxpy mosek_conif.py::MOSEK::handle_options mosek_conif.py::MOSEK::invert
test_that("MOSEK bfs = TRUE returns a vertex of the optimal face", {
  skip_if_not_installed("Rmosek")
  ## Degenerate LP: every point of the segment x1 + x2 = 1, x >= 0 is
  ## optimal. With bfs = TRUE basis identification runs
  ## (MSK_IPAR_INTPNT_BASIS = always, mosek_conif.py:678-679) and the basic
  ## solution is a vertex. CVXPY: (0, 1) with bfs; without it CVXPY reports
  ## the analytic center (0.5, 0.5) because it hands MOSEK the DUAL of every
  ## continuous problem (mosek_conif.py:198, :533), whereas CVXR hands MOSEK
  ## the primal and its presolve already lands on a vertex, so the
  ## no-bfs iterate is not compared here. Until 1.9.2.9146 `bfs` was refused
  ## as an invalid option (C0796).
  deg <- function() {
    y <- Variable(2)
    list(y = y, prob = Problem(Minimize(sum(y)), list(y >= 0, sum(y) >= 1)))
  }
  d <- deg()
  val <- psolve(d$prob, solver = "MOSEK", bfs = TRUE)
  expect_equal(status(d$prob), "optimal")
  expect_equal(val, 1, tolerance = 1e-6)
  expect_equal(sort(as.numeric(value(d$y))), c(0, 1), tolerance = 1e-6)

  ## With cones there is no basic solution: bfs_active requires
  ## task.getnumcone() == 0 (mosek_conif.py:469), so the interior-point
  ## solution is read. CVXPY: optimal, 0.7071068, x = (0.5, 0.5).
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2)), list(sum(x) >= 1))
  val <- psolve(prob, solver = "MOSEK", bfs = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(val, sqrt(0.5), tolerance = 1e-6)
  expect_equal(as.numeric(value(x)), c(0.5, 0.5), tolerance = 1e-5)
})

## @cvxpy mosek_conif.py::MOSEK::handle_options mosek_conif.py::MOSEK::solve_via_data
test_that("MOSEK save_file writes the task before optimizing", {
  skip_if_not_installed("Rmosek")
  ## mosek_conif.py:248-250 `task.writedata(save_file)`; Rmosek's
  ## `writebefore` is the same write. Was refused as an invalid option.
  f <- tempfile(fileext = ".task")
  on.exit(unlink(f), add = TRUE)
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  val <- psolve(prob, solver = "MOSEK", save_file = f)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 2, tolerance = 1e-6)
  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)
})

## @cvxpy mosek_conif.py::MOSEK::handle_options
test_that("MOSEK basis identification is off unless bfs is requested", {
  skip_if_not_installed("Rmosek")
  ## mosek_conif.py:678-681 sets MSK_IPAR_INTPNT_BASIS to never (always under
  ## bfs) AFTER the user's parameters. MOSEK's own default (IF_FEASIBLE) ran
  ## basis identification on every LP and returned a sol$bas nobody read.
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1, x[1] + 2 * x[2] >= 4))
  pd <- problem_data(prob, solver = "MOSEK")
  solve_raw <- function(opts) {
    solve_via_data(pd$chain, pd$data, warm_start = FALSE, verbose = FALSE,
                   solver_opts = opts)
  }
  raw <- solve_raw(list())
  expect_identical(names(raw$sol), "itr")
  expect_identical(raw[[CVXR:::MOSEK_SOL_TYPE]], "itr")

  raw <- solve_raw(list(bfs = TRUE))
  expect_setequal(names(raw$sol), c("itr", "bas"))
  expect_identical(raw[[CVXR:::MOSEK_SOL_TYPE]], "bas")

  ## A user's own INTPNT_BASIS is overridden, as upstream.
  raw <- solve_raw(list(iparam = list(INTPNT_BASIS = "MSK_BI_ALWAYS")))
  expect_identical(names(raw$sol), "itr")

  raw <- solve_raw(list(MSK_IPAR_OPTIMIZER = "MSK_OPTIMIZER_DUAL_SIMPLEX"))
  expect_identical(names(raw$sol), "bas")
  expect_identical(raw[[CVXR:::MOSEK_SOL_TYPE]], "bas")
})

## Claims-audit decision 28 (CVXR 1.9.2.9150). MOSEK's MSK_RES_TRM_* response
## codes (max iterations 100000, max time 100001, stall 100006, ...) are
## terminations, not errors: Rmosek returns them in response$code with the
## solution attached and its own solsta. CVXPY reads `rescode =
## task.optimize()` (mosek_conif.py:253) and never consults it again, except to
## warn on trm_max_time (:255-258); the status comes from solsta. Until 9149
## CVXR returned SOLVER_ERROR for any nonzero code. Measured against CVXPY 1.9.3
## + MOSEK 11.2.6 (notes/audit/claims_audit_item1.md, Decision 28): the
## iteration-limited solve is a SolverError by default and `optimal_inaccurate`
## under accept_unknown; the time-limited solve warns and then does the same.

## @cvxpy NONE
test_that("MOSEK termination codes are read through to the solution status", {
  skip_if_not_installed("Rmosek")
  set.seed(1)
  n <- 6; m <- 3
  A <- matrix(round(rnorm(m * n), 2), m, n)
  b <- as.numeric(A %*% (runif(n) + 0.5))
  cvec <- runif(n) + 1
  mk <- function() {
    x <- Variable(n)
    list(x = x, prob = Problem(Minimize(sum(cvec * x)), list(A %*% x == b, x >= 0)))
  }
  ## MSK_RES_TRM_MAX_ITERATIONS with solsta UNKNOWN: a solver error, as the
  ## UNKNOWN entry of the status map says, not because the code was nonzero.
  p <- mk()
  ## class: problem.py:1453-1456 raises SolverError at unpack (decision 19, 9151).
  expect_error(psolve(p$prob, solver = "MOSEK", MSK_IPAR_INTPNT_MAX_ITERATIONS = 1L),
               "MOSEK.*failed", class = "SolverError")
  ## The same termination under accept_unknown keeps the iterate, as CVXPY.
  p <- mk()
  expect_warning(
    val <- psolve(p$prob, solver = "MOSEK", MSK_IPAR_INTPNT_MAX_ITERATIONS = 1L,
                  accept_unknown = TRUE),
    "inaccurate")
  expect_identical(status(p$prob), "optimal_inaccurate")
  expect_true(is.finite(val))
  expect_length(as.numeric(value(p$x)), n)
  expect_true(all(is.finite(value(p$x))))
})

## @cvxpy NONE
test_that("MOSEK time limit warns as CVXPY and the terminated solution is kept", {
  skip_if_not_installed("Rmosek")
  x <- Variable(3)
  mk <- function() Problem(Minimize(sum(x)), list(x >= 1, x[1] + 2 * x[2] >= 4))
  ## mosek_conif.py:255-258: the warning fires on MSK_RES_TRM_MAX_TIME whatever
  ## the status; without accept_unknown the UNKNOWN solsta is still an error.
  prob <- mk()
  expect_warning(
    expect_error(psolve(prob, solver = "MOSEK", MSK_DPAR_OPTIMIZER_MAX_TIME = 0),
                 "MOSEK.*failed", class = "SolverError"),
    "terminated by time limit")
  prob <- mk()
  expect_warning(
    expect_warning(
      val <- psolve(prob, solver = "MOSEK", MSK_DPAR_OPTIMIZER_MAX_TIME = 0,
                    accept_unknown = TRUE),
      "inaccurate"),
    "terminated by time limit")
  expect_identical(status(prob), "optimal_inaccurate")
  expect_true(is.finite(val))
})
