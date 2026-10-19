## Claims audit decision 24b (notes/audit/claims_audit_item1.md): Gurobi and HiGHS
## interface options, CVXR 1.9.2.9159. Gurobi: `init_value` MIP start from the
## variables' current values (C0772/C0859), `save_file` (C0857), `env` abort;
## HiGHS: `highs_options` unpack and kError on a rejected option (C0780/C0865),
## the 1.9.3 column-name pattern (C0774). Expected behavior was measured on
## cvxpy 1.9.3 in a pinned `uv run` process (scratch d24b_cvxpy.py, d24b_ws.py);
## every test here was seen failing on 1.9.2.9158 first.

# -- helpers --------------------------------------------------------------------

## A 60-item knapsack with a greedy feasible start; big enough that Gurobi's
## log distinguishes "Loaded user MIP start" from its own heuristics.
d24b_knapsack <- function(n = 60L) {
  set.seed(7)
  w <- sample(10:40, n, replace = TRUE)
  v <- sample(10:50, n, replace = TRUE)
  cap <- floor(sum(w) / 3)
  o <- order(v / w, decreasing = TRUE)
  z0 <- numeric(n)
  tot <- 0
  for (i in o) {
    if (tot + w[i] <= cap) {
      z0[i] <- 1
      tot <- tot + w[i]
    }
  }
  list(w = w, v = v, cap = cap, z0 = z0)
}

## Solve with Gurobi logging to a file (OutputFlag = 0 would silence the file
## too, so verbose = TRUE with LogToConsole = 0) and return the MIP-start lines.
d24b_mip_start_lines <- function(prob, ...) {
  lf <- tempfile(fileext = ".log")
  on.exit(unlink(lf), add = TRUE)
  suppressMessages(psolve(prob, solver = "GUROBI", verbose = TRUE,
                          LogToConsole = 0, LogFile = lf, ...))
  grep("MIP start", readLines(lf), value = TRUE)
}

d24b_loaded_start_objective <- function(hits) {
  loaded <- grep("Loaded user MIP start with objective", hits, value = TRUE)
  if (length(loaded) == 0L) return(NA_real_)
  as.numeric(sub(".*objective (-?[0-9.eE+-]+).*", "\\1", loaded[1L]))
}

# -- C0772: Gurobi conic init_value warm start -----------------------------------

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_warm_start
test_that("Gurobi conic: warm_start = TRUE seeds the MIP start from the variables' current values on a cold solve", {
  require_solver("GUROBI")
  k <- d24b_knapsack()
  z <- Variable(length(k$w), boolean = TRUE)
  value(z) <- k$z0
  prob <- Problem(Maximize(sum(k$v * z)), list(sum(k$w * z) <= k$cap))
  hits <- d24b_mip_start_lines(prob, warm_start = TRUE)
  ## gurobi_conif.py:211-214: `elif warm_start: start_vals = data['init_value'][:n]`.
  ## Maximize is minimized negated, so the loaded start objective is -v'z0;
  ## matching it proves the start landed on the right columns.
  expect_equal(d24b_loaded_start_objective(hits), -sum(k$v * k$z0))
  expect_equal(status(prob), "optimal")
  expect_true(sum(k$v * value(z)) >= sum(k$v * k$z0))

  ## Negative control: the branch is gated on warm_start, as upstream.
  z2 <- Variable(length(k$w), boolean = TRUE)
  value(z2) <- k$z0
  prob2 <- Problem(Maximize(sum(k$v * z2)), list(sum(k$w * z2) <= k$cap))
  expect_length(d24b_mip_start_lines(prob2, warm_start = FALSE), 0L)
})

# -- C0859: Gurobi QP init_value warm start ---------------------------------------

## @cvxpy test_qp_solvers.py::TestQp::test_gurobi_warmstart
test_that("Gurobi QP: warm_start = TRUE seeds the MIP start from the variables' current values on a cold solve", {
  require_solver("GUROBI")
  k <- d24b_knapsack()
  z <- Variable(length(k$w), boolean = TRUE)
  value(z) <- k$z0
  ## Quadratic objective without a constant term, so Gurobi's reported start
  ## objective is exactly z0'z0 - v'z0 / 50 (gurobi_qpif.py:194-197).
  prob <- Problem(Minimize(sum_squares(z) - sum(k$v * z) / 50),
                  list(sum(k$w * z) <= k$cap))
  hits <- d24b_mip_start_lines(prob, warm_start = TRUE)
  expect_equal(d24b_loaded_start_objective(hits),
               sum(k$z0^2) - sum(k$v * k$z0) / 50, tolerance = 1e-6)
  expect_equal(status(prob), "optimal")

  z2 <- Variable(length(k$w), boolean = TRUE)
  value(z2) <- k$z0
  prob2 <- Problem(Minimize(sum_squares(z2) - sum(k$v * z2) / 50),
                   list(sum(k$w * z2) <= k$cap))
  expect_length(d24b_mip_start_lines(prob2, warm_start = FALSE), 0L)
})

## @cvxpy NONE
test_that("stack_vals: column-major values, default for unset variables", {
  a <- Variable(c(2, 2))
  b <- Variable(3)
  value(a) <- matrix(1:4, 2, 2)
  ## utilities.py:24-43 ravels in Fortran order and fills unset variables.
  expect_equal(CVXR:::stack_vals(list(a, b), NA_real_),
               c(1, 2, 3, 4, NA, NA, NA))
  expect_equal(CVXR:::stack_vals(list(b, a), 0), c(0, 0, 0, 1, 2, 3, 4))
  expect_equal(CVXR:::stack_vals(list(), NA_real_), numeric(0))
})

# -- C0857: Gurobi save_file and env ---------------------------------------------

## @cvxpy NONE
test_that("Gurobi: save_file writes the model before solving, on both interfaces", {
  require_solver("GUROBI")
  x <- Variable(2)
  lp <- Problem(Minimize(sum(x)), list(x >= 1, sum(x) <= 5))
  f1 <- tempfile(fileext = ".lp")
  f2 <- tempfile(fileext = ".lp")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  ## gurobi_conif.py:235-236 `model.write(solver_opts['save_file'])`; the key
  ## is an INTERFACE_ARGS entry (:43) and must not reach Gurobi as a parameter.
  expect_equal(psolve(lp, solver = "GUROBI", save_file = f1), 2, tolerance = 1e-6)
  expect_true(file.exists(f1))
  expect_true(any(grepl("Minimize", readLines(f1))))

  y <- Variable(2)
  qp <- Problem(Minimize(sum_squares(y - c(1, 2))), list(y >= 0))
  expect_equal(psolve(qp, solver = "GUROBI", save_file = f2), 0, tolerance = 1e-6)
  expect_true(file.exists(f2))
})

## @cvxpy NONE
test_that("Gurobi: env is rejected with a message, since the R gurobi API has no environment object", {
  require_solver("GUROBI")
  x <- Variable(2)
  lp <- Problem(Minimize(sum(x)), list(x >= 1))
  expect_error(psolve(lp, solver = "GUROBI", env = list()), "environment")
  y <- Variable(2)
  qp <- Problem(Minimize(sum_squares(y)), list(y >= 1))
  expect_error(psolve(qp, solver = "GUROBI", env = list()), "environment")
})

# -- C0780 / C0865: HiGHS options --------------------------------------------------

## test_conic_solvers.py:3005-3036 runs on test_lp_0 and test_mi_lp_0; in CVXR
## the LP takes the QP interface (highs_qpif.R) and the MILP the conic one
## (highs_conif.R), so the pair covers both files.
d24b_highs_problems <- function() list(sth_lp_0(), sth_mi_lp_0())

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_options
test_that("HiGHS: an option HiGHS rejects aborts with kError instead of being dropped", {
  skip_if_not_installed("highs")
  for (h in d24b_highs_problems()) {
    ## highs_conif.py:313-320 / highs_qpif.py:212-219.
    expect_error(psolve(h$prob, solver = "HIGHS",
                        highs_options = list(invalid_highs_option = NULL)),
                 "kError")
    expect_error(psolve(h$prob, solver = "HIGHS", invalid_highs_option = NULL),
                 "kError")
    expect_error(psolve(h$prob, solver = "HIGHS",
                        highs_options = list(threads = "_invalid_value_")),
                 "kError")
    expect_error(psolve(h$prob, solver = "HIGHS", threads = "_invalid_value_"),
                 "kError")
  }
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_options
test_that("HiGHS: an option given both directly and inside highs_options aborts", {
  skip_if_not_installed("highs")
  for (h in d24b_highs_problems()) {
    ## highs_conif.py:103: `dict(**solver_opts, **highs_options)` raises
    ## TypeError on a duplicate keyword.
    expect_error(psolve(h$prob, solver = "HIGHS", presolve = "off",
                        highs_options = list(presolve = "off")),
                 "presolve")
  }
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_options
test_that("HiGHS: highs_options carries every option type, including the name-clashing `solver`", {
  skip_if_not_installed("highs")
  ## The dict at test_conic_solvers.py:3024-3035, verbatim.
  highs_options <- list(
    solver = "simplex",
    output_flag = TRUE,
    time_limit = 0.123,
    random_seed = 1234,
    mip_rel_gap = 1.0,
    primal_feasibility_tolerance = 1e-3,
    dual_feasibility_tolerance = 1e-3,
    presolve_rule_off = 2^10
  )
  for (h in d24b_highs_problems()) {
    val <- psolve(h$prob, solver = "HIGHS", highs_options = highs_options)
    expect_equal(val, h$expect_obj, tolerance = 1e-4)
    expect_equal(as.numeric(value(variables(h$prob)[[1L]])), h$expect_x, tolerance = 1e-4)
  }
})

# -- C0774: HiGHS column-name pattern ---------------------------------------------

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_validate_column_name
test_that("HiGHS: validate_column_name follows the 1.9.3 pattern, which has no brackets in the tail", {
  ## highs_conif.py:33-35: the tail class ends `|~.=()<>]{,254}`; 1.8.2 had `[]` too.
  expect_error(validate_column_name("x[1]"), "Invalid column name")
  expect_error(validate_column_name("x]"), "Invalid column name")
  expect_silent(validate_column_name("x(1)"))
})
