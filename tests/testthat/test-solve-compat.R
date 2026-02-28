## Tests for solve() backward-compatibility shim and verbose psolve() output

# ═══════════════════════════════════════════════════════════════════
# Deprecation shim tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("solve() returns cvxr_result S3 class", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  result <- solve(prob)
  expect_s3_class(result, "cvxr_result")
})

## @cvxpy NONE
test_that("result$value matches psolve() return", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  pval <- psolve(prob)
  result <- solve(prob)
  expect_equal(result$value, pval, tolerance = 1e-6)
})

## @cvxpy NONE
test_that("result$status is optimal", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  result <- solve(prob)
  expect_equal(result$status, OPTIMAL)
})

## @cvxpy NONE
test_that("result$solver is a character string", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  result <- solve(prob)
  expect_type(result$solver, "character")
  expect_true(nzchar(result$solver))
})

## @cvxpy NONE
test_that("result$getValue emits deprecation warning", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  result <- solve(prob)
  ## Reset the once-per-session frequency guard for testing
  rlang::reset_warning_verbosity("cvxr_getValue_deprecated")
  expect_warning(val <- result$getValue(x), "deprecated")
  expect_equal(as.numeric(val), 5, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("result$getDualValue emits deprecation warning", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  con <- x >= 5
  prob <- Problem(Minimize(x), list(con))
  result <- solve(prob)
  ## Reset the once-per-session frequency guard for testing
  rlang::reset_warning_verbosity("cvxr_getDualValue_deprecated")
  expect_warning(dv <- result$getDualValue(con), "deprecated")
  expect_true(is.numeric(dv) || is.matrix(dv))
})

## @cvxpy NONE
test_that("print.cvxr_result works without error", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  result <- solve(prob)
  expect_output(print(result), "Status")
})

## @cvxpy NONE
test_that("psolve() still returns numeric scalar", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  val <- psolve(prob)
  expect_true(is.numeric(val))
  expect_equal(val, 5, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("infeasible problem: result$status and result$value", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5, x <= 3))
  result <- solve(prob)
  expect_s3_class(result, "cvxr_result")
  expect_equal(result$status, INFEASIBLE)
  expect_equal(result$value, Inf)
})

# ═══════════════════════════════════════════════════════════════════
# Verbose output tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("verbose=FALSE produces no messages", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  expect_no_message(psolve(prob, verbose = FALSE))
})

## @cvxpy NONE
test_that("verbose=TRUE produces CVXR header message", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  expect_message(psolve(prob, verbose = TRUE), "CVXR v")
})

## @cvxpy NONE
test_that("verbose=TRUE produces Compilation message", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  expect_message(psolve(prob, verbose = TRUE), "Compilation")
})

## @cvxpy NONE
test_that("verbose=TRUE produces Summary message", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  expect_message(psolve(prob, verbose = TRUE), "Summary")
})

## @cvxpy NONE
test_that("compile_time is stored in problem cache", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  psolve(prob, verbose = FALSE)
  ct <- prob@.cache$compile_time
  expect_true(is.numeric(ct))
  expect_true(ct >= 0)
})

## @cvxpy NONE
test_that("solve_time is stored in problem cache", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  psolve(prob, verbose = FALSE)
  st <- prob@.cache$solve_time
  expect_true(is.numeric(st))
  expect_true(st >= 0)
})

# ═══════════════════════════════════════════════════════════════════
# Backward-compatibility alias tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("tv() is deprecated alias for total_variation", {
  rlang::reset_warning_verbosity("cvxr_tv_deprecated")
  x <- Variable(3)
  expect_warning(e <- tv(x), "deprecated")
  expect_true(S7::S7_inherits(e, Expression))
})

## @cvxpy NONE
test_that("norm2() is deprecated alias for p_norm", {
  rlang::reset_warning_verbosity("cvxr_norm2_deprecated")
  x <- Variable(3)
  expect_warning(e <- norm2(x), "deprecated")
  expect_true(S7::S7_inherits(e, Pnorm))
})

## @cvxpy NONE
test_that("norm2() solves correctly (with deprecation warning)", {
  skip_if_not_installed("clarabel")
  rlang::reset_warning_verbosity("cvxr_norm2_deprecated")
  x <- Variable(3)
  expect_warning(obj <- norm2(x), "deprecated")
  prob <- Problem(Minimize(obj), list(x >= c(3, 4, 0)))
  val <- psolve(prob)
  expect_equal(val, 5, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("installed_solvers() returns character vector", {
  solvers <- installed_solvers()
  expect_type(solvers, "character")
  expect_true(length(solvers) >= 1)
  ## At least clarabel should be installed for tests to run
  expect_true("CLARABEL" %in% solvers)
})

# ═══════════════════════════════════════════════════════════════════
# Axis convention migration tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("axis=0 gives informative migration error", {
  x <- Variable(c(3, 2))
  expect_error(sum_entries(x, axis = 0L), "1-based axis indexing")
})

## @cvxpy NONE
test_that("axis=1 and axis=2 work correctly", {
  x <- Variable(c(3, 2))
  s1 <- SumEntries(x, axis = 1L)
  s2 <- SumEntries(x, axis = 2L)
  expect_equal(s1@shape, c(3L, 1L))
  expect_equal(s2@shape, c(1L, 2L))
})
