## Tests for DGP + DPP integration (Phase DGP-4)
## CVXPY reference values verified via `uv run python`

## @cvxpy NONE
test_that("parametric GP uses DPP fast path on re-solve", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  a <- Parameter(pos = TRUE)
  b <- Parameter(pos = TRUE)
  value(a) <- 2.0
  value(b) <- 3.0

  prob <- Problem(Minimize(x + y), list(a * x * b * y >= 1))

  ## Solve 1: a=2, b=3 → 6xy >= 1, min x+y = 2/sqrt(6) ≈ 0.8165
  val1 <- psolve(prob, gp = TRUE)
  expect_equal(val1, 0.816497, tolerance = 1e-3)

  ## Solve 2: change params → a=4, b=1 → 4xy >= 1, min x+y = 1.0
  value(a) <- 4.0
  value(b) <- 1.0
  val2 <- psolve(prob, gp = TRUE)
  expect_equal(val2, 1.0, tolerance = 1e-3)

  ## Verify DPP fast path was used (param_prog cached)
  expect_false(is.null(prob@.cache$param_prog))

  ## Solve 3: a=1, b=1 → xy >= 1, min x+y = 2.0
  value(a) <- 1.0
  value(b) <- 1.0
  val3 <- psolve(prob, gp = TRUE)
  expect_equal(val3, 2.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("DGP+DPP: is_dgp_dpp detects parametric GP", {
  x <- Variable(pos = TRUE)
  p <- Parameter(pos = TRUE)
  value(p) <- 1.0

  ## Parametric GP: DGP-DPP compliant
  prob <- Problem(Minimize(x), list(x >= p))
  ## DGP DPP should be TRUE (parameter in positive position)
  expect_true(CVXR:::.is_dgp_dpp(prob))

  ## But standard is_dpp is FALSE (not DCP because of param*var product)
  prob2 <- Problem(Minimize(x), list(p * x >= 1))
  expect_true(CVXR:::.is_dgp_dpp(prob2))
})

## @cvxpy NONE
test_that("DGP+DPP: no EvalParams in chain for DPP-compliant DGP", {
  x <- Variable(pos = TRUE)
  p <- Parameter(pos = TRUE)
  value(p) <- 2.0

  prob <- Problem(Minimize(x), list(p * x >= 1))
  chain <- CVXR:::construct_solving_chain(prob, gp = TRUE)

  ## Chain should NOT contain EvalParams (DPP fast path enabled)
  has_eval_params <- any(vapply(chain@reductions, function(r) {
    S7_inherits(r, EvalParams)
  }, logical(1L)))
  expect_false(has_eval_params)

  ## Chain SHOULD contain Dgp2Dcp
  has_dgp2dcp <- any(vapply(chain@reductions, function(r) {
    S7_inherits(r, Dgp2Dcp)
  }, logical(1L)))
  expect_true(has_dgp2dcp)
})

## @cvxpy NONE
test_that("DGP+DPP: parametric power constraint", {
  x <- Variable(pos = TRUE)
  p <- Parameter(pos = TRUE)
  value(p) <- 2.0

  ## Minimize x^p subject to x >= 3
  prob <- Problem(Minimize(power(x, p)), list(x >= 3))
  val1 <- psolve(prob, gp = TRUE)
  ## x=3, x^2 = 9
  expect_equal(val1, 9.0, tolerance = 1e-3)

  ## Re-solve with p=3 (DPP fast path)
  value(p) <- 3.0
  val2 <- psolve(prob, gp = TRUE)
  ## x=3, x^3 = 27
  expect_equal(val2, 27.0, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("DGP+DPP: verbose confirms fast path", {
  x <- Variable(pos = TRUE)
  p <- Parameter(pos = TRUE)
  value(p) <- 1.0

  prob <- Problem(Minimize(x + p / x))

  ## First solve: cold path
  val1 <- psolve(prob, gp = TRUE)
  expect_equal(val1, 2.0, tolerance = 1e-3)

  ## Second solve: fast path
  value(p) <- 4.0
  output <- capture.output(val2 <- psolve(prob, gp = TRUE, verbose = TRUE),
                           type = "message")
  expect_equal(val2, 4.0, tolerance = 1e-3)
  expect_true(any(grepl("fast path", output, ignore.case = TRUE)))
})
