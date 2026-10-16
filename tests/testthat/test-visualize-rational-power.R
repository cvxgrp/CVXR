## Regression tests: visualize() on atoms whose Power exponent is a gmp
## rational (bigq). Since the 1.9.2 sqrt/inv_pos -> SOC parity change these
## atoms build PowerApprox with `p_used` as a bigq; gmp's format() method
## misreads `digits` and returns a length-16 vector for a scalar rational,
## which vectorized the smith/relaxed annotation strings and crashed
## .viz_build_model()'s vapply() ("values must be length 1, but FUN(X[[5]])
## result is length 16"). Shipped broken in CVXR 1.9.2; fixed for 1.9.3.
## Found via the cvxr_tutorial CI (03-dcp.qmd visualize chunk).

## @cvxpy NONE
test_that("visualize: sqrt (bigq exponent 1/2) builds without error", {
  x <- Variable()
  prob <- Problem(Minimize(sqrt(1 + square(x))))
  out <- expect_no_error(visualize(prob))
  txt <- paste(capture.output(print(out)), collapse = "\n")
  ## The exponent must render as the scalar rational, not vectorize.
  expect_match(txt, "1/2", fixed = TRUE)
})

## @cvxpy NONE
test_that("visualize: inv_pos (bigq exponent -1) builds without error", {
  x <- Variable()
  prob <- Problem(Minimize(inv_pos(x)))
  expect_no_error(visualize(prob))
})

## @cvxpy NONE
test_that("visualize: fractional power (bigq exponent 3/2) builds without error", {
  x <- Variable()
  prob <- Problem(Minimize(power(x, 1.5)))
  expect_no_error(visualize(prob))
})

## @cvxpy NONE
test_that("visualize: integer power and fractional p_norm still render", {
  ## Guards for the two neighboring paths that were NOT broken: the p == 2
  ## special case and the Pnorm annotation (whose @p is plain numeric).
  x <- Variable()
  y <- Variable(3)
  expect_no_error(visualize(Problem(Minimize(power(x, 2)))))
  expect_no_error(visualize(Problem(Minimize(p_norm(y, 1.5)))))
})
