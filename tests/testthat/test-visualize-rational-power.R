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

## smith_annotation(Constant) formats a 1x1 value; for a sparse constant that
## is a dgCMatrix, and base format() prints its S4 slot summary unless
## Matrix's format generic is imported (1.9.2.9510, sweep item 4 pass I).
## @cvxpy NONE
test_that("visualize: a 1x1 sparse constant renders as its number", {
  x <- Variable(1, name = "x")
  dense  <- Problem(Minimize(square(Constant(2.5) * x)))
  sparse <- Problem(Minimize(square(Constant(Matrix::Matrix(2.5, 1, 1, sparse = TRUE)) * x)))
  txt_dense  <- paste(capture.output(print(visualize(dense))), collapse = "\n")
  txt_sparse <- paste(capture.output(print(visualize(sparse))), collapse = "\n")
  expect_false(grepl("S4 class", txt_sparse, fixed = TRUE))
  expect_match(txt_sparse, "2.5", fixed = TRUE)
  ## format() of a 1x1 value is a 1x1 character matrix on both sides (Matrix's
  ## adds dimnames); the annotation strings are consumed through
  ## sprintf()/paste(), so compare the text.
  ann <- CVXR:::smith_annotation
  lbl <- ann(Constant(Matrix::Matrix(2.5, 1, 1, sparse = TRUE)))$latex_name
  expect_type(lbl, "character")
  expect_length(lbl, 1L)
  expect_identical(as.character(lbl), as.character(ann(Constant(2.5))$latex_name))
})
