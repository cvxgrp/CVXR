## CVXPY 1.9 DNLP affine/structural-atom parity (nlp_tests/test_affine_atoms.py).
##
## sin/cos composed with affine structural atoms (reshape, vstack, upper_tri,
## diag, division, conv/convolve, 1-D matmul) are smooth-over-affine, hence
## DNLP. CVXPY checks only that the solve reaches an optimum and runs a Python
## DerivativeChecker (no R equivalent, dropped). Variables use box bounds
## [-1, 1] (CVXR `bounds = list(-1, 1)`). CVXPY's cp.nlp.sin/cos map to CVXR's
## sin/cos.

.aff_nlp_optimal <- function(prob) {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
  output <- capture.output(psolve(prob, nlp = TRUE), type = "output")
  expect_equal(status(prob), "optimal")
}

## @cvxpy nlp_tests/test_affine_atoms.py::TestAffineAtoms::test_reshape_C_order
test_that("DNLP affine: reshape(sin(X)) in C order solves", {
  set.seed(0)
  m <- 3; n <- 4
  X <- Variable(c(m, n), bounds = list(-1, 1))
  value(X) <- matrix(runif(m * n), m, n)
  obj <- sum_entries(reshape_expr(sin(X), c(m * n, 1L), order = "C"))
  .aff_nlp_optimal(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/test_affine_atoms.py::TestAffineAtoms::test_vstack
test_that("DNLP affine: vstack(sin(x), cos(x)) solves", {
  set.seed(0)
  n <- 5
  x <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n)
  obj <- sum_entries(vstack(sin(x), cos(x)))
  .aff_nlp_optimal(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/test_affine_atoms.py::TestAffineAtoms::test_upper_tri
test_that("DNLP affine: upper_tri(sin(X)) solves", {
  set.seed(0)
  n <- 4
  X <- Variable(c(n, n), bounds = list(-1, 1))
  value(X) <- matrix(runif(n * n), n, n)
  obj <- sum_entries(upper_tri(sin(X)))
  .aff_nlp_optimal(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/test_affine_atoms.py::TestAffineAtoms::test_diag_mat
test_that("DNLP affine: diag(sin(X)) solves", {
  set.seed(0)
  n <- 4
  X <- Variable(c(n, n), bounds = list(-1, 1))
  value(X) <- matrix(runif(n * n), n, n)
  obj <- sum_entries(diag(sin(X)))
  .aff_nlp_optimal(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/test_affine_atoms.py::TestAffineAtoms::test_div_by_scalar
test_that("DNLP affine: sin(x) / scalar solves", {
  set.seed(0)
  n <- 5
  x <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n)
  .aff_nlp_optimal(Problem(Minimize(sum_entries(sin(x) / 3.0))))
})

## @cvxpy nlp_tests/test_affine_atoms.py::TestAffineAtoms::test_div_by_vector
test_that("DNLP affine: sin(x) / vector solves", {
  set.seed(0)
  n <- 5
  c_vec <- as.numeric(1:n)
  x <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n)
  .aff_nlp_optimal(Problem(Minimize(sum_entries(sin(x) / c_vec))))
})

## @cvxpy nlp_tests/test_affine_atoms.py::TestAffineAtoms::test_conv_constant_kernel
test_that("DNLP affine: conv(kernel, sin(x)) solves", {
  set.seed(0)
  n <- 6
  kernel <- c(0.5, -1.0, 0.25)
  x <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n)
  .aff_nlp_optimal(Problem(Minimize(sum_entries(conv(kernel, sin(x))))))
})

## @cvxpy nlp_tests/test_affine_atoms.py::TestAffineAtoms::test_convolve_alias
test_that("DNLP affine: convolve(kernel, sin(x)) solves (conv alias)", {
  set.seed(0)
  n <- 6
  kernel <- c(0.3, 0.6, -0.2)
  x <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n)
  .aff_nlp_optimal(Problem(Minimize(sum_entries(convolve(kernel, sin(x))))))
})

## @cvxpy nlp_tests/test_affine_atoms.py::TestAffineAtoms::test_matmul_1d_dot_product
test_that("DNLP affine: sin(x)^T sin(x) (1-D dot product) solves", {
  set.seed(0)
  n <- 5
  x <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n)
  dot <- t(sin(x)) %*% sin(x)
  expect_equal(dim(dot), c(1L, 1L))
  .aff_nlp_optimal(Problem(Minimize(dot)))
})
