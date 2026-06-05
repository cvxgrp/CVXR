## Tests for get_bounds() — the CVXPY 1.9.0 bounds API used by the DNLP solve
## path. As of Wave 2 (#3080) get_bounds returns bounds SHAPED to the expression
## dimensions (matrices matching shape, column-major), matching CVXPY -- not
## flat vectors. Sparse #3103 and symbolic/parametric bounds are active here;
## get_bounds() skips symbolic bounds, while CvxAttr2Constr enforces them when
## bound attributes are lowered.

# ── get_bounds: unbounded ────────────────────────────────────────────
## @cvxpy test_bounds.py::TestBoundsPropagation::test_variable_attribute_bounds
test_that("get_bounds returns +-Inf for an unbounded variable", {
  x <- Variable(3)                       # shape (3,1) -> 3x1 matrix
  b <- get_bounds(x)
  expect_equal(b[[1]], matrix(-Inf, 3, 1))
  expect_equal(b[[2]], matrix(Inf, 3, 1))
})

# ── get_bounds: sign attributes ──────────────────────────────────────
## @cvxpy test_bounds.py::TestBoundsPropagation::test_variable_attribute_bounds
test_that("get_bounds folds in nonneg / nonpos / boolean", {
  xn <- Variable(2, nonneg = TRUE)
  expect_equal(get_bounds(xn)[[1]], matrix(c(0, 0), 2, 1))
  expect_equal(get_bounds(xn)[[2]], matrix(c(Inf, Inf), 2, 1))

  xp <- Variable(2, nonpos = TRUE)
  expect_equal(get_bounds(xp)[[1]], matrix(c(-Inf, -Inf), 2, 1))
  expect_equal(get_bounds(xp)[[2]], matrix(c(0, 0), 2, 1))

  xb <- Variable(2, boolean = TRUE)
  expect_equal(get_bounds(xb)[[1]], matrix(c(0, 0), 2, 1))
  expect_equal(get_bounds(xb)[[2]], matrix(c(1, 1), 2, 1))
})

# ── get_bounds: explicit bounds attribute ────────────────────────────
## @cvxpy test_bounds.py::TestHIGHSBoundsPropagation::test_get_problem_data_array_bounds
test_that("get_bounds returns explicit vector bounds", {
  x <- Variable(2, bounds = list(c(-1, -2), c(3, 4)))
  expect_equal(get_bounds(x)[[1]], matrix(c(-1, -2), 2, 1))
  expect_equal(get_bounds(x)[[2]], matrix(c(3, 4), 2, 1))
})

## @cvxpy test_bounds.py::TestBoundsPropagation::test_variable_with_bounds
test_that("get_bounds broadcasts scalar bounds across all entries", {
  x <- Variable(3, bounds = list(0, 10))
  expect_equal(get_bounds(x)[[1]], matrix(0, 3, 1))
  expect_equal(get_bounds(x)[[2]], matrix(10, 3, 1))
})

## @cvxpy test_bounds.py::TestExtractBounds::test_attribute_with_explicit_bounds
test_that("get_bounds combines explicit bounds with sign attributes", {
  # nonneg tightens the lower bound from -1 up to 0; upper stays 5
  x <- Variable(2, bounds = list(-1, 5), nonneg = TRUE)
  expect_equal(get_bounds(x)[[1]], matrix(c(0, 0), 2, 1))
  expect_equal(get_bounds(x)[[2]], matrix(c(5, 5), 2, 1))
})

## @cvxpy test_bounds.py::TestBoundsPropagation::test_matrix_variable_bounds
test_that("get_bounds is shaped (column-major) for matrix variables", {
  x <- Variable(c(2, 2), bounds = list(c(1, 2, 3, 4), c(5, 6, 7, 8)))
  expect_equal(get_bounds(x)[[1]], matrix(c(1, 2, 3, 4), 2, 2))   # column-major
  expect_equal(get_bounds(x)[[2]], matrix(c(5, 6, 7, 8), 2, 2))
})

# ── Deferred to Wave 2 (written now, dormant) ────────────────────────
## @cvxpy test_bounds.py::TestSparseBounds::test_sparse_variable_with_matching_sparse_bounds
test_that("get_bounds preserves sparse bounds", {
  lb <- Matrix::sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(-1, -2), dims = c(2, 2))
  ub <- Matrix::sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(3, 4), dims = c(2, 2))
  x <- Variable(c(2, 2), bounds = list(lb, ub))
  b <- get_bounds(x)
  expect_s4_class(b[[1]], "sparseMatrix")
  expect_s4_class(b[[2]], "sparseMatrix")
  expect_equal(as.matrix(b[[1]]), as.matrix(lb))
  expect_equal(as.matrix(b[[2]]), as.matrix(ub))
})

## @cvxpy test_parametric_bounds.py::TestParametricBoundsCreation::test_mixed_numeric_and_param test_parametric_bounds.py::TestParametricBoundsCreation::test_parameters_discovered
test_that("get_bounds skips symbolic Parameter bounds (enforced at solve time)", {
  lo <- Parameter(2)
  x <- Variable(2, bounds = list(lo, 5))
  b <- get_bounds(x)
  expect_equal(b[[1]], matrix(c(-Inf, -Inf), 2, 1))
  expect_equal(b[[2]], matrix(c(5, 5), 2, 1))

  prob <- Problem(Minimize(sum_entries(x)), list())
  reduced <- reduction_apply(CvxAttr2Constr(), prob)
  new_constraints <- reduced[[1]]@constraints
  expect_length(new_constraints, 2)
  expect_true(any(vapply(new_constraints, function(con) S7::S7_inherits(con@args[[1L]], Parameter), logical(1))))
})

## @cvxpy test_parametric_bounds.py::TestParametricBoundsCreation::test_basic_creation test_parametric_bounds.py::TestParametricBoundsCreation::test_invalid_shape_raises test_parametric_bounds.py::TestParametricBoundsCreation::test_parameter_rejects_expression_bounds test_parametric_bounds.py::TestParametricBoundsCreation::test_variable_rejects_variable_dependent_bounds test_parametric_bounds.py::TestParametricBoundsCreation::test_parameters_discovered
test_that("parametric bounds validate shape, reject variable-dependent bounds, and expose parameters", {
  lb <- Parameter(name = "lb")
  ub <- Parameter(name = "ub")
  x <- Variable(bounds = list(lb, ub))
  expect_true(S7::S7_inherits(x@attributes$bounds[[1L]], Expression))
  expect_true(S7::S7_inherits(x@attributes$bounds[[2L]], Expression))

  expect_error(Variable(c(2, 2), bounds = list(Parameter(3), 10)),
               "Expression bounds must be scalar")
  y <- Variable()
  expect_error(Variable(bounds = list(y, 10)),
               "must not depend on Variables")
  expect_error(Parameter(bounds = list(Parameter(), 10)),
               "Parametric bounds")

  prob <- Problem(Minimize(x), list())
  ids <- vapply(parameters(prob), function(p) p@id, integer(1))
  expect_true(lb@id %in% ids)
  expect_true(ub@id %in% ids)
})

## @cvxpy test_parametric_bounds.py::TestParametricBoundsSolving::test_re_solve
test_that("parametric bounds solve and re-solve through lowered constraints", {
  lb <- Parameter(value = 2)
  ub <- Parameter(value = 5)
  x <- Variable(bounds = list(lb, ub))
  prob <- Problem(Minimize(x), list())

  expect_equal(psolve(prob, solver = "CLARABEL"), 2, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 2, tolerance = 1e-5)

  value(lb) <- 3
  expect_equal(psolve(prob, solver = "CLARABEL"), 3, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 3, tolerance = 1e-5)
})

## @cvxpy matrix_stuffing.py:48-52,84-88
test_that("HiGHS keeps sparse bounds on lowered constraints", {
  skip_if_not_installed("highs")

  lb <- Matrix::sparseMatrix(i = 1, j = 1, x = -1, dims = c(2, 1))
  sx <- Variable(2, bounds = list(lb, 5))
  sprob <- Problem(Minimize(sum_entries(sx)), list())
  sdata <- problem_data(sprob, HIGHS_SOLVER)$data
  expect_null(sdata[[LOWER_BOUNDS]])
  expect_null(sdata[[UPPER_BOUNDS]])
  expect_true(sdata[[SD_DIMS]]@nonneg >= 2L)
})

## @cvxpy matrix_stuffing.py::extract_bounds_tensor test_parametric_bounds.py::TestParametricBoundsSolving::test_bounded_solver_basic test_parametric_bounds.py::TestParametricBoundsSolving::test_bounded_solver_re_solve
test_that("HiGHS receives parametric bounds through 2D bound tensors", {
  skip_if_not_installed("highs")

  lo <- Parameter(value = 2)
  hi <- Parameter(value = 5)
  px <- Variable(bounds = list(lo, hi))
  pprob <- Problem(Minimize(px), list())
  pdata <- problem_data(pprob, HIGHS_SOLVER)$data
  expect_equal(pdata[[LOWER_BOUNDS]], 2)
  expect_equal(pdata[[UPPER_BOUNDS]], 5)
  expect_equal(pdata[[SD_DIMS]]@nonneg, 0L)

  expect_equal(psolve(pprob, solver = HIGHS_SOLVER), 2, tolerance = 1e-5)
  expect_equal(as.numeric(value(px)), 2, tolerance = 1e-5)

  value(lo) <- 3
  expect_equal(psolve(pprob, solver = HIGHS_SOLVER), 3, tolerance = 1e-5)
  expect_equal(as.numeric(value(px)), 3, tolerance = 1e-5)
  expect_false(is.null(pprob@.cache$param_prog@lb_tensor))
  expect_false(is.null(pprob@.cache$param_prog@ub_tensor))
})

## @cvxpy test_parametric_bounds.py::TestParametricBoundsSolving::test_matrix_variable
test_that("HiGHS parametric bound tensors support exact 2D shapes", {
  skip_if_not_installed("highs")

  lo <- Parameter(c(2, 2), value = matrix(c(1, 2, 3, 4), 2, 2))
  X <- Variable(c(2, 2), bounds = list(lo, 10))
  prob <- Problem(Minimize(sum_entries(X)), list())
  data <- problem_data(prob, HIGHS_SOLVER)$data
  expect_equal(data[[LOWER_BOUNDS]], c(1, 2, 3, 4))
  expect_equal(data[[UPPER_BOUNDS]], rep(10, 4))

  expect_equal(psolve(prob, solver = HIGHS_SOLVER), 10, tolerance = 1e-5)
  expect_equal(as.numeric(value(X)), c(1, 2, 3, 4), tolerance = 1e-5)
})

## @cvxpy leaf.py:607-648 (#3087)
test_that("Parameter values may contain matching infinities", {
  p <- Parameter(2)
  value(p) <- c(-Inf, Inf)
  expect_equal(as.numeric(value(p)), c(-Inf, Inf))

  pn <- Parameter(2, nonneg = TRUE)
  value(pn) <- c(0, Inf)
  expect_equal(as.numeric(value(pn)), c(0, Inf))
})

## @cvxpy test_bounds.py::TestBoundsPropagation::test_negation test_bounds.py::TestBoundsPropagation::test_addition test_bounds.py::TestBoundsPropagation::test_scalar_multiplication test_bounds.py::TestBoundsPropagation::test_sum test_bounds.py::TestBoundsPropagation::test_indexing_1d test_bounds.py::TestBoundsPropagation::test_transpose test_bounds.py::TestBoundsPropagation::test_reshape test_bounds.py::TestBoundsPropagation::test_matmul_constant_matrix
## #3080 affine bounds propagation. Every expected value below was computed from
## CVXPY 1.9.0's expr.get_bounds() on the same bounded variables.
test_that("get_bounds propagates through affine expression trees", {
  gb <- function(e) lapply(get_bounds(e), as.numeric)
  x <- Variable(2, bounds = list(1, 3))     # x in [1,3]^2
  y <- Variable(2, bounds = list(-2, 4))

  expect_equal(gb(-x),         list(c(-3, -3), c(-1, -1)))   # neg
  expect_equal(gb(x + 1),      list(c(2, 2),   c(4, 4)))     # add constant
  expect_equal(gb(x + y),      list(c(-1, -1), c(7, 7)))     # add intervals
  expect_equal(gb(2 * x),      list(c(2, 2),   c(6, 6)))     # scale
  expect_equal(gb(sum_entries(x)), list(2, 6))               # sum reduction
  expect_equal(gb(x[1]),       list(1, 3))                   # index

  M <- Variable(c(2, 2), bounds = list(1, 2))
  expect_equal(gb(t(M)), list(rep(1, 4), rep(2, 4)))         # transpose
  expect_equal(gb(reshape_expr(M, c(4L, 1L), order = "C")),
               list(rep(1, 4), rep(2, 4)))                    # reshape (C-order)

  ## matmul by a constant matrix uses the exact split formula.
  A <- matrix(c(1, 3, -2, 4), 2, 2)         # [[1,-2],[3,4]]
  expect_equal(gb(A %*% x), list(c(-5, 7), c(1, 21)))
})

## @cvxpy test_bounds.py::TestBoundsPropagation::test_abs_bounds test_bounds.py::TestBoundsPropagation::test_monotonic_atom_bounds test_bounds.py::TestBoundsPropagation::test_power test_bounds.py::TestBoundsPropagation::test_max_reduction test_bounds.py::TestBoundsPropagation::test_min_reduction test_bounds.py::TestBoundsPropagation::test_norm1 test_bounds.py::TestBoundsPropagation::test_norm_inf
## #3080 nonlinear / PWL bounds propagation. Expected values from CVXPY 1.9.0.
test_that("get_bounds propagates through nonlinear / PWL atoms", {
  gb <- function(e) lapply(get_bounds(e), as.numeric)
  x <- Variable(3, bounds = list(-2, 3))    # x in [-2,3]^3
  y <- Variable(3, bounds = list(1, 4))

  expect_equal(gb(abs(x)),       list(rep(0, 3), rep(3, 3)))
  expect_equal(gb(exp(x)),       list(rep(exp(-2), 3), rep(exp(3), 3)))
  expect_equal(gb(log(y)),       list(rep(0, 3), rep(log(4), 3)))
  expect_equal(gb(power(x, 2)),  list(rep(0, 3), rep(9, 3)))     # even power spans 0
  expect_equal(gb(max_entries(x)), list(-2, 3))
  expect_equal(gb(min_entries(x)), list(-2, 3))
  expect_equal(gb(p_norm(x, 1)),   list(0, 9))                   # sum|x|, |x| in [0,..]
  expect_equal(gb(p_norm(x, Inf)), list(0, 3))                   # max|x|
})

## @cvxpy test_bounds.py::TestBoundsPropagation::test_constant_bounds test_bounds.py::TestBoundsPropagation::test_parameter_with_bounds test_bounds.py::TestBoundsPropagation::test_sign_refinement test_bounds.py::TestBoundsPropagation::test_elementwise_multiplication test_bounds.py::TestBoundsPropagation::test_multiply_spanning_zero test_bounds.py::TestBoundsPropagation::test_division_positive test_bounds.py::TestBoundsPropagation::test_indexing_2d test_bounds.py::TestBoundsPropagation::test_promote test_bounds.py::TestBoundsPropagation::test_maximum test_bounds.py::TestBoundsPropagation::test_minimum test_bounds.py::TestBoundsPropagation::test_matmul_variable_intervals test_bounds.py::TestBoundsPropagation::test_composed_expression test_bounds.py::TestBoundsPropagation::test_complex_expression
test_that("get_bounds covers additional CVXPY propagation cases", {
  gb <- function(e) lapply(get_bounds(e), as.numeric)

  cst <- Constant(c(1, 2, 3))
  expect_equal(gb(cst), list(c(1, 2, 3), c(1, 2, 3)))

  p <- Parameter(3, bounds = list(-1, 1))
  expect_equal(gb(p), list(rep(-1, 3), rep(1, 3)))

  xp <- Variable(3, nonneg = TRUE)
  expect_equal(gb(abs(xp)), list(rep(0, 3), rep(Inf, 3)))

  x <- Variable(3, bounds = list(1, 2))
  y <- Variable(3, bounds = list(3, 4))
  expect_equal(gb(x * y), list(rep(3, 3), rep(8, 3)))

  xs <- Variable(3, bounds = list(-1, 2))
  ys <- Variable(3, bounds = list(-3, 4))
  expect_equal(gb(xs * ys), list(rep(-6, 3), rep(8, 3)))

  xd <- Variable(3, bounds = list(2, 4))
  yd <- Variable(3, bounds = list(1, 2))
  expect_equal(gb(xd / yd), list(rep(1, 3), rep(4, 3)))

  X <- Variable(c(3, 4), bounds = list(1, 5))
  expect_equal(gb(X[1:2, 2:3]), list(rep(1, 4), rep(5, 4)))

  s <- Variable(bounds = list(2, 5))
  z <- Variable(3, bounds = list(0, 0))
  expect_equal(gb(s + z), list(rep(2, 3), rep(5, 3)))

  xm <- Variable(3, bounds = list(1, 2))
  ym <- Variable(3, bounds = list(0, 3))
  expect_equal(gb(Maximum(xm, ym)), list(rep(1, 3), rep(3, 3)))
  expect_equal(gb(Minimum(xm, ym)), list(rep(0, 3), rep(2, 3)))

  A <- Variable(c(2, 3), bounds = list(0, 1))
  B <- Variable(c(3, 4), bounds = list(0, 1))
  mm <- get_bounds(A %*% B)
  expect_equal(dim(mm[[1]]), c(2L, 4L))
  expect_true(all(is.infinite(mm[[1]]) & mm[[1]] < 0))
  expect_true(all(is.infinite(mm[[2]]) & mm[[2]] > 0))

  xc <- Variable(bounds = list(0, 1))
  expect_equal(gb(2 * xc + 1), list(1, 3))

  xcx <- Variable(3, bounds = list(0, 1))
  ycx <- Variable(3, bounds = list(1, 2))
  expect_equal(gb(abs(xcx - ycx)), list(rep(0, 3), rep(2, 3)))
})

## @cvxpy test_bounds.py::TestBoundsPropagation::test_scalar_bounds_broadcast
test_that("scalar bounds broadcast through bounded atoms", {
  x <- Variable(c(10, 10), bounds = list(-1, 1))
  for (expr in list(exp(x), abs(x), -x, t(x))) {
    b <- get_bounds(expr)
    expect_equal(dim(b[[1L]]), expr@shape)
    expect_equal(dim(b[[2L]]), expr@shape)
  }

  xp <- Variable(c(10, 10), bounds = list(1, 4))
  for (expr in list(log(xp), sqrt(xp))) {
    b <- get_bounds(expr)
    expect_equal(dim(b[[1L]]), expr@shape)
    expect_equal(dim(b[[2L]]), expr@shape)
  }
})

## @cvxpy test_bounds.py::TestBoundsUtilityFunctions::test_unbounded test_bounds.py::TestBoundsUtilityFunctions::test_scalar_bounds test_bounds.py::TestBoundsUtilityFunctions::test_arithmetic_bounds test_bounds.py::TestBoundsUtilityFunctions::test_mul_bounds test_bounds.py::TestBoundsUtilityFunctions::test_mul_bounds_zero_times_inf test_bounds.py::TestBoundsUtilityFunctions::test_div_bounds test_bounds.py::TestBoundsUtilityFunctions::test_abs_bounds test_bounds.py::TestBoundsUtilityFunctions::test_monotonic_bounds test_bounds.py::TestBoundsUtilityFunctions::test_power_bounds test_bounds.py::TestBoundsUtilityFunctions::test_maximum_bounds test_bounds.py::TestBoundsUtilityFunctions::test_minimum_bounds test_bounds.py::TestBoundsUtilityFunctions::test_reduction_bounds test_bounds.py::TestBoundsUtilityFunctions::test_norm1_bounds test_bounds.py::TestBoundsUtilityFunctions::test_norm_inf_bounds test_bounds.py::TestBoundsUtilityFunctions::test_broadcast_bounds test_bounds.py::TestBoundsUtilityFunctions::test_reshape_bounds test_bounds.py::TestBoundsUtilityFunctions::test_transpose_bounds test_bounds.py::TestBoundsUtilityFunctions::test_index_bounds test_bounds.py::TestBoundsUtilityFunctions::test_matmul_bounds_constant_left test_bounds.py::TestBoundsUtilityFunctions::test_matmul_bounds_both_intervals test_bounds.py::TestBoundsUtilityFunctions::test_refine_bounds_from_sign
test_that("bounds utility functions mirror CVXPY 1.9", {
  bp <- function(b) lapply(b, as.numeric)

  expect_equal(bp(unbounded(c(2L, 3L))), list(rep(-Inf, 6), rep(Inf, 6)))
  expect_equal(bp(scalar_bounds(-1, 2)), list(-1, 2))
  expect_equal(bp(add_bounds(c(1, 2), c(3, 4), c(5, 6), c(7, 8))),
               list(c(6, 8), c(10, 12)))
  expect_equal(bp(neg_bounds(c(1, 2), c(3, 4))), list(c(-3, -4), c(-1, -2)))
  expect_equal(bp(scale_bounds(c(1, 2), c(3, 4), -2)), list(c(-6, -8), c(-2, -4)))

  expect_equal(bp(mul_bounds(c(1, 2), c(3, 4), c(5, 6), c(7, 8))),
               list(c(5, 12), c(21, 32)))
  expect_equal(bp(mul_bounds(0, Inf, 0, Inf)), list(0, Inf))
  expect_equal(bp(mul_bounds(0, 1, -Inf, Inf)), list(-Inf, Inf))
  expect_equal(bp(mul_bounds(-Inf, 0, 0, Inf)), list(-Inf, 0))

  expect_equal(bp(div_bounds(c(2, 4), c(6, 8), c(1, 2), c(2, 4))),
               list(c(1, 1), c(6, 4)))
  expect_equal(bp(div_bounds(1, 2, -1, 1)), list(-Inf, Inf))

  expect_equal(bp(abs_bounds(c(-2, 1), c(3, 2))), list(c(0, 1), c(3, 2)))
  expect_equal(bp(exp_bounds(c(0, 1), c(1, 2))), list(c(1, exp(1)), c(exp(1), exp(2))))
  expect_equal(bp(log_bounds(c(1, exp(1)), c(exp(1), exp(2)))), list(c(0, 1), c(1, 2)))
  expect_equal(bp(sqrt_bounds(c(1, 4), c(4, 9))), list(c(1, 2), c(2, 3)))
  expect_equal(bp(power_bounds(c(2, -3), c(4, -1), 2)), list(c(4, 1), c(16, 9)))

  bl <- list(list(c(1, 2), c(3, 4)), list(c(0, 5), c(2, 6)))
  expect_equal(bp(maximum_bounds(bl)), list(c(1, 5), c(3, 6)))
  expect_equal(bp(minimum_bounds(bl)), list(c(0, 2), c(2, 4)))

  lb <- matrix(c(1, 3, 2, 4), 2, 2)
  ub <- matrix(c(5, 7, 6, 8), 2, 2)
  expect_equal(bp(sum_bounds(lb, ub)), list(10, 26))
  expect_equal(bp(max_reduction_bounds(lb, ub)), list(4, 8))
  expect_equal(bp(min_reduction_bounds(lb, ub)), list(1, 5))
  expect_equal(bp(norm1_bounds(c(1, 2, 3), c(2, 3, 4))), list(6, 9))
  expect_equal(bp(norm_inf_bounds(c(1, 2, 3), c(2, 5, 4))), list(3, 5))

  bb <- broadcast_bounds(1, 2, c(3L, 4L))
  expect_equal(dim(bb[[1]]), c(3L, 4L))
  expect_equal(bp(bb), list(rep(1, 12), rep(2, 12)))
  expect_equal(bp(reshape_bounds(matrix(c(1, 3, 2, 4), 2, 2),
                                 matrix(c(5, 7, 6, 8), 2, 2), c(4L, 1L))),
               list(c(1, 3, 2, 4), c(5, 7, 6, 8)))
  tb <- transpose_bounds(matrix(1:6, 2, 3), matrix(7:12, 2, 3))
  expect_equal(dim(tb[[1]]), c(3L, 2L))
  expect_equal(bp(index_bounds(1:5, 6:10, function(a) a[2:4])), list(2:4, 7:9))

  A <- matrix(c(1, 3, 2, 4), 2, 2)
  bmat <- matrix(1, 2, 2)
  expect_equal(bp(matmul_bounds(A, A, bmat, 2 * bmat)), list(as.numeric(A %*% bmat),
                                                             as.numeric(A %*% (2 * bmat))))
  mb <- matmul_bounds(matrix(c(1, 3, 2, 4), 2, 2), matrix(c(2, 4, 3, 5), 2, 2),
                      bmat, 2 * bmat)
  expect_true(all(is.infinite(mb[[1]]) & mb[[1]] < 0))
  expect_true(all(is.infinite(mb[[2]]) & mb[[2]] > 0))

  expect_equal(bp(refine_bounds_from_sign(c(-1, -2), c(3, 4), TRUE, FALSE)),
               list(c(0, 0), c(3, 4)))
  expect_equal(bp(refine_bounds_from_sign(c(-1, -2), c(3, 4), FALSE, TRUE)),
               list(c(-1, -2), c(0, 0)))
})

## @cvxpy test_bounds.py::TestHIGHSBoundsPropagation::test_get_problem_data_scalar_bounds test_bounds.py::TestHIGHSBoundsPropagation::test_get_problem_data_nonneg_bounds test_bounds.py::TestHIGHSBoundsPropagation::test_get_problem_data_nonpos_bounds test_bounds.py::TestHIGHSBoundsPropagation::test_get_problem_data_combined_nonneg_and_bounds test_bounds.py::TestHIGHSBoundsPropagation::test_get_problem_data_combined_nonpos_and_bounds test_bounds.py::TestHIGHSBoundsPropagation::test_get_problem_data_matrix_bounds test_bounds.py::TestHIGHSBoundsPropagation::test_get_problem_data_multiple_variables
test_that("HiGHS problem_data receives native dense variable bounds", {
  skip_if_not_installed("highs")

  x <- Variable(3, bounds = list(-2, 5))
  data <- problem_data(Problem(Minimize(sum_entries(x)), list()), HIGHS_SOLVER)$data
  expect_equal(data[[LOWER_BOUNDS]], rep(-2, 3))
  expect_equal(data[[UPPER_BOUNDS]], rep(5, 3))

  xn <- Variable(3, nonneg = TRUE)
  data <- problem_data(Problem(Minimize(sum_entries(xn)), list()), HIGHS_SOLVER)$data
  expect_equal(data[[LOWER_BOUNDS]], rep(0, 3))
  expect_null(data[[UPPER_BOUNDS]])

  xp <- Variable(3, nonpos = TRUE)
  data <- problem_data(Problem(Maximize(sum_entries(xp)), list()), HIGHS_SOLVER)$data
  expect_null(data[[LOWER_BOUNDS]])
  expect_equal(data[[UPPER_BOUNDS]], rep(0, 3))

  xt <- Variable(3, nonneg = TRUE, bounds = list(2, 10))
  data <- problem_data(Problem(Minimize(sum_entries(xt)), list()), HIGHS_SOLVER)$data
  expect_equal(data[[LOWER_BOUNDS]], rep(2, 3))
  expect_equal(data[[UPPER_BOUNDS]], rep(10, 3))

  xu <- Variable(3, nonpos = TRUE, bounds = list(-10, -2))
  data <- problem_data(Problem(Maximize(sum_entries(xu)), list()), HIGHS_SOLVER)$data
  expect_equal(data[[LOWER_BOUNDS]], rep(-10, 3))
  expect_equal(data[[UPPER_BOUNDS]], rep(-2, 3))

  X <- Variable(c(2, 3), bounds = list(-1, 1))
  data <- problem_data(Problem(Minimize(sum_entries(X)), list()), HIGHS_SOLVER)$data
  expect_equal(length(data[[LOWER_BOUNDS]]), 6L)
  expect_equal(data[[LOWER_BOUNDS]], rep(-1, 6))
  expect_equal(data[[UPPER_BOUNDS]], rep(1, 6))

  y <- Variable(2, bounds = list(c(1, 2), c(10, 20)))
  data <- problem_data(Problem(Minimize(sum_entries(x[1:2] + y)), list()), HIGHS_SOLVER)$data
  expect_true(all(c(1, 2) %in% data[[LOWER_BOUNDS]]))
  expect_true(all(c(10, 20) %in% data[[UPPER_BOUNDS]]))
})

## @cvxpy test_parametric_bounds.py::TestParametricBoundsCreation::test_expression_bounds test_parametric_bounds.py::TestParametricBoundsCreation::test_mixed_none_and_param test_parametric_bounds.py::TestParametricBoundsCreation::test_scalar_param_broadcast test_parametric_bounds.py::TestParametricBoundsCreation::test_has_lower_upper_bounds
test_that("parametric bounds creation covers expressions and one-sided bounds", {
  p <- Parameter(value = 2)
  x <- Variable(bounds = list(p + 1, 2 * p))
  expect_false(is.null(x@attributes$bounds))

  ub <- Parameter(value = 5)
  one_sided <- Variable(bounds = list(NULL, ub))
  b <- get_bounds(one_sided)
  expect_equal(as.numeric(b[[1]]), -Inf)
  expect_equal(as.numeric(b[[2]]), Inf)
  expect_true(S7::S7_inherits(one_sided@attributes$bounds[[2L]], Expression))

  v <- Variable(3, bounds = list(Parameter(), Parameter()))
  expect_false(is.null(v@attributes$bounds))

  lower <- Variable(bounds = list(Parameter(), 10))
  upper <- Variable(bounds = list(0, Parameter()))
  expect_true(S7::S7_inherits(lower@attributes$bounds[[1L]], Expression))
  expect_true(S7::S7_inherits(upper@attributes$bounds[[2L]], Expression))
})

## @cvxpy test_parametric_bounds.py::TestParametricBoundsSolving::test_minimize_maximize test_parametric_bounds.py::TestParametricBoundsSolving::test_expression_bounds test_parametric_bounds.py::TestParametricBoundsSolving::test_mixed_numeric_and_param test_parametric_bounds.py::TestParametricBoundsSolving::test_one_sided_param_bound test_parametric_bounds.py::TestParametricBoundsSolving::test_vector_variable_scalar_param test_parametric_bounds.py::TestParametricBoundsSolving::test_multiple_variables test_parametric_bounds.py::TestParametricBoundsEdgeCases::test_equal_bounds_pins_value test_parametric_bounds.py::TestParametricBoundsEdgeCases::test_infeasible_bounds test_parametric_bounds.py::TestParametricBoundsEdgeCases::test_unset_parameter_raises test_parametric_bounds.py::TestParametricBoundsEdgeCases::test_infinite_bounds_unbounded
test_that("parametric bounds solve edge cases through Clarabel", {
  lb <- Parameter(value = 2)
  ub <- Parameter(value = 5)
  x <- Variable(bounds = list(lb, ub))
  prob_min <- Problem(Minimize(x), list())
  expect_equal(psolve(prob_min, solver = "CLARABEL"), 2, tolerance = 1e-5)
  prob_max <- Problem(Maximize(x), list())
  expect_equal(psolve(prob_max, solver = "CLARABEL"), 5, tolerance = 1e-5)

  scale <- Parameter(nonneg = TRUE, value = 10)
  xe <- Variable(bounds = list(-scale, scale))
  expect_equal(psolve(Problem(Minimize(xe), list()), solver = "CLARABEL"), -10,
               tolerance = 1e-5)

  up <- Parameter(value = 5)
  xu <- Variable(bounds = list(0, up))
  expect_equal(psolve(Problem(Maximize(xu), list()), solver = "CLARABEL"), 5,
               tolerance = 1e-5)

  xo <- Variable(bounds = list(NULL, up))
  expect_equal(psolve(Problem(Maximize(xo), list()), solver = "CLARABEL"), 5,
               tolerance = 1e-5)

  lv <- Parameter(value = -2)
  xv <- Variable(3, bounds = list(lv, NULL))
  prob_vec <- Problem(Minimize(sum_entries(xv)), list(xv <= 10))
  expect_equal(psolve(prob_vec, solver = "CLARABEL"), -6, tolerance = 1e-5)
  expect_equal(as.numeric(value(xv)), rep(-2, 3), tolerance = 1e-5)

  l1 <- Parameter(value = 1)
  l2 <- Parameter(value = 2)
  x1 <- Variable(bounds = list(l1, NULL))
  x2 <- Variable(bounds = list(l2, NULL))
  expect_equal(psolve(Problem(Minimize(x1 + x2), list()), solver = "CLARABEL"), 3,
               tolerance = 1e-5)

  b <- Parameter(value = 7)
  pinned <- Variable(bounds = list(b, b))
  expect_equal(psolve(Problem(Minimize(pinned), list()), solver = "CLARABEL"), 7,
               tolerance = 1e-5)

  bad_lb <- Parameter(value = 10)
  bad_ub <- Parameter(value = 5)
  infeas <- Variable(bounds = list(bad_lb, bad_ub))
  infeas_prob <- Problem(Minimize(infeas), list())
  psolve(infeas_prob, solver = "CLARABEL")
  expect_true(status(infeas_prob) %in% c("infeasible", "infeasible_inaccurate"))

  unset <- Variable(bounds = list(Parameter(), 10))
  expect_error(psolve(Problem(Minimize(unset), list()), solver = "CLARABEL"))

  ub_inf <- Parameter(value = 5)
  unbounded <- Variable(bounds = list(-Inf, ub_inf))
  unbounded_prob <- Problem(Minimize(unbounded), list())
  psolve(unbounded_prob, solver = "CLARABEL")
  expect_true(status(unbounded_prob) %in% c("unbounded", "unbounded_inaccurate"))
})

## @cvxpy test_parametric_bounds.py::TestParametricBoundsDPP::test_is_dpp_and_chain
test_that("parametric bounds participate in DPP", {
  lb <- Parameter()
  ub <- Parameter()
  x <- Variable(bounds = list(lb, ub))
  prob <- Problem(Minimize(x), list())
  expect_true(is_dpp(prob))

  value(lb) <- 1
  value(ub) <- 10
  expect_equal(psolve(prob, solver = "CLARABEL"), 1, tolerance = 1e-5)
  value(lb) <- 5
  expect_equal(psolve(prob, solver = "CLARABEL"), 5, tolerance = 1e-5)

  scale <- Parameter(nonneg = TRUE, value = 10)
  xe <- Variable(bounds = list(-scale, scale))
  pe <- Problem(Minimize(xe), list())
  expect_true(is_dpp(pe))
  expect_equal(psolve(pe, solver = "CLARABEL"), -10, tolerance = 1e-5)

  p1 <- Parameter(value = 1)
  p2 <- Parameter(value = 1)
  bad <- Variable(bounds = list(p1 * p2, NULL))
  bad_prob <- Problem(Minimize(bad), list(bad <= 10))
  expect_false(is_dpp(bad_prob))
})

## @cvxpy test_bounds.py::TestDiagBounds::test_diag_variable_with_scalar_bounds test_bounds.py::TestDiagBounds::test_diag_variable_scalar_bounds_must_contain_zero test_bounds.py::TestDiagBounds::test_diag_variable_rejects_dense_array_bounds test_bounds.py::TestDiagBounds::test_diag_variable_rejects_expression_bounds test_bounds.py::TestDiagBounds::test_diag_variable_with_scalar_bounds_solves
test_that("diagonal variables validate and solve scalar bounds", {
  x <- Variable(c(3, 3), diag = TRUE, bounds = list(-1, 1))
  expect_false(is.null(x@attributes$bounds))

  expect_error(Variable(c(3, 3), diag = TRUE, bounds = list(1, 2)),
               "Scalar lower bound")
  expect_error(Variable(c(3, 3), diag = TRUE, bounds = list(-2, -1)),
               "Scalar upper bound")
  expect_false(is.null(Variable(c(3, 3), diag = TRUE, bounds = list(0, 1))@attributes$bounds))
  expect_false(is.null(Variable(c(3, 3), diag = TRUE, bounds = list(-1, 0))@attributes$bounds))

  expect_error(Variable(c(3, 3), diag = TRUE,
                        bounds = list(-matrix(1, 3, 3), matrix(1, 3, 3))),
               "Dense array bounds")

  p <- Parameter(value = 1)
  expect_error(Variable(c(3, 3), diag = TRUE, bounds = list(p, 1)),
               "Expression bounds")

  xd <- Variable(c(3, 3), diag = TRUE, bounds = list(-1, 2))
  prob <- Problem(Minimize(matrix_trace(xd)), list())
  expect_equal(psolve(prob, solver = "CLARABEL"), -3, tolerance = 1e-5)
})
