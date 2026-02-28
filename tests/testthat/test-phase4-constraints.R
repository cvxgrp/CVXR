## Tests for Phase 4: Cone Constraints
## SOC, PSD, ExpCone, PowCone3D, PowConeND

library(testthat)

# ══════════════════════════════════════════════════════════════════
# Utility: .reshape_c_order
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that(".reshape_c_order produces row-major fill", {
  ## C-order: 1,2,3 fill across columns first (row-major)
  m <- CVXR:::.reshape_c_order(1:6, 2L, 3L)
  expect_equal(m, matrix(1:6, 2, 3, byrow = TRUE))
})

## @cvxpy NONE
test_that(".reshape_c_order handles single element", {
  m <- CVXR:::.reshape_c_order(42, 1L, 1L)
  expect_equal(m, matrix(42, 1, 1))
})

## @cvxpy NONE
test_that(".reshape_c_order handles rectangular matrix", {
  m <- CVXR:::.reshape_c_order(1:12, 3L, 4L)
  expect_equal(m, matrix(1:12, 3, 4, byrow = TRUE))
})

# ══════════════════════════════════════════════════════════════════
# Settings constants
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Cone dimension keys are defined", {
  expect_equal(EQ_DIM, "f")
  expect_equal(LEQ_DIM, "l")
  expect_equal(SOC_DIM, "q")
  expect_equal(PSD_DIM, "s")
  expect_equal(EXP_DIM, "ep")
})

## @cvxpy NONE
test_that("Constraint type IDs are defined", {
  expect_equal(CVXR:::CONSTR_ID_EQ, 0L)
  expect_equal(CVXR:::CONSTR_ID_SOC, 2L)
  expect_equal(CVXR:::CONSTR_ID_PSD, 4L)
  expect_equal(CVXR:::CONSTR_ID_EXP, 5L)
})

# ══════════════════════════════════════════════════════════════════
# Cone base class
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Cone inherits from Constraint", {
  x <- Variable(c(2L, 1L))
  constr <- SOC(Variable(1L), x)
  expect_true(S7::S7_inherits(constr, Cone))
  expect_true(S7::S7_inherits(constr, Constraint))
})

# ══════════════════════════════════════════════════════════════════
# SOC
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("SOC construction with scalar t", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  constr <- SOC(t_var, x_var)
  expect_true(S7::S7_inherits(constr, SOC))
  expect_true(S7::S7_inherits(constr, Cone))
  expect_equal(constr@axis, 2L)
  expect_equal(length(constr@args), 2L)
})

## @cvxpy NONE
test_that("SOC construction with vector t", {
  t_var <- Variable(c(2L, 1L))
  X_var <- Variable(c(3L, 2L))
  constr <- SOC(t_var, X_var, axis = 2L)
  expect_true(S7::S7_inherits(constr, SOC))
  expect_equal(num_cones(constr), 2L)
})

## @cvxpy NONE
test_that("SOC with axis=1", {
  t_var <- Variable(c(3L, 1L))
  X_var <- Variable(c(3L, 2L))
  constr <- SOC(t_var, X_var, axis = 1L)
  expect_equal(constr@axis, 1L)
  expect_equal(num_cones(constr), 3L)
})

## @cvxpy NONE
test_that("SOC rejects matrix t", {
  t_var <- Variable(c(2L, 2L))
  X_var <- Variable(c(3L, 2L))
  expect_error(SOC(t_var, X_var), "must be a vector")
})

## @cvxpy NONE
test_that("SOC rejects dimension mismatch", {
  t_var <- Variable(c(3L, 1L))
  X_var <- Variable(c(3L, 2L))
  ## t has 3 entries but X has 2 columns (axis=0)
  expect_error(SOC(t_var, X_var, axis = 2L), "incompatible")
})

## @cvxpy NONE
test_that("SOC is_dcp with affine args", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  constr <- SOC(t_var, x_var)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("SOC is_dgp always FALSE", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  constr <- SOC(t_var, x_var)
  expect_false(is_dgp(constr))
})

## @cvxpy NONE
test_that("SOC expr_name", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  constr <- SOC(t_var, x_var)
  nm <- expr_name(constr)
  expect_true(grepl("SOC", nm))
})

## @cvxpy NONE
test_that("SOC get_data returns axis and id", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  constr <- SOC(t_var, x_var, axis = 2L)
  d <- get_data(constr)
  expect_equal(d[[1]], 2L)
  expect_equal(d[[2]], constr@id)
})

## @cvxpy NONE
test_that("SOC num_cones for scalar t", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  constr <- SOC(t_var, x_var)
  expect_equal(num_cones(constr), 1L)
})

## @cvxpy NONE
test_that("SOC cone_sizes", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  constr <- SOC(t_var, x_var)
  ## cone_size = 1 + X_dim = 1 + 3 = 4
  expect_equal(cone_sizes(constr), 4L)
})

## @cvxpy NONE
test_that("SOC constr_size", {
  t_var <- Variable(c(2L, 1L))
  X_var <- Variable(c(3L, 2L))
  constr <- SOC(t_var, X_var, axis = 2L)
  ## cone_size = 1 + 3 = 4; num_cones = 2; size = 8
  expect_equal(constr_size(constr), 8L)
})

## @cvxpy NONE
test_that("SOC residual when feasible", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  value(t_var) <- 5
  value(x_var) <- matrix(c(1, 2, 2), 3, 1)  ## norm = 3 < 5
  constr <- SOC(t_var, x_var)
  r <- residual(constr)
  expect_true(!is.null(r))
  expect_equal(r, 0, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("SOC residual when infeasible", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  value(t_var) <- 1
  value(x_var) <- matrix(c(3, 4, 0), 3, 1)  ## norm = 5 > 1
  constr <- SOC(t_var, x_var)
  r <- residual(constr)
  expect_true(r > 0)
})

## @cvxpy NONE
test_that("SOC residual is NULL when no value", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  constr <- SOC(t_var, x_var)
  expect_null(residual(constr))
})

## @cvxpy NONE
test_that("SOC dual_cone is self-dual", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  constr <- SOC(t_var, x_var)
  dc <- dual_cone(constr)
  expect_true(S7::S7_inherits(dc, SOC))
  expect_equal(dc@axis, constr@axis)
})

## @cvxpy NONE
test_that("SOC save_dual_value dispatch", {
  t_var <- Variable(1L)
  x_var <- Variable(c(2L, 1L))
  constr <- SOC(t_var, x_var)
  ## cone_size = 1 + 2 = 3; num_cones = 1; total = 3 entries
  save_dual_value(constr, c(10, 20, 30))
  tv <- value(constr@dual_variables[[1L]])
  xv <- value(constr@dual_variables[[2L]])
  expect_equal(as.numeric(tv), 10)
})

# ══════════════════════════════════════════════════════════════════
# PSD
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PSD construction with square matrix", {
  X <- Variable(c(3L, 3L))
  constr <- PSD(X)
  expect_true(S7::S7_inherits(constr, PSD))
  expect_true(S7::S7_inherits(constr, Cone))
  expect_true(S7::S7_inherits(constr, Constraint))
  expect_equal(constr@shape, c(3L, 3L))
})

## @cvxpy NONE
test_that("PSD rejects non-square matrix", {
  X <- Variable(c(2L, 3L))
  expect_error(PSD(X), "Non-square")
})

## @cvxpy NONE
test_that("PSD is_dcp with affine arg", {
  X <- Variable(c(3L, 3L))
  constr <- PSD(X)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("PSD is_dgp always FALSE", {
  X <- Variable(c(3L, 3L))
  constr <- PSD(X)
  expect_false(is_dgp(constr))
})

## @cvxpy NONE
test_that("PSD expr_name contains >>", {
  X <- Variable(c(3L, 3L))
  constr <- PSD(X)
  nm <- expr_name(constr)
  expect_true(grepl(">>", nm))
})

## @cvxpy NONE
test_that("PSD num_cones is 1", {
  X <- Variable(c(3L, 3L))
  constr <- PSD(X)
  expect_equal(num_cones(constr), 1L)
})

## @cvxpy NONE
test_that("PSD cone_sizes returns n", {
  X <- Variable(c(3L, 3L))
  constr <- PSD(X)
  expect_equal(cone_sizes(constr), 3L)
})

## @cvxpy NONE
test_that("PSD constr_size returns n*n", {
  X <- Variable(c(3L, 3L))
  constr <- PSD(X)
  expect_equal(constr_size(constr), 9L)
})

## @cvxpy NONE
test_that("PSD residual for PSD matrix is 0", {
  X <- Variable(c(2L, 2L))
  ## diag(2,2) is PSD
  value(X) <- diag(2)
  constr <- PSD(X)
  r <- residual(constr)
  expect_equal(r, 0, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("PSD residual for non-PSD matrix is > 0", {
  X <- Variable(c(2L, 2L))
  ## [[1,0],[0,-1]] has eigenvalue -1
  value(X) <- matrix(c(1, 0, 0, -1), 2, 2)
  constr <- PSD(X)
  r <- residual(constr)
  expect_true(r > 0)
  expect_equal(r, 1, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("PSD residual is NULL when no value", {
  X <- Variable(c(2L, 2L))
  constr <- PSD(X)
  expect_null(residual(constr))
})

## @cvxpy NONE
test_that("PSD value satisfied for PSD matrix", {
  X <- Variable(c(2L, 2L))
  value(X) <- diag(2)
  constr <- PSD(X)
  expect_true(value(constr))
})

## @cvxpy NONE
test_that("PSD value violated for non-PSD matrix", {
  X <- Variable(c(2L, 2L))
  value(X) <- matrix(c(1, 0, 0, -1), 2, 2)
  constr <- PSD(X)
  expect_false(value(constr))
})

## @cvxpy NONE
test_that("PSD dual_cone is self-dual", {
  X <- Variable(c(3L, 3L))
  constr <- PSD(X)
  dc <- dual_cone(constr)
  expect_true(S7::S7_inherits(dc, PSD))
})

# ══════════════════════════════════════════════════════════════════
# ExpCone
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ExpCone construction", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  z <- Variable(c(2L, 1L))
  constr <- ExpCone(x, y, z)
  expect_true(S7::S7_inherits(constr, ExpCone))
  expect_true(S7::S7_inherits(constr, Cone))
  expect_equal(length(constr@args), 3L)
  expect_equal(length(constr@dual_variables), 3L)
})

## @cvxpy NONE
test_that("ExpCone shape override", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  z <- Variable(c(2L, 1L))
  constr <- ExpCone(x, y, z)
  ## shape = c(3, prod(x_shape)) = c(3, 2)
  expect_equal(constr@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("ExpCone shape override scalar", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- ExpCone(x, y, z)
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("ExpCone rejects shape mismatch", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(3L, 1L))
  z <- Variable(c(2L, 1L))
  expect_error(ExpCone(x, y, z), "same shapes")
})

## @cvxpy NONE
test_that("ExpCone rejects non-affine", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  z <- square(Variable(c(2L, 1L)))  ## non-affine
  expect_error(ExpCone(x, y, z), "affine and real")
})

## @cvxpy NONE
test_that("ExpCone is_dcp with affine args", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- ExpCone(x, y, z)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("ExpCone is_dgp always FALSE", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- ExpCone(x, y, z)
  expect_false(is_dgp(constr))
})

## @cvxpy NONE
test_that("ExpCone num_cones", {
  x <- Variable(c(3L, 1L))
  y <- Variable(c(3L, 1L))
  z <- Variable(c(3L, 1L))
  constr <- ExpCone(x, y, z)
  expect_equal(num_cones(constr), 3L)
})

## @cvxpy NONE
test_that("ExpCone cone_sizes", {
  x <- Variable(c(3L, 1L))
  y <- Variable(c(3L, 1L))
  z <- Variable(c(3L, 1L))
  constr <- ExpCone(x, y, z)
  expect_equal(cone_sizes(constr), rep(3L, 3L))
})

## @cvxpy NONE
test_that("ExpCone constr_size", {
  x <- Variable(c(3L, 1L))
  y <- Variable(c(3L, 1L))
  z <- Variable(c(3L, 1L))
  constr <- ExpCone(x, y, z)
  expect_equal(constr_size(constr), 9L)
})

## @cvxpy NONE
test_that("ExpCone residual feasible", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  value(x) <- 0
  value(y) <- 1
  value(z) <- 2  ## y*exp(x/y) = 1*exp(0) = 1 <= 2
  constr <- ExpCone(x, y, z)
  r <- residual(constr)
  expect_equal(r, 0, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("ExpCone residual infeasible", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  value(x) <- 2
  value(y) <- 1
  value(z) <- 1  ## y*exp(x/y) = exp(2) ~ 7.39 > 1
  constr <- ExpCone(x, y, z)
  r <- residual(constr)
  expect_true(r > 0)
})

## @cvxpy NONE
test_that("ExpCone residual NULL when no value", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- ExpCone(x, y, z)
  expect_null(residual(constr))
})

## @cvxpy NONE
test_that("ExpCone expr_name", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- ExpCone(x, y, z)
  nm <- expr_name(constr)
  expect_true(grepl("ExpCone", nm))
})

## @cvxpy NONE
test_that("ExpCone dual_cone mapping", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- ExpCone(x, y, z)
  dc <- dual_cone(constr)
  expect_true(S7::S7_inherits(dc, ExpCone))
})

## @cvxpy NONE
test_that("ExpCone save_dual_value dispatch", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- ExpCone(x, y, z)
  ## 3 values for 1 cone
  save_dual_value(constr, c(10, 20, 30))
  dv <- value(constr@dual_variables[[1L]])
  expect_equal(as.numeric(dv), 10)
  dv2 <- value(constr@dual_variables[[2L]])
  expect_equal(as.numeric(dv2), 20)
  dv3 <- value(constr@dual_variables[[3L]])
  expect_equal(as.numeric(dv3), 30)
})

# ══════════════════════════════════════════════════════════════════
# PowCone3D
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PowCone3D construction", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  z <- Variable(c(2L, 1L))
  constr <- PowCone3D(x, y, z, 0.5)
  expect_true(S7::S7_inherits(constr, PowCone3D))
  expect_true(S7::S7_inherits(constr, Cone))
  expect_equal(length(constr@args), 3L)
})

## @cvxpy NONE
test_that("PowCone3D scalar alpha promotion", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  z <- Variable(c(2L, 1L))
  constr <- PowCone3D(x, y, z, 0.3)
  ## alpha should be promoted to match x shape
  alpha_val <- value(constr@alpha)
  expect_equal(dim(alpha_val), c(2L, 1L))
  expect_equal(as.numeric(alpha_val), c(0.3, 0.3))
})

## @cvxpy NONE
test_that("PowCone3D alpha out of (0,1) errors", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  expect_error(PowCone3D(x, y, z, 0), "open interval")
  expect_error(PowCone3D(x, y, z, 1), "open interval")
  expect_error(PowCone3D(x, y, z, -0.5), "open interval")
  expect_error(PowCone3D(x, y, z, 1.5), "open interval")
})

## @cvxpy NONE
test_that("PowCone3D rejects shape mismatch", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(3L, 1L))
  z <- Variable(c(2L, 1L))
  expect_error(PowCone3D(x, y, z, 0.5), "same shapes")
})

## @cvxpy NONE
test_that("PowCone3D is_dcp", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- PowCone3D(x, y, z, 0.5)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("PowCone3D is_dgp FALSE", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- PowCone3D(x, y, z, 0.5)
  expect_false(is_dgp(constr))
})

## @cvxpy NONE
test_that("PowCone3D get_data", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- PowCone3D(x, y, z, 0.5)
  d <- get_data(constr)
  expect_true(S7::S7_inherits(d[[1]], CVXR:::Expression))  ## alpha
  expect_equal(d[[2]], constr@id)
})

## @cvxpy NONE
test_that("PowCone3D shape override", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  z <- Variable(c(2L, 1L))
  constr <- PowCone3D(x, y, z, 0.5)
  expect_equal(constr@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("PowCone3D num_cones", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  z <- Variable(c(2L, 1L))
  constr <- PowCone3D(x, y, z, 0.5)
  expect_equal(num_cones(constr), 2L)
})

## @cvxpy NONE
test_that("PowCone3D cone_sizes", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  z <- Variable(c(2L, 1L))
  constr <- PowCone3D(x, y, z, 0.5)
  expect_equal(cone_sizes(constr), c(3L, 3L))
})

## @cvxpy NONE
test_that("PowCone3D constr_size", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  z <- Variable(c(2L, 1L))
  constr <- PowCone3D(x, y, z, 0.5)
  expect_equal(constr_size(constr), 6L)
})

## @cvxpy NONE
test_that("PowCone3D residual feasible", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  value(x) <- 4
  value(y) <- 4
  value(z) <- 2  ## 4^0.5 * 4^0.5 = 4 >= |2|
  constr <- PowCone3D(x, y, z, 0.5)
  r <- residual(constr)
  expect_equal(r, 0, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("PowCone3D residual infeasible", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  value(x) <- 1
  value(y) <- 1
  value(z) <- 5  ## 1^0.5 * 1^0.5 = 1 < |5|
  constr <- PowCone3D(x, y, z, 0.5)
  r <- residual(constr)
  expect_true(r > 0)
})

## @cvxpy NONE
test_that("PowCone3D dual_cone scaling", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- PowCone3D(x, y, z, 0.5)
  dc <- dual_cone(constr)
  expect_true(S7::S7_inherits(dc, PowCone3D))
})

## @cvxpy NONE
test_that("PowCone3D save_dual_value", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- PowCone3D(x, y, z, 0.5)
  ## (3, num_cones=1) → reshape (3, 1) C-order
  save_dual_value(constr, c(10, 20, 30))
  expect_equal(as.numeric(value(constr@dual_variables[[1L]])), 10)
  expect_equal(as.numeric(value(constr@dual_variables[[2L]])), 20)
  expect_equal(as.numeric(value(constr@dual_variables[[3L]])), 30)
})

## @cvxpy NONE
test_that("PowCone3D expr_name", {
  x <- Variable(1L)
  y <- Variable(1L)
  z <- Variable(1L)
  constr <- PowCone3D(x, y, z, 0.5)
  nm <- expr_name(constr)
  expect_true(grepl("PowCone3D", nm))
})

# ══════════════════════════════════════════════════════════════════
# PowConeND
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PowConeND construction", {
  W <- Variable(c(3L, 2L))
  z <- Variable(c(2L, 1L))
  alpha <- Constant(matrix(c(0.2, 0.3, 0.5, 0.4, 0.4, 0.2), 3, 2))
  constr <- PowConeND(W, z, alpha, axis = 2L)
  expect_true(S7::S7_inherits(constr, PowConeND))
  expect_true(S7::S7_inherits(constr, Cone))
})

## @cvxpy NONE
test_that("PowConeND alpha not summing to 1 errors", {
  W <- Variable(c(3L, 1L))
  z <- Variable(1L)
  alpha <- Constant(matrix(c(0.2, 0.3, 0.4), 3, 1))  ## sum = 0.9 != 1
  expect_error(PowConeND(W, z, alpha, axis = 2L), "sum to 1")
})

## @cvxpy NONE
test_that("PowConeND alpha shape != W shape errors", {
  W <- Variable(c(3L, 2L))
  z <- Variable(c(2L, 1L))
  alpha <- Constant(matrix(rep(1/3, 6), 2, 3))  ## 2x3 != 3x2
  expect_error(PowConeND(W, z, alpha), "not equal")
})

## @cvxpy NONE
test_that("PowConeND invalid axis errors", {
  W <- Variable(c(3L, 1L))
  z <- Variable(1L)
  alpha <- Constant(matrix(c(0.2, 0.3, 0.5), 3, 1))
  expect_error(PowConeND(W, z, alpha, axis = 1L), "incompatible")
})

## @cvxpy NONE
test_that("PowConeND is_dcp always TRUE", {
  W <- Variable(c(3L, 2L))
  z <- Variable(c(2L, 1L))
  alpha <- Constant(matrix(c(0.2, 0.3, 0.5, 0.4, 0.4, 0.2), 3, 2))
  constr <- PowConeND(W, z, alpha, axis = 2L)
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("PowConeND is_dgp FALSE", {
  W <- Variable(c(3L, 2L))
  z <- Variable(c(2L, 1L))
  alpha <- Constant(matrix(c(0.2, 0.3, 0.5, 0.4, 0.4, 0.2), 3, 2))
  constr <- PowConeND(W, z, alpha, axis = 2L)
  expect_false(is_dgp(constr))
})

## @cvxpy NONE
test_that("PowConeND get_data", {
  W <- Variable(c(3L, 2L))
  z <- Variable(c(2L, 1L))
  alpha <- Constant(matrix(c(0.2, 0.3, 0.5, 0.4, 0.4, 0.2), 3, 2))
  constr <- PowConeND(W, z, alpha, axis = 2L)
  d <- get_data(constr)
  expect_true(S7::S7_inherits(d[[1]], CVXR:::Expression))  ## alpha
  expect_equal(d[[2]], 2L)  ## axis
  expect_equal(d[[3]], constr@id)
})

## @cvxpy NONE
test_that("PowConeND num_cones", {
  W <- Variable(c(3L, 2L))
  z <- Variable(c(2L, 1L))
  alpha <- Constant(matrix(c(0.2, 0.3, 0.5, 0.4, 0.4, 0.2), 3, 2))
  constr <- PowConeND(W, z, alpha, axis = 2L)
  expect_equal(num_cones(constr), 2L)
})

## @cvxpy NONE
test_that("PowConeND cone_sizes", {
  W <- Variable(c(3L, 2L))
  z <- Variable(c(2L, 1L))
  alpha <- Constant(matrix(c(0.2, 0.3, 0.5, 0.4, 0.4, 0.2), 3, 2))
  constr <- PowConeND(W, z, alpha, axis = 2L)
  ## cone_size = 1 + W.shape[axis+1] = 1 + 3 = 4; 2 cones
  expect_equal(cone_sizes(constr), c(4L, 4L))
})

## @cvxpy NONE
test_that("PowConeND shape override (m+1, n)", {
  W <- Variable(c(3L, 2L))
  z <- Variable(c(2L, 1L))
  alpha <- Constant(matrix(c(0.2, 0.3, 0.5, 0.4, 0.4, 0.2), 3, 2))
  constr <- PowConeND(W, z, alpha, axis = 2L)
  ## axis=0: m=3, n=2; shape = (4, 2)
  expect_equal(constr@shape, c(4L, 2L))
})

## @cvxpy NONE
test_that("PowConeND residual feasible", {
  W <- Variable(c(2L, 1L))
  z <- Variable(1L)
  alpha <- Constant(matrix(c(0.5, 0.5), 2, 1))
  value(W) <- matrix(c(4, 4), 2, 1)  ## 4^0.5 * 4^0.5 = 4
  value(z) <- 2  ## |2| <= 4
  constr <- PowConeND(W, z, alpha, axis = 2L)
  r <- residual(constr)
  expect_equal(r, 0, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("PowConeND residual infeasible", {
  W <- Variable(c(2L, 1L))
  z <- Variable(1L)
  alpha <- Constant(matrix(c(0.5, 0.5), 2, 1))
  value(W) <- matrix(c(1, 1), 2, 1)  ## 1^0.5 * 1^0.5 = 1
  value(z) <- 5  ## |5| > 1
  constr <- PowConeND(W, z, alpha, axis = 2L)
  r <- residual(constr)
  expect_true(r > 0)
})

## @cvxpy NONE
test_that("PowConeND expr_name", {
  W <- Variable(c(3L, 2L))
  z <- Variable(c(2L, 1L))
  alpha <- Constant(matrix(c(0.2, 0.3, 0.5, 0.4, 0.4, 0.2), 3, 2))
  constr <- PowConeND(W, z, alpha, axis = 2L)
  nm <- expr_name(constr)
  expect_true(grepl("PowConeND", nm))
})

# ══════════════════════════════════════════════════════════════════
# Cross-cutting
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("All cones inherit from Cone and Constraint", {
  t_var <- Variable(1L)
  x_var <- Variable(c(3L, 1L))
  soc <- SOC(t_var, x_var)
  psd <- PSD(Variable(c(2L, 2L)))
  ec <- ExpCone(Variable(1L), Variable(1L), Variable(1L))
  pc3 <- PowCone3D(Variable(1L), Variable(1L), Variable(1L), 0.5)

  for (c in list(soc, psd, ec, pc3)) {
    expect_true(S7::S7_inherits(c, Cone))
    expect_true(S7::S7_inherits(c, Constraint))
  }
})

## @cvxpy NONE
test_that("constr_size generic dispatches on Constraint base", {
  ## NonPos also gets constr_size via Constraint method
  x <- Variable(c(3L, 1L))
  constr <- NonPos(x)
  expect_equal(constr_size(constr), 3L)
})

## @cvxpy NONE
test_that("is_dgp generic dispatches on Constraint base (FALSE)", {
  x <- Variable(c(3L, 1L))
  constr <- NonPos(x)
  expect_false(is_dgp(constr))
})

## @cvxpy NONE
test_that("save_dual_value generic dispatches on Constraint base", {
  x <- Variable(c(2L, 1L))
  constr <- NonPos(x)
  save_dual_value(constr, matrix(c(1, 2), 2, 1))
  dv <- value(constr@dual_variables[[1L]])
  expect_equal(as.numeric(dv), c(1, 2))
})
