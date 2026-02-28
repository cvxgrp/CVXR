## Tests for DQCP (Disciplined Quasiconvex Programming) curvature queries
## Phase 1: is_quasiconvex, is_quasiconcave, is_quasilinear, is_dqcp

## @cvxpy NONE
test_that("convex expressions are quasiconvex", {
  x <- Variable()
  expect_true(is_quasiconvex(square(x)))
  expect_true(is_quasiconvex(p_norm(x, 2)))
  expect_true(is_quasiconvex(x))  # affine is convex, hence quasiconvex
})

## @cvxpy NONE
test_that("concave expressions are quasiconcave", {
  x <- Variable(nonneg = TRUE)
  expect_true(is_quasiconcave(sqrt(x)))
  expect_true(is_quasiconcave(-square(x)))
  expect_true(is_quasiconcave(x))  # affine is concave, hence quasiconcave
})

## @cvxpy NONE
test_that("affine expressions are quasilinear", {
  x <- Variable(3)
  expect_true(is_quasilinear(x))
  expect_true(is_quasilinear(2 * x + 1))
  expect_true(is_quasiconvex(x))
  expect_true(is_quasiconcave(x))
})

## @cvxpy NONE
test_that("constants are quasiconvex and quasiconcave", {
  c <- Constant(5)
  expect_true(is_quasiconvex(c))
  expect_true(is_quasiconcave(c))
  expect_true(is_quasilinear(c))
  expect_true(is_dqcp(c))
})

## @cvxpy NONE
test_that("Maximum composition: quasiconvex if all args quasiconvex", {
  x <- Variable()
  y <- Variable()
  ## max_elemwise of convex args → quasiconvex
  expr <- max_elemwise(square(x), square(y))
  expect_true(is_quasiconvex(expr))
})

## @cvxpy NONE
test_that("Minimum composition: quasiconcave if all args quasiconcave", {
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  expr <- min_elemwise(sqrt(x), sqrt(y))
  expect_true(is_quasiconcave(expr))
})

## @cvxpy NONE
test_that("MaxEntries composition: quasiconvex if arg quasiconvex", {
  x <- Variable(3)
  expr <- max_entries(square(x))
  expect_true(is_quasiconvex(expr))
})

## @cvxpy NONE
test_that("MinEntries composition: quasiconcave if arg quasiconcave", {
  x <- Variable(3, nonneg = TRUE)
  expr <- min_entries(sqrt(x))
  expect_true(is_quasiconcave(expr))
})

## @cvxpy NONE
test_that("increasing function of quasiconvex arg is quasiconvex (real fn)", {
  x <- Variable()
  ## exp is increasing and scalar → scalar
  ## exp(convex) is convex, hence quasiconvex
  expr <- exp(x)
  expect_true(is_quasiconvex(expr))
})

## @cvxpy NONE
test_that("decreasing function of quasiconcave arg is quasiconvex (real fn)", {
  x <- Variable(nonneg = TRUE)
  ## neg is decreasing, scalar → scalar
  ## neg(concave) → convex → quasiconvex
  expr <- neg(sqrt(x))
  expect_true(is_quasiconvex(expr))
})

## @cvxpy NONE
test_that("is_dqcp on Expression: quasiconvex or quasiconcave", {
  x <- Variable()
  expect_true(is_dqcp(square(x)))   # convex → quasiconvex
  expect_true(is_dqcp(x))           # affine → both
  expect_true(is_dqcp(-square(x)))  # concave → quasiconcave
})

## @cvxpy NONE
test_that("NonPos constraint is_dqcp when arg quasiconvex", {
  x <- Variable()
  con <- NonPos(square(x) - 1)
  expect_true(is_dqcp(con))
  expect_true(is_dcp(con))  # also DCP
})

## @cvxpy NONE
test_that("NonNeg constraint is_dqcp when arg quasiconcave", {
  x <- Variable(nonneg = TRUE)
  con <- NonNeg(sqrt(x) - 1)
  expect_true(is_dqcp(con))
  expect_true(is_dcp(con))  # also DCP
})

## @cvxpy NONE
test_that("Inequality constraint is_dqcp", {
  x <- Variable()
  ## Standard DCP inequality: convex <= constant
  con1 <- Inequality(square(x), Constant(1))
  expect_true(is_dqcp(con1))

  ## Also DCP: quasiconvex <= constant → true
  con2 <- Inequality(x, Constant(1))
  expect_true(is_dqcp(con2))

  ## constant <= quasiconcave (concave is quasiconcave)
  y <- Variable(nonneg = TRUE)
  con3 <- Inequality(Constant(0), sqrt(y))
  expect_true(is_dqcp(con3))
})

## @cvxpy NONE
test_that("Zero and Equality constraints is_dqcp = is_dcp", {
  x <- Variable()
  con_z <- Zero(x - 1)
  expect_true(is_dqcp(con_z))

  con_eq <- Equality(x, Constant(1))
  expect_true(is_dqcp(con_eq))

  ## Non-affine zero constraint: not DCP, not DQCP
  con_z2 <- Zero(square(x) - 1)
  expect_false(is_dqcp(con_z2))
})

## @cvxpy NONE
test_that("Minimize objective is_dqcp when quasiconvex", {
  x <- Variable()
  obj <- Minimize(square(x))
  expect_true(is_dqcp(obj))
})

## @cvxpy NONE
test_that("Maximize objective is_dqcp when quasiconcave", {
  x <- Variable(nonneg = TRUE)
  obj <- Maximize(sqrt(x))
  expect_true(is_dqcp(obj))
})

## @cvxpy NONE
test_that("Problem is_dqcp when all parts are DQCP", {
  x <- Variable()
  prob <- Problem(Minimize(square(x)), list(x >= 1))
  expect_true(is_dqcp(prob))
  expect_true(is_dcp(prob))  # also DCP
})

## @cvxpy NONE
test_that("Non-DCP quasiconvex expressions", {
  x <- Variable()
  ## max_elemwise(x, -x) is quasiconvex (Maximum of quasiconvex args)
  ## but x and -x are affine, so max_elemwise is actually convex
  expr <- max_elemwise(x, -x)
  expect_true(is_quasiconvex(expr))
  expect_true(is_convex(expr))
})

## @cvxpy NONE
test_that("is_atom_quasiconvex defaults to is_atom_convex", {
  x <- Variable()
  sq <- square(x)
  expect_equal(is_atom_quasiconvex(sq), is_atom_convex(sq))
})

## @cvxpy NONE
test_that("is_atom_quasiconcave defaults to is_atom_concave", {
  x <- Variable(nonneg = TRUE)
  sq_root <- sqrt(x)
  expect_equal(is_atom_quasiconcave(sq_root), is_atom_concave(sq_root))
})

## @cvxpy NONE
test_that("cone constraints is_dqcp falls back to is_dcp", {
  x <- Variable(3)
  t_var <- Variable()
  ## SOC constraint: DCP (affine args), so also DQCP
  con <- SOC(t_var, x)
  expect_true(is_dqcp(con))
  expect_equal(is_dqcp(con), is_dcp(con))
})

## @cvxpy NONE
test_that("PSD constraint is_dqcp falls back to is_dcp", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  con <- PSD(X)
  expect_true(is_dqcp(con))
  expect_equal(is_dqcp(con), is_dcp(con))
})
