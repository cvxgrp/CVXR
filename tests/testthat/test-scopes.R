## CVXPY 1.9.0 scope-manager parity (utilities/scopes.py).
##
## CVXPY uses `with dpp_scope(): ...` context managers backed by thread-local
## state. CVXR is single-threaded and uses R's on.exit-based equivalents:
##   with_dpp_scope(expr) / dpp_scope_active()
##   with_quad_form_dpp_scope(expr) / quad_form_dpp_scope_active()
## (the quad_form pair is the CVXPY v1.9.0 #3142 addition).
##
## These verify the manager mechanics: the flag is set inside the scope,
## cleared on normal exit, and cleared even when the body raises.

## @cvxpy test_scopes.py::TestScopes::test_dpp_scope
test_that("dpp_scope sets the active flag inside and clears it after", {
  expect_false(dpp_scope_active())
  with_dpp_scope(expect_true(dpp_scope_active()))
  expect_false(dpp_scope_active())
})

## @cvxpy test_scopes.py::TestScopes::test_dpp_scope_active
test_that("dpp_scope_active() reflects the scope state", {
  with_dpp_scope(expect_true(dpp_scope_active()))
  expect_false(dpp_scope_active())
})

## @cvxpy test_scopes.py::TestScopes::test_dpp_scope_exception
test_that("dpp_scope clears the active flag even on error", {
  expect_error(with_dpp_scope(stop("boom")), "boom")
  expect_false(dpp_scope_active())
})

## @cvxpy test_scopes.py::TestScopes::test_quad_form_dpp_scope
test_that("quad_form_dpp_scope sets the active flag inside and clears it after", {
  expect_false(quad_form_dpp_scope_active())
  with_quad_form_dpp_scope(expect_true(quad_form_dpp_scope_active()))
  expect_false(quad_form_dpp_scope_active())
})

## @cvxpy test_scopes.py::TestScopes::test_quad_form_dpp_scope_active
test_that("quad_form_dpp_scope_active() reflects the scope state", {
  with_quad_form_dpp_scope(expect_true(quad_form_dpp_scope_active()))
  expect_false(quad_form_dpp_scope_active())
})

## @cvxpy test_scopes.py::TestScopes::test_quad_form_dpp_scope_exception
test_that("quad_form_dpp_scope clears the active flag even on error", {
  expect_error(with_quad_form_dpp_scope(stop("boom")), "boom")
  expect_false(quad_form_dpp_scope_active())
})

## The two scopes are independent: entering one must not flip the other.
## @cvxpy NONE
test_that("dpp_scope and quad_form_dpp_scope are independent", {
  with_dpp_scope({
    expect_true(dpp_scope_active())
    expect_false(quad_form_dpp_scope_active())
  })
  with_quad_form_dpp_scope({
    expect_true(quad_form_dpp_scope_active())
    expect_false(dpp_scope_active())
  })
})
