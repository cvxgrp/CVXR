context("test-g05-conic_solvers-ecos_bb")

test_that("test ecos bb explicit only", {
  x <- Variable(1, name = "x", integer = TRUE)
  objective <- Minimize(sum(x))
  prob <- Problem(objective, list(x >= 0))
  if(!identical(INSTALLED_MI_SOLVERS, list("ECOS_BB"))) {
    result <- solve(prob)
    expect_false(result$solver_stats$solver_name == "ECOS_BB")
  } else
    expect_error(solve(prob), "You need a mixed-integer solver for this model")
})

test_that("test ecos bb lp 0", {
  StandardTestLPs.test_lp_0(solver = "ECOS_BB")
})

test_that("test ecos bb lp 1", {
  # Default settings.
  StandardTestLPs.test_lp_1(solver = "ECOS_BB")
  # Require a basic feasible solution.
  StandardTestLPs.test_lp_1(solver = "ECOS_BB")
})

test_that("test ecos bb lp 2", {
  StandardTestLPs.test_lp_2(solver = "ECOS_BB")
})

test_that("test ecos bb lp 3", {
  StandardTestLPs.test_lp_3(solver = "ECOS_BB")
})

test_that("test ecos bb lp 4", {
  StandardTestLPs.test_lp_4(solver = "ECOS_BB")
})

test_that("test ecos bb lp 5", {
  StandardTestLPs.test_lp_5(solver = "ECOS_BB")
})

test_that("test ecos bb socp 0", {
  StandardTestSOCPs.test_socp_0(solver = "ECOS_BB")
})

test_that("test ecos bb socp 1", {
  StandardTestSOCPs.test_socp_1(solver = "ECOS_BB")
})

test_that("test ecos bb socp 2", {
  StandardTestSOCPs.test_socp_2(solver = "ECOS_BB")
})

test_that("test ecos bb socp 3", {
  # Axis 1.
  StandardTestSOCPs.test_socp_3ax1(solver = "ECOS_BB")
  # Axis 2.
  StandardTestSOCPs.test_socp_3ax2(solver = "ECOS_BB")
})

test_that("test ecos bb expcone 1", {
  StandardTestECPs.test_expcone_1(solver = "ECOS_BB")
})

test_that("test ecos bb exp soc 1", {
  StandardTestMixedCPs.test_exp_soc_1(solver = "ECOS_BB")
})

test_that("test ecos bb mi lp 0", {
  StandardTestLPs.test_mi_lp_0(solver = "ECOS_BB")
})

test_that("test ecos bb mi lp 2", {
  print("Known bug in ECOS BB. Skipping test.")
  return()
  StandardTestLPs.test_mi_lp_2(solver = "ECOS_BB")
})

test_that("test ecos bb mi lp 3", {
  StandardTestLPs.test_mi_lp_3(solver = "ECOS_BB")
})

test_that("test ecos bb mi lp 5", {
  StandardTestLPs.test_mi_lp_5(solver = "ECOS_BB")
})

test_that("test ecos bb mi socp 1", {
  print("Known bug in ECOS BB. Skipping test.")
  return()
  StandardTestSOCPs.test_mi_socp_1(solver = "ECOS_BB")
})
