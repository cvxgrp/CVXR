context("test-g01-perspective")

quad_example <- function(request) {
  # Reference expected output.
  x <- Variable()
  s <- Variable()
  
  r <- request$param
  
  obj <- quad_over_lin(x, s) + r*x - 4*s
  constraints <- list(x >= 2, s <= 0.5)
  prob_ref <- Problem(Minimize(obj), constraints)
  result <- solve(prob_ref, solver = "ECOS")
  
  return(list(result$getValue(prob_ref), result$getValue(s), result$getValue(x), r))
}

lse_example <- function() {
  # Reference problem.
  ref_x <- Variable(3)
  ref_s <- Variable()
  ref_z <- Variable(3)
  ref_t <- Variable()
  
  ref_constraints <- list(ref_s >= sum(ref_z), 
                          c(1,2,3) <= ref_x, 1 <= ref_s, ref_s <= 2)
  ref_constraints <- c(ref_constraints, lapply(1:3, function(i) { ExpCone(ref_x[i] - ref_t, ref_s, ref_z[i]) }))
  ref_prob <- Problem(Minimize(ref_t), ref_constraints)
  result <- solve(ref_prob, solver = "ECOS")
  return(list(result$value, result$getValue(ref_x), result$getValue(ref_s)))
}

test_that("test monotonicity", {
  x <- Variable(nonneg = TRUE)
  f <- exp(x)
  s <- Variable(nonneg = TRUE)
  p <- perspective(f, s)
  
  expect_true(is_nonneg(p))
  expect_false(is_nonpos(p))
  
  expect_false(is_incr(p, 1))
  expect_false(is_incr(p, 2))
  
  expect_false(is_decr(p, 1))
  expect_false(is_decr(p, 2))
})

test_that("test exp", {
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- exp(x)
  obj <- perspective(f, s)
  constraints <- list(s >= 1, 1 <= x)
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "ECOS")
  
  # Reference problem.
  ref_x <- Variable()
  ref_s <- Variable()
  ref_z <- Variable()
  
  obj <- ref_z
  ref_constraints <- list(ExpCone(ref_x, ref_s, ref_z), ref_x >= 1, ref_s >= 1)
  ref_prob <- Problem(Minimize(obj), ref_constraints)
  result_ref <- solve(ref_prob, solver = "ECOS")
  
  expect_true(np.isclose(result$value, result_ref$value))
  expect_true(np.isclose(result$getValue(x), result_ref$getValue(ref_x)))
  expect_true(np.isclose(result$getValue(s), result_ref$getValue(ref_s)))
})

test_that("test lse", {
  x <- Variable(3)
  s <- Variable(nonneg = TRUE)
  f <- log_sum_exp(x)
  
  obj <- perspective(f, s)
  constraints <- list(1 <= s, s <= 2, c(1,2,3) <= x)
  prob <- Problem(Minimize(obj), constraints)
  result <- solve(prob, solver = "ECOS")
  
  ref <- lse_example()
  ref_prob <- ref[[1]]
  ref_x <- ref[[2]]
  ref_s <- ref[[3]]
  
  expect_true(np.isclose(result$value, ref_prob))
  expect_true(is.allclose(result$getValue(x), ref_x))
  expect_true(np.isclose(result$getValue(s), ref_s))
})

test_that("test parameter", {
  p <- Parameter(nonneg = TRUE)
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- p*square(x)
  
  obj <- perspective(f, s)
  prob <- Problem(Minimize(obj), list(s <= 1, x >= 2))
  value(p) <- 99
  
  result <- solve(prob)
  
  expect_true(np.isclose(result$value, 4*value(p)))
})

test_that("test affine s", {
  # Test requiring affine s to be nonneg.
  x <- Variable()
  s <- Variable(2)
  expect_error(perspective(square(x), sum(s)), "s must be a variable")
})

test_that("test dpp", {
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  a <- Parameter()
  
  obj <- perspective(square(a + x), s)
  expect_false(is_dpp(obj))
  
  obj <- perspective(log(a + x), s)
  expect_false(is_dpp(obj))
})
