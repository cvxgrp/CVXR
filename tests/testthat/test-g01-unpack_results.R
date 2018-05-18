TOL <- 1e-6

test_that("Test that unpack_results unpacks things properly for solvers", {
    x <- Variable(2)
    obj <- Minimize(x[1] + cvxr_norm(x, 1))
    constraints <- list(x >= 2)
    prob1 <- Problem(obj, constraints)
    ## Solve with ECOS.
    ecos_data <- get_problem_data(prob1, "ECOS")
    ## Call ECOS solver interface directly
    ecos_output <- ECOSolveR::ECOS_csolve(
                                  c = ecos_data[["c"]],
                                  G = ecos_data[["G"]],
                                  h = ecos_data[["h"]],
                                  dims = ecos_data[["dims"]],
                                  A = ecos_data[["A"]],
                                  b = ecos_data[["b"]]
                              )
    ## Unpack raw solver output.
    res1 <- unpack_results(prob1, "ECOS", ecos_output)
    ## Without DCP validation (so be sure of your math), above is equivalent to:
    ## res1 <- solve(prob1, solver = "ECOS")
    expect_equal(res1$value, 6.0, tolerance = TOL)
    expect_equal(res1$getValue(x), matrix(c(2,2), nrow=2), tolerance = TOL)

    X <- Semidef(2)
    Fmat <- rbind(c(1,0), c(0,-1))
    obj <- Minimize(sum_squares(X - Fmat))
    prob2 <- Problem(obj)
    scs_data <- get_problem_data(prob2, "SCS")
    scs_output <- scs::scs(
                           A = scs_data[['A']],
                           b = scs_data[['b']],
                           obj = scs_data[['c']],
                           cone = scs_data[['dims']]
                       )
    res2 <- unpack_results(prob2, "SCS", scs_output)
    expect_equal(res2$value, 1.0, tolerance = TOL)
    expect_equal(as.numeric(res2$getValue(X)), c(1.0,0,0,0), tolerance = TOL)
})

