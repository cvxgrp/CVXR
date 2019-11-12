context("test-g04-scs")
TOL <- 1e-3

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("Test log problem", {
    ## Log in objective
    skip_on_cran()
    obj <- Maximize(sum(log(x)))
    constr <- list(x <= matrix(c(1, exp(1))))
    p <- Problem(obj, constr)
    result <- solve(p, solver = "SCS")
    expect_equal(result$value, 1, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(1, exp(1))), tolerance = TOL)

    ## Log in constraint
    obj <- Minimize(sum(x))
    constr <- list(log(x) >= 0, x <= matrix(c(1,1)))
    p <- Problem(obj, constr)
    result <- solve(p, solver = "SCS")
    expect_equal(result$value, 2, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(1, 1)), tolerance = TOL)

    ## Index into log
    obj <- Maximize(log(x)[2])
    constr <- list(x <= c(1, exp(1)))
    p <- Problem(obj, constr)
    result <- solve(p, solver = "SCS")
    expect_equal(result$value, 1, tolerance = TOL)
})

test_that("Test sigma_max", {
    skip_on_cran()
    const <- Constant(rbind(c(1,2), c(3,4), c(5,6)))
    constr <- list(C == const)
    prob <- Problem(Minimize(norm(C, "F")), constr)
    result <- solve(prob, solver = "SCS")
    expect_equal(result$value, value(sqrt(sum(const^2))), tolerance = TOL)
    expect_equal(result$getValue(C), value(const), tolerance = TOL)
})

test_that("Test sdp variable", {
    skip_on_cran()
    const <- Constant(rbind(c(1,2,3), c(4,5,6), c(7,8,9)))
    X <- Variable(3, 3, PSD = TRUE)
    prob <- Problem(Minimize(0), list(X == const))
    result <- solve(prob, verbose = TRUE, solver = "SCS")
    expect_equal(tolower(result$status), "infeasible")
})

test_that("Test complex matrices", {
    skip_on_cran()
    # Complex-valued matrix.
    K <- matrix(runif(4), nrow = 2, ncol = 2) + 1i*matrix(runif(4), nrow = , ncol = 2)   # Example matrix.
    n1 <- sum(svd(K)$d)   # Trace norm of K.
    
    # Dual problem.
    X <- Variable(2, 2, complex = TRUE)
    Y <- Variable(2, 2, complex = TRUE)
    # X, Y >= 0 so trace is real.
    objective <- Minimize(Re(0.5*matrix_trace(X) + 0.5*matrix_trace(Y)))
    constraints <- list(bmat(list(list(X, -t(Conj(K))), list(-K, Y))) %>>% 0, X %>>% 0, Y %>>% 0)
    problem <- Problem(objective, constraints)
    
    sol_scs <- solve(problem, solver = "SCS")
    expect_equal(dim(sol_scs$getDualValue(constraints[[1]])), c(4, 4))
    expect_equal(dim(sol_scs$getDualValue(constraints[[2]])), c(2, 2))
    expect_equal(dim(sol_scs$getDualValue(constraints[[3]])), c(2, 2))
    expect_equal(sol_scs$value, n1, tolerance = TOL)
})

test_that("Test a problem with kl_div", {
    skip_on_cran()
    kK <- 50
    set.seed(10)
    
    # Generate a random reference distribution.
    npSPriors <- runif(kK, 0, 1)
    npSPriors <- npSPriors/sum(npSPriors)
    
    # Reference distribution.
    p_refProb <- Parameter(kK, 1, nonneg = TRUE)
    value(p_refProb) <- npSPriors
    
    # Distributio nto be estimated.
    v_prob <- Variable(kK, 1)
    objkl <- 0.0
    for(k in seq_len(kK))
        objkl <- objkl + kl_div(v_prob[k, 1], p_refProb[k, 1])
    
    constrs <- list(do.call("sum", lapply(seq_len(kK), function(k) { v_prob[k, 1] })) == 1)
    klprob <- Problem(Minimize(objkl), constrs)
    result <- solve(klprob, solver = "SCS", verbose = TRUE)
    expect_equal(result$getValue(v_prob), matrix(npSPriors), tolerance = TOL)
})

test_that("Test a problem with entr", {
    skip_on_cran()
    for(n in c(5, 10, 25)) {
        print(n)
        x <- Variable(n)
        obj <- Maximize(sum(entr(x)))
        p <- Problem(obj, list(sum(x) == 1))
        result <- solve(p, solver = "SCS", verbose = TRUE)
        expect_equal(result$getValue(x), matrix(rep(1.0/n, times = n)), tolerance = TOL)
    }
})

test_that("Test a problem with exp", {
    skip_on_cran()
    for(n in c(5, 10, 25)) {
        print(n)
        x <- Variable(n)
        obj <- Minimize(sum(exp(x)))
        p <- Problem(obj, list(sum(x) == 1))
        result <- solve(p, solver = "SCS", verbose = TRUE)
        expect_equal(result$getValue(x), matrix(rep(1.0/n, times = n)), tolerance = TOL)
    }
})

test_that("Test a problem with log", {
    skip_on_cran()
    for(n in c(5, 10, 25)) {
        print(n)
        x <- Variable(n)
        obj <- Maximize(sum(log(x)))
        p <- Problem(obj, list(sum(x) == 1))
        result <- solve(p, solver = "SCS", verbose = TRUE)
        expect_equal(result$getValue(x), matrix(rep(1.0/n, times = n)), tolerance = TOL)
    }
})

test_that("Test warm starting", {
    skip_on_cran()
    x <- Variable(10)
    obj <- Minimize(sum(exp(x)))
    prob <- Problem(obj, list(sum(x) == 1))
    result <- solve(prob, solver = "SCS")
    expect_equal(solve(prob, solver = "SCS")$value, result$value)
    ## expect_false(solve(prob, solver = "SCS", warm_start = TRUE, verbose = TRUE) == result)
})

test_that("Test PSD constraint", {
    skip_on_cran()
    s <- Variable(2, 2)
    obj <- Maximize(min_elemwise(s[1,2], 10))
    const <- list(s %>>% 0, diag(s) == c(1,1))
    prob <- Problem(obj, const)
    r <- solve(prob, solver = "SCS")
    s <- r$getValue(s)
    # TODO_NARAS_13: print(residual(const[[1]]))   TODO: Need to implement the residual function.
    cat("value", r$value)
    cat("s", s)
    eigs <- eigen(s + t(s), only.values = TRUE)$values
    cat("eigs", eigs)
    expect_true(all(eigs >= 0))
})
