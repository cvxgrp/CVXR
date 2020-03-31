context("test-g02-examples")
TOL <- 1e-6

ECOS.dims_to_solver_dict <- CVXR:::ECOS.dims_to_solver_dict

test_that("Find the largest Euclidean ball in the polyhedron", {
    skip_on_cran()
    ## Create the input data
    a1 <- matrix(c(2,1), nrow = 2, ncol = 1)
    a2 <- matrix(c(2,-1), nrow = 2, ncol = 1)
    a3 <- matrix(c(-1,2), nrow = 2, ncol = 1)
    a4 <- matrix(c(-1,-2), nrow = 2, ncol = 1)
    b <- rep(1,4)
    
    ## Create and solve the model
    r <- Variable(name = "r")
    x_c <- Variable(2, name = "x_c")
    obj <- Maximize(r)
    constraints <- list(
        t(a1) %*% x_c + norm(a1,"F")*r <= b[1],
        t(a2) %*% x_c + norm(a2,"F")*r <= b[2],
        t(a3) %*% x_c + norm(a3,"F")*r <= b[3],
        t(a4) %*% x_c + norm(a4,"F")*r <= b[4]
    )
    
    p <- Problem(obj, constraints)
    result <- solve(p)
    expect_equal(result$value, 0.447214, tolerance = TOL)
    expect_equal(result$getValue(r), result$value, tolerance = TOL)
    expect_equal(result$getValue(x_c), matrix(c(0,0)), tolerance = TOL)
})

test_that("Test issue with scalars", {
    skip_on_cran()
    n <- 6
    eps <- 1e-6
    set.seed(10)
    P0 <- matrix(rnorm(n*n), nrow = n, ncol = n)
    eye <- diag(n)
    P0 <- t(P0) %*% P0 + eps*eye
    
    print(P0)
    
    P1 <- matrix(rnorm(n*n), nrow = n, ncol = n)
    P1 <- t(P1) %*% P1
    P2 <- matrix(rnorm(n*n), nrow = n, ncol = n)
    P2 <- t(P2) %*% P2
    P3 <- matrix(rnorm(n*n), nrow = n, ncol = n)
    P3 <- t(P3) %*% P3
    
    q0 <- matrix(rnorm(n), nrow = n)
    q1 <- matrix(rnorm(n), nrow = n)
    q2 <- matrix(rnorm(n), nrow = n)
    q3 <- matrix(rnorm(n), nrow = n)
    
    r0 <- matrix(rnorm(1))
    r1 <- matrix(rnorm(1))
    r2 <- matrix(rnorm(1))
    r3 <- matrix(rnorm(1))
    
    slack <- Variable()
    # Form the problem
    x <- Variable(n)
    objective <- Minimize(0.5*quad_form(x, P0) + t(q0) %*% x + r0 + slack)
    constraints <- list(0.5*quad_form(x, P1) + t(q1) %*% x + r1 <= slack,
                        0.5*quad_form(x, P2) + t(q2) %*% x + r2 <= slack,
                        0.5*quad_form(x, P3) + t(q3) %*% x + r3 <= slack)
    
    # We now find the primal result and compare it to the dual result to check if strong duality holds, i.e., the duality gap is effectively zero
    p <- Problem(objective, constraints)
    result <- solve(p)
    
    # Note that since our data is random, we may need to run this program multiple times to get a feasible primal.
    # When feasible, we can print out the following values.
    print(result$getValue(x))   # Solution
    lam1 <- result$getDualValue(constraints[[1]])
    lam2 <- result$getDualValue(constraints[[2]])
    lam3 <- result$getDualValue(constraints[[3]])
    print(class(lam1))
    
    P_lam <- P0 + lam1*P1 + lam2*P2 + lam3*P3
    q_lam <- q0 + lam1*q1 + lam2*q2 + lam3*q3
    r_lam <- r0 + lam1*r1 + lam2*r2 + lam3*r3
    dual_result <- -0.5*t(q_lam) %*% P_lam %*% q_lam + r_lam
    print(dim(dual_result))
    expect_equal(dim(dual_result), c(1,1))
})

test_that("Test advanced tutorial", {
    ## Solving a problem with different solvers
    skip_on_cran()
    x <- Variable(2)
    obj <- Minimize(x[1] + norm1(x))
    constraints <- list(x >= 2)
    prob <- Problem(obj, constraints)

    ## Solve with ECOS
    result <- solve(prob, solver = "ECOS")
    print(paste("optimal value with ECOS:", result$value))
    expect_equal(result$value, 6, tolerance = TOL)

    ## Solve with ECOS_BB
    result <- solve(prob, solver = "ECOS_BB")
    print(paste("optimal value with ECOS_BB:", result$value))
    expect_equal(result$value, 6, tolerance = TOL)

    ## Solve with CVXOPT
    if("CVXOPT" %in% installed_solvers()) {
        result <- solve(prob, solver = "CVXOPT")
        print(paste("optimal value with CVXOPT:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)
    }

    ## Solve with SCS
    result <- solve(prob, solver = "SCS")
    print(paste("optimal value with SCS:", result$value))
    expect_equal(result$value, 6, tolerance = 1e-2)

    ## Solve with CPLEX
    if("CPLEX" %in% installed_solvers()) {
        result <- solve(prob, solver = "CPLEX")
        print(paste("optimal value with CPLEX:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)
    }

    ## Solve with GLPK
    if("GLPK" %in% installed_solvers()) {
        result <- solve(prob, solver = "GLPK")
        print(paste("optimal value with GLPK:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)

        result <- solve(prob, solver = "GLPK_MI")
        print(paste("optimal value with GLPK_MI:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)
    }

    ## Solve with GUROBI
    if("GUROBI" %in% installed_solvers()) {
        result <- solve(prob, solver = "GUROBI")
        print(paste("optimal value with GUROBI:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)
    }

    print(installed_solvers())
})

test_that("Test log determinant", {
    skip_on_cran()
    x <- t(data.frame(c(0.55, 0.25, -0.2, -0.25, -0.0, 0.4),
                      c(0.0, 0.35, 0.2, -0.1, -0.3, -0.2)))
    n <- nrow(x)
    m <- ncol(x)

    ## Create and solve the model
    A <- Variable(n, n)
    b <- Variable(n)
    obj <- Maximize(log_det(A))
    constraints <- lapply(1:m, function(i) { p_norm(A %*% as.matrix(x[,i]) + b) <= 1 })
    p <- Problem(obj, constraints)
    result <- solve(p)
    expect_equal(result$value, 1.9746, tolerance = 1e-2)
})

test_that("Test portfolio problem", {
    skip_on_cran()
    set.seed(5)
    n <- 100
    m <- 10
    pbar <- matrix(1, nrow = n, ncol = 1) * 0.03 +
        c(matrix(stats::rnorm(n-1), nrow = n-1, ncol = 1), 0) * 0.12

    Fmat <- Matrix::rsparsematrix(m, n, density = 0.01)
    Fmat@x <- rep(1, length(Fmat@x))
    D <- Matrix::sparseMatrix(1:n, 1:n, x = rep(1,n))
    D@x <- stats::rnorm(length(D@x))^2
    Z <- matrix(stats::rnorm(m), nrow = m, ncol = 1)
    Z <- Z %*% t(Z)

    x <- Variable(n)
    y <- Fmat %*% x
    mu <- 1
    ret <- t(pbar) %*% x
    ## DCP attr causes error because not all the curvature
    ## matrices are reduced to constants when an atom is scalar.
    risk <- p_norm(D %*% x)^2 + (Z %*% y)^2
})

test_that("Test examples from CVXR introduction", {
    skip_on_cran()
    m <- 30
    n <- 20
    set.seed(1)
    A <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
    b <- matrix(stats::rnorm(m), nrow = m, ncol = 1)

    ## Construct the problem
    x <- Variable(n)
    objective <- Minimize(sum_squares(A %*% x - b))
    constraints <- list(0 <= x, x <= 1)
    prob <- Problem(objective, constraints)

    ## The optimal objective is returned by solve(prob)
    result <- solve(prob)
    ## The optimal value for x
    print(result$getValue(x))
    ## The optimal Lagrange multiplier for a constraint is stored in result$getDualValue(constraints[[1]])
    print(result$getDualValue(constraints[[1]]))

###########################################
    ## Create two scalar variables
    x <- Variable()
    y <- Variable()

    ## Create two constraints
    constraints <- list(x + y == 1, x - y >= 1)

    ## Form objective
    obj <- Minimize((x - y)^2)

    ## Form and solve problem
    prob <- Problem(obj, constraints)
    result <- solve(prob)
    cat("status:", result$status)
    cat("\noptimal value:", result$value)
    cat("\noptimal var:", result$getValue(x), ", ", result$getValue(y))

    expect_equal(tolower(result$status), "optimal")
    expect_equal(result$value, 1.0, tolerance = TOL)
    expect_equal(result$getValue(x), 1.0, tolerance = TOL)
    expect_equal(result$getValue(y), 0, tolerance = TOL)

###########################################
    ## Replace the objective
    prob@objective <- Maximize(x + y)
    result <- solve(prob)
    cat("optimal value:", result$value)

    expect_equal(result$value, 1.0, tolerance = TOL)

    ## Replace the constraints (x + y == 1)
    prob@constraints[[1]] <- (x + y <= 3)
    result <- solve(prob)
    cat("optimal value:", result$value)

    expect_equal(result$value, 3.0, tolerance = TOL)

###########################################
    x <- Variable()

    ## An infeasible problem
    prob <- Problem(Minimize(x), list(x >= 1, x <= 0))
    # TODO_NARAS_11: result <- solve(prob)
    result <- solve(prob, solver = "ECOS")
    cat("status:", result$status)
    cat("optimal value:", result$value)

    expect_equal(result$status, "infeasible")
    expect_equal(result$value, Inf)

    ## An unbounded problem
    prob <- Problem(Minimize(x))
    result <- solve(prob)
    cat("status:", result$status)
    cat("optimal value:", result$value)

    expect_equal(tolower(result$status), "unbounded")
    expect_equal(result$value, -Inf)

###########################################
    ## A scalar variable
    a <- Variable()

    ## Column vector variable of length 5
    x <- Variable(5)

    ## Matrix variable with 4 rows and 7 columns
    A <- Variable(4, 7)

###########################################
    m <- 10
    n <- 5
    set.seed(1)
    A <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
    b <- matrix(stats::rnorm(m), nrow = m, ncol = 1)

    ## Construct the problem
    x <- Variable(n)
    objective <- Minimize(sum((A %*% x - b)^2))
    constraints <- list(0 <= x, x <= 1)
    prob <- Problem(objective, constraints)

    result <- solve(prob)
    cat("Optimal value:", result$value)
    cat("Optimal var:", result$getValue(x))

    expect_equal(result$value, 7.244277, tolerance = TOL)

###########################################
    ## Positive scalar parameter
    m <- Parameter(nonneg = TRUE)

    ## Column vector parameter with unknown sign (by default)
    c <- Parameter(5)

    ## Matrix parameter with negative entries
    G <- Parameter(4, 7, nonpos = TRUE)

    ## Assigns a constant value to G
    value(G) <- -matrix(1, nrow = 4, ncol = 7)

###########################################
    ## Create parameter, then assign value
    rho <- Parameter(nonneg = TRUE)
    value(rho) <- 2

    ## Initialize parameter with a value
    rho <- Parameter(nonneg = TRUE, value = 2)

###########################################
    ## Problem data
    n <- 15
    m <- 10
    set.seed(1)
    A <- matrix(stats::rnorm(n*m), nrow = n, ncol = m)
    b <- matrix(stats::rnorm(n), nrow = n, ncol = 1)
    ## gamma must be positive due to DCP rules
    gamma <- Parameter(nonneg = TRUE)

    ## Construct the problem
    x <- Variable(m)
    error <- sum_squares(A %*% x - b)

    ## Construct a trade-off curve of ||Ax-b||^2 vs. ||x||_1
    sq_penalty <- c()
    l1_penalty <- c()
    x_values <- list()
    gammas <- 10^seq(-4, 6, length.out = 50)
    for(gamma in gammas) {
        obj <- Minimize(error + gamma*p_norm(x,1))
        prob <- Problem(obj)

        result <- solve(prob)
        sq_penalty <- c(sq_penalty, result$getValue(error))
        l1_penalty <- c(l1_penalty, result$getValue(norm1(x)))
        x_values <- c(x_values, result$getValue(x))
    }

###########################################
    X <- Variable(5, 4)
    A <- matrix(1, nrow = 3, ncol = 5)

    ## Use size(expr) to get the dimensions
    cat("dimensions of X:", size(X))
    cat("\ndimensions of sum_entries(X):", size(sum_entries(X)))
    cat("\ndimensions of A %*% X:", size(A %*% X))

    ## ValueError raised for invalid dimensions
    expect_error(A + X)
})

test_that("Test advanced tutorial 2", {
    skip_on_cran()
    x <- Variable()
    prob <- Problem(Minimize(x^2), list(x == 2))

    # Get ECOS arguments.
    tmp <- get_problem_data(prob, "ECOS")

    # Get ECOS_BB arguments.
    tmp <- get_problem_data(prob, "ECOS_BB")

    # Get CVXOPT arguments.
    if("CVXOPT" %in% installed_solvers())
        tmp <- get_problem_data(prob, "CVXOPT")

    # Get SCS arguments.
    tmp <- get_problem_data(prob, "SCS")

    # Get ECOS arguments.
    library(ECOSolveR)
    tmp <- get_problem_data(prob, "ECOS")
    data <- tmp[[1]]
    chain <- tmp[[2]]
    inverse <- tmp[[3]]

    # Call ECOS solver.
    solution <- ECOSolveR::ECOS_csolve(c = data$c, G = data$G, h = data$h, dims = ECOS.dims_to_solver_dict(data$dims), A = data$A, b = data$b)

    # Unpack raw solver output.
    unpack_results(prob, solution, chain, inverse)
})
